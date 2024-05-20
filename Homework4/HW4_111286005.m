%% Simulation problem setting
clear; clc; close all;
c0 = 299792458;         % Free-space phase velocity (m/s)
f_max = 5e9;            % Highest frequncy (Hz)
lambda_0 = c0/2.4e9;    % Wavelength (m)
Slab_THK = 0.3048;      % Dummy dielectric slab (m)
er_slab = 12;           % Relative permittivity of slab
ur_slab = 1;            % Relative permeability of slab
er_ARC  = sqrt(12*1);   % Relative permittivity of anti reflection coating
ur_ARC  = 1;            % Relative permeability of anti reflection coating
ARC_THK = lambda_0/4/sqrt(er_ARC);      % Thickness of anti reflection coating
d_min = min([Slab_THK,ARC_THK]);      % Critical dimension
nbc = 1;                % Refraction index at the boundary
%% Initial Grid Resolution :
dz_Wavelength = c0/f_max/10;          % Δ = λ/N; N>=10
dz_Structure = d_min/10;               % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dz = min(dz_Wavelength,dz_Structure); % The initial Grid size
%% "Snap" Grid to Critical Dimensions
Mz = ceil(ARC_THK/dz);
dz = ARC_THK/Mz;                      % Adjust the grid size
%% The Courant Stability Condition
dt = nbc*dz/2/c0;                     % Δt < n*Δ/2/c0 (s)
%% Build permittivity & permeability matrix of simulation space
% Air + ARC + Slab + ARC + Air
er_space = [ones(1,200),er_ARC*ones(1,ARC_THK/dz),er_slab*ones(1,round(Slab_THK/dz)),er_ARC*ones(1,ARC_THK/dz),ones(1,200)];
ur_space = [ones(1,200),ur_ARC*ones(1,ARC_THK/dz),ur_slab*ones(1,round(Slab_THK/dz)),ur_ARC*ones(1,ARC_THK/dz),ones(1,200)];
z = dz*linspace(0,length(er_space)-1,length(er_space));
h2 = figure(2);
set(h2,'Name','Permittivity & Permeability','color','w','units','normalized','outerposition',[0.25 0.6 0.5 0.4])
subplot(2,1,1)
plot(z,er_space,'color','b','linewidth',2)
xlabel("\itz (m)"),ylabel("\it\epsilon_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);axis tight;grid on;
subplot(2,1,2)
plot(z,ur_space,'color','r','linewidth',2)
xlabel("\itz (m)"),ylabel("\it\mu_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);axis tight;grid on;
%% Compute the Number of Steps
Nz = length(er_space);                % Number of grids in space
n_max = max([sqrt(er_slab),sqrt(er_ARC),nbc]);     % Maximum of refractive index in space
t_prop = n_max*Nz*dz/c0;              % Time it takes a wave to propagate across the grid
tau = 0.5/f_max;                      % tau ~= 0.5/f_max;
T = 12*tau + 5*t_prop;                % Rule of thumb : T = 12tao + 5t_prop
STEPS = round(T/dt);                  % Number of iterations
%% Gaussian Pulse Source
nsrc = 1;                             % material refractive index where source is injectd
tau = 0.5/f_max;                      % tau ~= 0.5/f_max;
t0 = tau*6;                           % t0 > 3*tau
t = (0:STEPS-1)*dt;                   % Simulation time (s)
Source_E = exp(-((t-t0)/tau).^2);     % E field Gaussian Pulse
Source_H = -sqrt(1/1)*exp(-((t-t0+nsrc*dz/2/c0+0.5*dt)/tau).^2);       % H field Gaussian Pulse
h1 = figure(1);
set(h1,'Name','Gaussian Pulse Source','color','w','units','normalized','outerposition',[0 0.2 0.2 0.4])

plot(t*1e9,Source_E,'color','b','linewidth',2),hold on;
plot(t*1e9,Source_E.*sin(2*pi*f_max*t),'color','r','linewidth',2),hold off;axis tight
set(gca,'Fontname','times new roman');set(gca,'fontsize',12)
xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\frac{(t-t_{0})^{2}}{\tau^{2}}}}}$','Interpreter','latex','fontsize',20)
%% Source Parameters
f2 = 5.0 * 1e9;
NFREQ = 1000;
FREQ = linspace(0,f2,NFREQ);
%% Initialize Fourier Transform
K = exp(-1i*2*pi*dt.*FREQ);
REF = zeros(1,NFREQ);
TRN = zeros(1,NFREQ);
SRC = zeros(1,NFREQ);
%% Initialized FDTD parameters
Nz = length(ur_space);
mHx = c0*dt./ur_space;      % Update coefficient
mEy = c0*dt./er_space;      % Update coefficient
Hx = zeros(1,Nz);
Ey = zeros(1,Nz);
H1 = Hx(1);E1 = Ey(Nz);E2 = 0;H2 = 0;

nz_src = 60;                % Source injection point
for T = 1:STEPS
    
    H2 = H1; H1 = Hx(1);              % Perfect boundary conditions
    Ey2 = [Ey(2:Nz),E2];
    Ey1 = Ey;
    
    Hx = Hx + mHx.*(Ey2 - Ey1)/dz;    % Magnetic field Update Equation
    Hx(nz_src-1)=Hx(nz_src-1)-mHx(nz_src-1)*Source_E(T)/dz; % TF/SF source
    
    E2 = E1; E1 = Ey(Nz);             % Perfect boundary conditions
    Hx1 = Hx;
    Hx2 = [H2,Hx(1:Nz-1)];
    Ey = Ey + mEy.*(Hx1 - Hx2)/dz;    % Electric field Update Equation
    Ey(nz_src)=Ey(nz_src)-mEy(nz_src)*Source_H(T)/dz; % TF/SF source
    
    % Update Fourier Transforms
    REF = REF + (K.^T)*Ey(1); 
    TRN = TRN + (K.^T)*Ey(Nz); 
    SRC = SRC + (K.^T)*Source_E(T); 
    
    % Visualize FDTD
    if mod(T,100) == 0
        h3 = figure(3);set(h3,'Name','FDTD Monitor','color','w','units','normalized','outerposition',[0.215 0.05 0.78 0.55]);clf;
        subplot(2,3,[1,2])
        % Draw the slab and anti-reflected coating layer
        rectangle('Position',[200*dz+ARC_THK -1 Slab_THK 2],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');hold on;
        rectangle('Position',[200*dz -1 ARC_THK 2],...
            'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','k');
        rectangle('Position',[200*dz+ARC_THK+Slab_THK -1 ARC_THK 2],...
            'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','k');
        plot(z,Ey,'color','b') % Draw electric field
        plot(z,Hx,'color','r') % Draw magnetic field
        plot([z(nz_src) z(nz_src)],[-1 1],'color','k','linewidth',3)%Draw source position
        hold off;box on;xlabel("\itz (m)");ylabel("\itE_y & H_x");title("Field at steps "+T+" of "+STEPS);
        ylim([-1,1]);
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',12);
        xlim([0 max(z)]);ylim([-1 1])
        subplot(2,3,[3 6])
        plot(t*1e9,Source_E,'color','b','linewidth',2)
        xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\left ( \frac{t-t_{0}}{\tau} \right )^{2}}}}$','Interpreter','latex')
        set(gca,'Fontname','times new roman');set(gca,'fontsize',12);hold on
        scatter(T*dt*1e9,Source_E(T),'o','markeredgecolor','k','markerfacecolor','r');hold off;axis tight
        
        
        subplot(2,3,[4 5])
        plot(FREQ/1e9,abs(REF./SRC).^2,'color','r')
        hold on
        plot(FREQ/1e9,abs(TRN./SRC).^2,'color','b')
        plot(FREQ/1e9,abs(REF./SRC).^2+abs(TRN./SRC).^2,'color','k')
        xlabel('Freq (GHz)')
        ylim([-0.125 1.4])
        legend({'Reflectance','Transmittance','T+R'},'Location','northwest','NumColumns',3)
        hold off
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',12);
        drawnow
    end
end