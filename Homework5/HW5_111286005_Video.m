%% Simulation problem setting
clear; clc; close all;
c0 = 299792458;             % Free-space phase velocity (m/s)
lambda = 500e-9;            % Wavelength (m)
f_max = c0/lambda;          % Highest frequncy (Hz)
n_SiN = 2;
n_SiO2 = 1.5;
d_SiN = 980e-9/4/2;         % SiN Dummy dielectric slab (m)
d_SiO2 = 980e-9/4/1.5;         % SiO2 Dummy dielectric slab (m)
periods = 10;
d_min = min([d_SiN,d_SiO2]);   % Critical dimension
er_SiN = n_SiN^2;            % Relative permittivity of slab
ur_SiN = 1;            % Relative permeability of slab
er_SiO2 = n_SiO2^2;            % Relative permittivity of slab
ur_SiO2 = 1;            % Relative permeability of slab
nbc = 1;                % Refraction index at the boundary
%% Initial Grid Resolution :
dz_Wavelength = lambda/20;            % Δ = λ/N; N>=10
dz_Structure = d_min/12;               % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dz = min(dz_Wavelength,dz_Structure); % The initial Grid size
%% "Snap" Grid to Critical Dimensions
Mz1 = ceil(d_SiN/dz);
Mz2 = ceil(d_SiO2/dz);
Mz = lcm(Mz1,Mz2);
dz = d_SiN/Mz;                     % Adjust the grid size
%% The Courant Stability Condition
dt = nbc*dz/2/c0;                     % Δt < n*Δ/2/c0 (s)
%% Compute the Number of Steps
FDTD_space = 102*dz + periods*(d_SiN+d_SiO2) + 102*dz;         % Problem space (m)
Nz = round(FDTD_space/dz) +1;                % Number of grids in space
n_max = max([n_SiN,n_SiO2,nbc]);     % Maximum of refractive index in space
t_prop = n_max*Nz*dz/c0;              % Time it takes a wave to propagate across the grid
tau = 0.5/f_max;                      % tau ~= 0.5/f_max;
T = 12*tau + 15*t_prop;               % Rule of thumb : T = 12tao + 15t_prop
STEPS = round(T/dt);                  % Number of iterations
%% Gaussian Pulse Source
nsrc = 1;
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

%% Build permittivity & permeability matrix of simulation space
A = ones(1,periods);
er_Pair = [er_SiN*ones(1,d_SiN/dz),er_SiO2*ones(1,d_SiO2/dz)];
er_MultiLayer = kron(A,er_Pair);
ur_Pair = [ur_SiN*ones(1,d_SiN/dz),ur_SiO2*ones(1,d_SiO2/dz)];
ur_MultiLayer = kron(A,ur_Pair);
er_space = [ones(1,101),er_MultiLayer,ones(1,101)];
ur_space = [ones(1,101),ur_MultiLayer,ones(1,101)];
z = dz*linspace(0,length(er_space)-1,length(er_space));
h2 = figure(2);
set(h2,'Name','Permittivity & Permeability','color','w','units','normalized','outerposition',[0.25 0.6 0.5 0.4])
subplot(2,1,1)
plot(z,er_space,'color','b','linewidth',2)
xlabel("\itz (\mum)"),ylabel("\it\epsilon_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);axis tight;xticks([0:0.2:2.2]);
subplot(2,1,2)
plot(z,ur_space,'color','r','linewidth',2)
xlabel("\itz (\mum)"),ylabel("\it\mu_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);axis tight;xticks([0:0.2:2.2]);
%% INITIALIZED FOURIER TRANSFORMS
NFREQ = 500;
FREQ = linspace(c0/(1100*1e-9),c0/(500*1e-9),NFREQ);
K = exp(-1i*2*pi*dt.*FREQ);
REF = zeros(1,NFREQ);
TRN = zeros(1,NFREQ);
SRC = zeros(1,NFREQ);
%%
Nz = length(er_space);
mHx = c0*dt./ur_space;
mEy = c0*dt./er_space;
Hx = zeros(1,Nz);
Ey = zeros(1,Nz);
H1 = Hx(1);
E1 = Ey(Nz);
E2 = 0;
H2 = 0;

writerObj = VideoWriter('HW5-2_Video','MPEG-4');
writerObj.FrameRate = 20;
open(writerObj);
nz_src = 2;
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

    REF = REF + (K.^T)*Ey(1); 
    TRN = TRN + (K.^T)*Ey(Nz); 
    SRC = SRC + (K.^T)*Source_E(T);  
    % draw FDTD simulation frame
    if mod(T,60) == 0
        h3 = figure(3);set(h3,'Name','FDTD Monitor','color','w','units','normalized','outerposition',[0.215 0.05 0.78 0.55]);clf;
        subplot(2,3,[1,2])
        % draw a semi-transparent green rectangle
        for N = 1:periods
        rectangle('Position',[dz*101 + d_SiN*(N-1) + d_SiO2*(N-1)  -1 d_SiN 2],...
            'FaceColor',[0.3 0.3 0.3 0.5],'EdgeColor','k');hold on;
        rectangle('Position',[dz*101 + d_SiN*(N) + d_SiO2*(N-1) -1 d_SiO2 2],...
            'FaceColor',[0.7 0.7 0.7 0.5],'EdgeColor','k');
        end
        plot(z,Ey,'color','b') % draw electric field
        plot(z,Hx,'color','r') % draw magnetic field
        plot([z(2) z(2)],[-1.5 1.5],'color','k','linewidth',3)%draw source position
        hold off;box on;xlabel("\itz (m)");ylabel("\itE_y & H_x");title("Field at steps "+T+" of "+STEPS);
        ylim([-1,1]);
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',12)
        
        subplot(2,3,[3 6])
        plot(t*1e9,Source_E,'color','b','linewidth',2)
        xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\left ( \frac{t-t_{0}}{\tau} \right )^{2}}}}$','Interpreter','latex')
        set(gca,'Fontname','times new roman');set(gca,'fontsize',12);hold on
        scatter(T*dt*1e9,Source_E(T),'o','markeredgecolor','k','markerfacecolor','r');hold off;axis tight
        
        subplot(2,3,[4 5])
        plot(c0./FREQ*1e9,abs(REF./SRC).^2,'color','r')
        hold on
        plot(c0./FREQ*1e9,abs(TRN./SRC).^2,'color','b')
        plot(c0./FREQ*1e9,abs(REF./SRC).^2+abs(TRN./SRC).^2,'--','color','k')
        xlabel('Wavelength (nm)')
        ylim([-0.05 1.19])
        legend on
        legend({'Reflectance','Transmittance','T+R'},'Location','northwest','NumColumns',3)
        legend('boxoff')
        hold off
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',12);
%         drawnow
        F = getframe(h3);
        writeVideo(writerObj,F)
    end
end
close(writerObj);