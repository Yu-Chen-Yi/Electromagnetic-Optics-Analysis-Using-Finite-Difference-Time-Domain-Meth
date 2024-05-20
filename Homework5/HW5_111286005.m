%% Simulation problem setting
clear; clc; close all;
c0 = 299792458;                % Free-space phase velocity (m/s)
lambda = 900e-9;               % Wavelength (m)
f_max = c0/lambda;             % Highest frequncy (Hz)
lambda0 = 980e-9;              % Operating wavelength of reflector (m)
n_SiN = 2;                     % Thin film with high refrative index
n_SiO2 = 1.5;                  % Thin film with low refrative index
d_SiN  = lambda0/4/2;          % SiN : quarter wavelength thickness (m)
d_SiO2 = lambda0/4/1.5;        % SiO2 : quarter wavelength thickness (m)
periods = 15;                  % Number of periods
d_min = min([d_SiN,d_SiO2]);   % Critical dimension (m)
er_SiN = n_SiN^2;              % Relative permittivity of SiN2 slab
ur_SiN = 1;                    % Relative permeability of SiN2 slab
er_SiO2 = n_SiO2^2;            % Relative permittivity of SiO2 slab
ur_SiO2 = 1;                   % Relative permeability of SiO2 slab
nbc = 1;                       % Refraction index at the boundary
N_Wavelength = 20;             % Wavelength sampling
N_Structure = 4;               % Feature of structure sampling
nz_src = 2;                    % Source injection point
visualize_rate = 100;          % Time steps per frame
NFREQ = 201;                   % Frequency point
lambda_min = 900*1e-9;         % Spectrum of short wavelength
lambda_max = 1100*1e-9;        % Spectrum of long wavelength
%% Initial Grid Resolution :
n_max = max([n_SiN,n_SiO2,nbc]);      % Maximum of refractive index in space
dz_Wavelength = lambda/N_Wavelength;  % Δ = λ/N; N>=10
dz_Structure = d_min/N_Structure;     % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dz = min(dz_Wavelength,dz_Structure); % The initial Grid size
%% "Snap" Grid to Critical Dimensions
Mz1 = ceil(d_SiN/dz);
Mz2 = ceil(d_SiO2/dz);
Mz = lcm(Mz1,Mz2);
dz = d_SiN/Mz;                       % Adjust the grid size
%% The Courant Stability Condition
dt = nbc*dz/2/c0;                    % Δt < n*Δ/2/c0 (s)
%% Compute the Number of Steps
FDTD_space = 102*dz + periods*(d_SiN+d_SiO2) + 102*dz;  % Problem space (m)
Nz = round(FDTD_space/dz) +1;                           % Number of grids in space
t_prop = n_max*Nz*dz/c0;                                % Time it takes a wave to propagate across the grid
tau = 0.5/f_max;                                        % tau ~= 0.5/f_max;
T = 12*tau + 15*t_prop;                                 % Rule of thumb : T = 12tao + 15t_prop
STEPS = round(T/dt);                                    % Number of iterations
%% Gaussian Pulse Source
t0 = tau*6;                           % t0 > 3*tau
t = (0:STEPS-1)*dt;                   % Simulation time (s)
Source_E = exp(-((t-t0)/tau).^2);     % E field Gaussian Pulse
Source_H = -sqrt(1/1)*exp(-((t-t0+nz_src*dz/2/c0+0.5*dt)/tau).^2);       % H field Gaussian Pulse
h1 = figure(1);
set(h1,'Name','Gaussian Pulse Source','color','w','units','normalized','outerposition',[0 0.2 0.2 0.4])

plot(t*1e9,Source_E,'color','b','linewidth',2),hold on;
plot(t*1e9,Source_E.*sin(2*pi*f_max*t),'color','r','linewidth',2),hold off;axis tight
set(gca,'Fontname','times new roman');set(gca,'fontsize',12)
xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\frac{(t-t_{0})^{2}}{\tau^{2}}}}}$','Interpreter','latex','fontsize',20)

%% Build permittivity & permeability matrix of simulation space
A = ones(1,periods);           %[HL1,HL2,......,HLn]
er_Pair = [er_SiN*ones(1,d_SiN/dz),er_SiO2*ones(1,d_SiO2/dz)];  % permittivity of thin film : (H L)
er_MultiLayer = kron(A,er_Pair);                                % (H L)^n
er_space = [ones(1,101),er_MultiLayer,ones(1,101)];             % spacer region + (H L)^n + spacer region
ur_Pair = [ur_SiN*ones(1,d_SiN/dz),ur_SiO2*ones(1,d_SiO2/dz)];  % permeability of thin films
ur_MultiLayer = kron(A,ur_Pair);                                % (H L)^n
ur_space = [ones(1,101),ur_MultiLayer,ones(1,101)];             % spacer region + (H L)^n + spacer region
z = dz*linspace(0,length(er_space)-1,length(er_space));
h2 = figure(2);
set(h2,'Name','Permittivity & Permeability','color','w','units','normalized','outerposition',[0.25 0.6 0.5 0.4])
subplot(2,1,1)
plot(z,er_space,'color','b','linewidth',2)
xlabel("\itz (\mum)"),ylabel("\it\epsilon_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);axis tight;xticks(0:0.2:2.2);
subplot(2,1,2)
plot(z,ur_space,'color','r','linewidth',2)
xlabel("\itz (\mum)"),ylabel("\it\mu_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);axis tight;xticks(0:0.2:2.2);
%% INITIALIZED FOURIER TRANSFORMS
FREQ = linspace(c0/(lambda_max),c0/(lambda_min),NFREQ);
K = exp(-1i*2*pi*dt.*FREQ);
REF = zeros(1,NFREQ);                       
TRN = zeros(1,NFREQ);
SRC = zeros(1,NFREQ);
%% INITIALIZED FDTD UPDATE PARAMETER
Nz = length(er_space);
mHx = c0*dt./ur_space;
mEy = c0*dt./er_space;
Hx = zeros(1,Nz);
Ey = zeros(1,Nz);
H1 = Hx(1);
E1 = Ey(Nz);
E2 = 0;
H2 = 0;
%%  Visualize a semi-transparent green rectangle
h3 = figure(3);set(h3,'Name','FDTD Monitor','color','w','units','normalized','outerposition',[0.215 0.05 0.78 0.55]);clf;
        subplot(2,3,[1,2])
        plot([z(2) z(2)],[-1.5 1.5],'color','k','linewidth',3)%Visualize source position
for N = 1:periods
        a1=rectangle('Position',[dz*101 + d_SiN*(N-1) + d_SiO2*(N-1)  -1.5 d_SiN 3],...
            'FaceColor',[0.3 0.3 0.3 0.5],'EdgeColor','k');hold on;
        a2=rectangle('Position',[dz*101 + d_SiN*(N) + d_SiO2*(N-1) -1.5 d_SiO2 3],...
            'FaceColor',[0.7 0.7 0.7 0.5],'EdgeColor','k');
end
        aE=plot(z,Ey,'color','b'); % Visualize electric field
        aH=plot(z,Hx,'color','r'); % Visualize magnetic field
        box on;xlabel("\itz (m)");ylabel("\itE_y & H_x");
        ylim([-1.5,1.5]);xlim('tight')
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',12)
        hold off
        
%%
        subplot(2,3,[3 6])
        aS = plot(t*1e9,Source_E,'color','b','linewidth',2);
        xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\left ( \frac{t-t_{0}}{\tau} \right )^{2}}}}$','Interpreter','latex')
        set(gca,'Fontname','times new roman');set(gca,'fontsize',12);hold on
        aP = scatter(1*dt*1e9,Source_E(1),'o','markeredgecolor','k','markerfacecolor','r');hold off;axis tight
        
%% MAIN FDTD LOOP
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
    % UPDATE FOURIER TRANSFORMS
    REF = REF + (K.^T)*Ey(1); 
    TRN = TRN + (K.^T)*Ey(Nz); 
    SRC = SRC + (K.^T)*Source_E(T);  
    % Visualize FDTD simulation frame
    if mod(T,visualize_rate) == 0
        
       
        
        figure(3)
        subplot(2,3,[1,2])
        title("Field at steps "+T+" of "+STEPS);
        set(aE,'ydata',Ey)  % Reset the data.
        set(aH,'ydata',Hx)  % Reset the data.
                
        set(aP,'xdata',T*dt*1e9)
        set(aP,'ydata',Source_E(T))
        
        subplot(2,3,[4 5])
        plot(c0./FREQ*1e9,log10(abs(REF./SRC).^2)*10,'color','r')
        hold on
        plot(c0./FREQ*1e9,log10(abs(TRN./SRC).^2)*10,'color','b')
        plot(c0./FREQ*1e9,log10(abs(REF./SRC).^2+abs(TRN./SRC).^2)*10,'--','color','k')
        xlabel('Wavelength (nm)');ylabel('dB')
        ylim([-40 9])
        legend on
        legend({'Reflectance','Transmittance','T+R'},'Location','northwest','NumColumns',3)
        legend('boxoff')
        hold off
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',12);
        drawnow
    end
end