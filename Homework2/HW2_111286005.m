%% Simulation problem setting
clear; clc; close all;
c0 = 299792458;         % Free-space phase velocity (m/s)
f_max = 1e9;            % Highest frequncy (Hz)
lambda = c0/f_max;      % Wavelength (m)
FDTD_space = 5;       % Problem space (m)
Slab_THK = 0.3;         % Dummy dielectric slab (m)
d_min = Slab_THK;       % Critical dimension
er_slab = 1;            % Relative permittivity of slab
ur_slab = 1;            % Relative permeability of slab
nbc = 1;                % Refraction index at the boundary
%% Initial Grid Resolution :
dz_Wavelength = lambda/20;            % Δ = λ/N; N>=10
dz_Structure = d_min/5;               % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dz = min(dz_Wavelength,dz_Structure); % The initial Grid size
%% "Snap" Grid to Critical Dimensions
Mz = ceil(Slab_THK/dz);
dz = Slab_THK/Mz;                     % Adjust the grid size
%% The Courant Stability Condition
dt = nbc*dz/2/c0;                     % Δt < n*Δ/2/c0 (s)
%% Compute the Number of Steps
Nz = FDTD_space/dz +1;                % Number of grids in space
n_max = max([sqrt(er_slab),nbc]);     % Maximum of refractive index in space
t_prop = n_max*Nz*dz/c0;              % Time it takes a wave to propagate across the grid
tau = 0.5/f_max;                      % tau ~= 0.5/f_max;
T = 12*tau + 0.5*t_prop;                % Rule of thumb : T = 12tao + 5t_prop
STEPS = round(T/dt);                  % Number of iterations
%% Gaussian Pulse Source
tau = 0.5/f_max;                      % tau ~= 0.5/f_max;
t0 = tau*6;                           % t0 > 3*tau
t = (0:STEPS-1)*dt;                   % Simulation time (s)
Source_E = exp(-((t-t0)/tau).^2);     % E field Gaussian Pulse
h1 = figure(1);
set(h1,'Name','Gaussian Pulse Source','color','w','units','normalized','outerposition',[0 0.2 0.2 0.4])

plot(t*1e9,Source_E,'color','b','linewidth',2),hold on;
plot(t*1e9,Source_E.*sin(2*pi*f_max*t),'color','r','linewidth',2),hold off;axis tight
set(gca,'Fontname','times new roman');set(gca,'Fontsize',14)
xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\frac{(t-t_{0})^{2}}{\tau^{2}}}}}$','Interpreter','latex','fontsize',20)

%% Build permittivity & permeability matrix of simulation space
z = (1:Nz)*dz;
slab = z>=0.5*FDTD_space-0.5*Slab_THK & z<=0.5*FDTD_space+0.5*Slab_THK;
er_space = ones(1,Nz)+(er_slab-nbc)*ones(1,Nz).*slab;
ur_space = ones(1,Nz)+(ur_slab-nbc)*ones(1,Nz).*slab;
h2 = figure(2);
set(h2,'Name','Permittivity & Permeability','color','w','units','normalized','outerposition',[0.25 0.6 0.5 0.4])
subplot(2,1,1)
plot(z,er_space,'color','b','linewidth',2)
xlabel("\itz (m)"),ylabel("\it\epsilon_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',14);axis tight;xticks([0:0.5:5]);
subplot(2,1,2)
plot(z,ur_space,'color','r','linewidth',2)
xlabel("\itz (m)"),ylabel("\it\mu_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',14);axis tight;xticks([0:0.5:5]);
%%
mHx = c0*dt./ur_space;
mEy = c0*dt./er_space;
Hx = zeros(1,Nz);
Ey = zeros(1,Nz);
H1 = Hx(1);
E1 = Ey(Nz);
E2 = 0;
H2 = 0;

for T = 1:STEPS
    
    H2 = H1; H1 = Hx(1);              % Perfect boundary conditions
    Ey2 = [Ey(2:Nz),E2];
    Ey1 = Ey;
    
    Hx = Hx + mHx.*(Ey2 - Ey1)/dz;    % Magnetic field Update Equation
    
    E2 = E1; E1 = Ey(Nz);             % Perfect boundary conditions
    Hx1 = Hx;
    Hx2 = [H2,Hx(1:Nz-1)];
    
    Ey = Ey + mEy.*(Hx1 - Hx2)/dz;    % Electric field Update Equation
    
    Ey(175) = Source_E(T) + Ey(175); % Add a E source at position z(25)

    % draw FDTD simulation frame
    if mod(T,10) == 0
        h3 = figure(3);set(h3,'Name','FDTD Monitor','color','w','units','normalized','outerposition',[0.215 0.05 0.78 0.55]);clf;
        subplot(1,3,[1,2])
        % draw a semi-transparent green rectangle
        rectangle('Position',[0.5*FDTD_space-0.5*Slab_THK -3 Slab_THK 6],...
            'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','k');hold on;
        plot(z,Ey,'.-','color','b') % draw electric field
        plot(z,Hx,'.-','color','r') % draw magnetic field
        plot([z(175) z(175)],[-1.5 1.5],'color','k','linewidth',3)%draw source position
        hold off;box on;xlabel("\itz (m)");ylabel("\itE_y & H_x");title("Field at steps "+T+" of "+STEPS);xticks([0:0.5:5])
        xlim([0,FDTD_space]);ylim([-1,1]);
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',14)
        subplot(1,3,3)
        plot(t*1e9,Source_E,'color','b','linewidth',2)
        xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\left ( \frac{t-t_{0}}{\tau} \right )^{2}}}}$','Interpreter','latex')
        set(gca,'Fontname','times new roman');set(gca,'fontsize',14);hold on
        scatter(T*dt*1e9,Source_E(T),'o','markeredgecolor','k','markerfacecolor','r');hold off;axis tight
        drawnow
    end
end