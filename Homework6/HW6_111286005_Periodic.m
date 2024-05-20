%% Simulation problem setting
clear; clc; close all;
e0 = 8.854*10^-12;      % Permittivity of vacuum [farad/meter]
u0 = 4*pi*10^-7;        % Permeability of vacuum [henry/meter]
c0 = 1/(e0*u0)^0.5;     % Speed of light         [meter/second]
f_max = 193.41448e12;   % Highest frequncy       [Hz]
lambda = c0/f_max;      % Wavelength             [meter]
Slab_THK = 1000e-9;     % Dummy dielectric slab  [meter]
er_slab = 1;            % Relative permittivity of slab
ur_slab = 1;            % Relative permeability of slab
d_min = min([Slab_THK]);% Critical dimension
nbc = 1;                % Refraction index at the boundary
N_Wavelength = 20;      % Wavelength sampling
N_Structure = 10;       % Feature of structure sampling
%% Initial Grid Resolution :
dx_Wavelength = c0/f_max/N_Wavelength;         % Δ = λ/N; N>=10
dx_Structure = d_min/N_Structure;              % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dx = min(dx_Wavelength,dx_Structure);          % The initial Grid size
%% "Snap" Grid to Critical Dimensions
Mx = ceil(Slab_THK/dx);                        
dx = Slab_THK/Mx;                              % Adjust the grid size
dy = dx;
%% The Courant Stability Condition
dt = nbc/2/c0/sqrt((1/dx)^2 + (1/dy)^2);       % Δt < n/2/c0/sqrt(Δx^-2 + Δy^-2 + Δz^-2) [second]
%% Build permittivity & permeability matrix of simulation space
Nx = round(3*lambda/dx) + Slab_THK/dx + round(3*lambda/dx);
Ny = round(3*lambda/dy) + Slab_THK/dy + round(3*lambda/dy);

er_space = zeros(Nx,Ny)+er_slab;               %er(i,j)
ur_space = zeros(Nx,Ny)+er_slab;               %ur(i,j)
x = 0:dx:dx*(Nx-1);
y = 0:dy:dy*(Ny-1);

h2 = figure(2);
set(h2,'Name','Permittivity & Permeability','color','w','units','normalized','outerposition',[0 0.55 0.5 0.4])
subplot(1,2,1)
imagesc(x*1e6,y*1e6,er_space);axis equal;axis tight;
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");title("\it\epsilon_r");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);
colormap jet;colorbar;
subplot(1,2,2)
imagesc(x*1e6,y*1e6,ur_space);axis equal;axis tight;
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");title("\it\mu_r");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);
colormap jet;colorbar;
%% Compute the Number of Steps
Nx = size(er_space,1);                % Number of grids in space
Ny = size(er_space,2);                % Number of grids in space
n_max = max([sqrt(er_slab),nbc]);     % Maximum of refractive index in space
t_prop = n_max*sqrt((Nx*dx)^2+(Ny*dy)^2)/c0;              % Time it takes a wave to propagate across the grid
tau = 1/f_max/pi;                     % tau ~= 0.5/f_max;
T = 12*tau + 2*t_prop;                % Rule of thumb : T = 12tao + 5t_prop
STEPS = round(T/dt);                  % Number of iterations
%% Gaussian Pulse Source
nsrc = 1;                             % material refractive index where source is injectd
t0 = tau*6;                           % t0 > 3*tau
t = (0:STEPS-1)*dt;                   % Simulation time (s)
Source_E = exp(-((t-t0)/tau).^2);     % E field Gaussian Pulse
Source_D = er_space(round(Nx/3),round(Nx/3)) * sqrt(e0/u0) *  Source_E;     % D field Gaussian Pulse

h1 = figure(1);
set(h1,'Name','Gaussian Pulse Source','color','w'...
    ,'units','normalized','outerposition',[0.5 0.025 0.2 0.4])
plot(t*1e9,Source_D,'color','b','linewidth',2),hold on;
plot(t*1e9,Source_D.*sin(2*pi*f_max*t),'color','r','linewidth',2),hold off;axis tight
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12)
xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\frac{(t-t_{0})^{2}}{\tau^{2}}}}}$','Interpreter','latex','fontsize',20)
%% Initialized FDTD parameters

Ez = zeros(Nx,Ny);
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);

CEx = zeros(Nx,Ny);
CEy = zeros(Nx,Ny);
CHz = zeros(Nx,Ny);
Dz = er_space.*Ez;

h3 = figure(3);set(h3,'Name','FDTD Monitor','color','w'...
    ,'units','normalized','outerposition',[0 0.025 0.5 0.5]);
subplot(1,2,1)
f1 = imagesc(x*1e6,y*1e6,Ez*120*pi);
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',18)
colormap jet; caxis([-0.08 0.08]);axis equal;axis tight
subplot(1,2,2)
plot((1:STEPS)*dt*1e15,Source_D,'color','b','linewidth',2);hold on;
xlabel("\itt (fs)"),title('$\mathbf{\mathit{e^{-\left ( \frac{t-t_{0}}{\tau} \right )^{2}}}}$','Interpreter','latex','fontsize',36)
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',18)
T = 1;
f2 = plot(T*dt*1e15,Source_D(T),'o','markeredgecolor','k','markerfacecolor','r');

for T = 1:STEPS
    %% Magnetic Field Update
    CEx =  diff([Ez,Ez(:,1)],1,2)/dy;
    CEy = -diff([Ez;Ez(1,:)],1,1)/dx;
    
    Hx = Hx + (-c0*dt./ur_space).*CEx;
    Hy = Hy + (-c0*dt./ur_space).*CEy;
    %% Electric displement Field Update
    CHz = diff([Hy(end,:);Hy],1,1)/dx - diff([Hx(:,end),Hx],1,2)/dy;
    
    Dz = Dz + (c0*dt)*CHz;
    Dz(round(Nx/3),round(Ny/3)) = Dz(round(Nx/3),round(Ny/3)) + Source_D(T);
    Ez = Dz./er_space;
    
    %% Visualize FDTD
    subplot(1,2,1)
    set(f1,'CData',Ez*120*pi)  % Reset the data.
    title("Field at steps "+T+" of "+STEPS);
    
    set(f2,'XData',T*dt*1e15)
    set(f2,'YData',Source_D(T))
    drawnow
end