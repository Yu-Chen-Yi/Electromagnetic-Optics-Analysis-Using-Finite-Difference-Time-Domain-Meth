%% Simulation problem setting
clear; clc; close all;
c0 = 299792458;         % Free-space phase velocity (m/s)
f_max = 1e9;            % Highest frequncy (Hz)
Slab_THK = 0.3;         % Dummy dielectric slab (m)
er_slab = 1;            % Relative permittivity of slab
ur_slab = 1;            % Relative permeability of slab
d_min = Slab_THK;       % Critical dimension
nbc = 1;                % Refraction index at the boundary
%% Initial Grid Resolution :
dz_Wavelength = c0/f_max/20;          % Δ = λ/N; N>=10
dz_Structure = d_min/5;              % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dz = min(dz_Wavelength,dz_Structure); % The initial Grid size
%% "Snap" Grid to Critical Dimensions
Mz = ceil(Slab_THK/dz);
dz = Slab_THK/Mz;                      % Adjust the grid size
%% The Courant Stability Condition
dt = nbc*dz/c0;                     % Δt < n*Δ/2/c0 (s)
%% Build permittivity & permeability matrix of simulation space
% Air + ARC + Slab + ARC + Air
er_space = [ones(1,165),er_slab*ones(1,Slab_THK/dz),ones(1,165)];
ur_space = [ones(1,165),ur_slab*ones(1,Slab_THK/dz),ones(1,165)];
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
n_max = max([sqrt(er_slab),nbc]);     % Maximum of refractive index in space
t_prop = n_max*Nz*dz/c0;              % Time it takes a wave to propagate across the grid
tau = 0.5/f_max;                      % tau ~= 0.5/f_max;
T = 12*tau + 6*t_prop;                % Rule of thumb : T = 12tao + 5t_prop
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
%% Initialized FDTD parameters
Nz = length(ur_space);
mHx = c0*dt./ur_space;      % Update coefficient
mEy = c0*dt./er_space;      % Update coefficient
Hx = zeros(1,Nz);
Ey = zeros(1,Nz);

writerObj = VideoWriter('HW3_Dirichlet_Video','MPEG-4');
writerObj.FrameRate = 20;
open(writerObj);
nz_src = 176;                % Source injection point
for T = 1:STEPS
    
    Ey2 = [Ey(2:Nz),0];
    Ey1 = Ey;
    
    Hx = Hx + mHx.*(Ey2 - Ey1)/dz;    % Magnetic field Update Equation
    Hx(nz_src-1)=Hx(nz_src-1)-mHx(nz_src-1)*Source_E(T)/dz; % TF/SF source
    
    Hx1 = Hx;
    Hx2 = [0,Hx(1:Nz-1)];
    Ey = Ey + mEy.*(Hx1 - Hx2)/dz;    % Electric field Update Equation
    Ey(nz_src)=Ey(nz_src)-mEy(nz_src)*Source_H(T)/dz; % TF/SF source
    
    % Visualize FDTD
    if mod(T,5) == 0
        h3 = figure(3);set(h3,'Name','FDTD Monitor','color','w','units','normalized','outerposition',[0.215 0.05 0.78 0.55]);clf;
        subplot(1,3,[1,2])
        % Draw the slab layer
        rectangle('Position',[164*dz -1 Slab_THK 2],...
            'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','k');hold on;
        plot(z,Ey,'color','b') % Draw electric field
        plot(z,Hx,'color','r') % Draw magnetic field
        plot([z(nz_src) z(nz_src)],[-1 1],'color','k','linewidth',3)%Draw source position
        hold off;box on;xlabel("\itz (m)");ylabel("\itE_y & H_x");title("Field at steps "+T+" of "+STEPS);
        ylim([-1,1]);
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',12);
        xlim([0 max(z)]);ylim([-1 1])
        subplot(1,3,3)
        plot(t*1e9,Source_E,'color','b','linewidth',2)
        xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\left ( \frac{t-t_{0}}{\tau} \right )^{2}}}}$','Interpreter','latex')
        set(gca,'Fontname','times new roman');set(gca,'fontsize',12);hold on
        scatter(T*dt*1e9,Source_E(T),'o','markeredgecolor','k','markerfacecolor','r');hold off;axis tight
%         drawnow
        F = getframe(h3);
        writeVideo(writerObj,F)
    end
end
close(writerObj);