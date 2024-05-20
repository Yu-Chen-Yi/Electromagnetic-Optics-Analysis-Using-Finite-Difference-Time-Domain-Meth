%% Simulation problem setting
clear; clc; close all;
c0 = 299792458*1e6;     % Free-space phase velocity (\mum/s)
lambda_min = 0.5;       % Minimum Wavelength (\mum)
lambda_max = 0.7;       % Maximum Wavelength (\mum)
STEPS = 1500;           % Numbor of iterations
FDTD_space = 1.65;      % Problem space (\mum)
Slab_THK = 0.25;        % Dummy dielectric slab (\mum)
d_min = Slab_THK;       % Critical dimension
er_slab = 1;            % Relative permittivity of slab
ur_slab = 1;            % Relative permeability of slab
nbc = 1;                % Refraction index at the boundary
%% Initial Grid Resolution :
dz_Wavelength = lambda_min/20;        % Δ = λ/N; N>=10
dz_Structure = d_min/5;               % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dz = min(dz_Wavelength,dz_Structure); % The initial Grid size
%% "Snap" Grid to Critical Dimensions
Mz = ceil(Slab_THK/dz);
dz = Slab_THK/Mz;                     % Adjust the grid size
%% The Courant Stability Condition
dt = nbc*dz/2/c0*1e15;                % Δt < n*Δ/2/]5231]c0 (fs)
%% Gaussian Pulse Source
f = linspace(c0/lambda_max-0.15*1e15,c0/lambda_min+0.15*1e15,1000);
f0   = (c0/lambda_max + c0/lambda_min)/2;% central frequency
fmax = (c0/lambda_min - c0/lambda_max)/2;% frequency content
Gf = exp(-(f-f0).^2/(fmax).^2);          % Frequncy domain Gaussian Pulse

tau = 1/fmax/pi*1e15;                 % tau = 1/f_max/pi;
t0 = tau*6;                           % t0 > 3*tau
t = (0:STEPS-1)*dt;                   % Simulation time (fs)
Source_E = exp(-((t-t0)/tau).^2);     % Time domain Gaussian Pulse
% Source_H = -sqrt(1/1)*exp(-((t-t0+dz/2/c0+0.5*dt)/tau).^2);       % Gaussian Pulse
h1 = figure(1);
set(h1,'Name','Gaussian Pulse Source','color','w','units','normalized','outerposition',[0 0.2 0.2 0.8])

subplot(2,1,1)
plot(f,Gf,'Color','b','linewidth',2);hold on;
line([f0,f0],[0 1],'Color','k','LineStyle','--')
line([c0/lambda_max,c0/lambda_max],[0 1/exp(1)],'Color','k','LineStyle','--')
line([c0/lambda_min,c0/lambda_min],[0 1/exp(1)],'Color','k','LineStyle','--')
line([min(f),max(f)],[1/exp(1) 1/exp(1)],'Color','k','LineStyle','--')
hold off;xlabel("\itfrequency (Hz)"),title('$\mathbf{\mathit{e^{-\frac{(f-f_{0})^{2}}{f_{max}^{2}}}}}$','Interpreter','latex','fontsize',20);axis tight
subplot(2,1,2)
plot(t,Source_E,'color','b','linewidth',2),hold on;
plot(t,Source_E.*sin(2*pi*c0/0.5*1e-15*t),'color','r','linewidth',2),hold off;axis tight
set(gca,'Fontname','times new roman');set(gca,'Fontsize',14)
xlabel("\itt (fs)"),title('$\mathbf{\mathit{e^{-\frac{(t-t_{0})^{2}}{\tau^{2}}}}}$','Interpreter','latex','fontsize',20)

%% Build permittivity & permeability matrix of simulation space
Nz = round(FDTD_space/dz) +1;
z = (1:Nz)*dz;
slab = z>=0.5*FDTD_space-0.5*Slab_THK & z<=0.5*FDTD_space+0.5*Slab_THK;
er_space = ones(1,Nz)+(er_slab-nbc)*ones(1,Nz).*slab;
ur_space = ones(1,Nz)+(ur_slab-nbc)*ones(1,Nz).*slab;
h2 = figure(2);
set(h2,'Name','Permittivity & Permeability','color','w','units','normalized','outerposition',[0.25 0.6 0.5 0.4])
subplot(2,1,1)
plot(z,er_space,'color','b','linewidth',2)
xlabel("\itz (\mum)"),ylabel("\it\epsilon_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',14);axis tight
subplot(2,1,2)
plot(z,ur_space,'color','r','linewidth',2)
xlabel("\itz (\mum)"),ylabel("\it\mu_r")
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',14);axis tight
%%
mHx = c0*dt*1e-15./ur_space;
mEy = c0*dt*1e-15./er_space;
Hx = zeros(1,Nz);
Ey = zeros(1,Nz);
H1 = Hx(1);
E1 = Ey(Nz);
E2 = 0;
H2 = 0;

for T = 1:STEPS
    
    H2 = H1; H1 = Hx(1);
    Ey2 = [Ey(2:Nz),E2];
    Ey1 = Ey;
    
    Hx = Hx + mHx.*(Ey2 - Ey1)/dz;    % Magnetic field Update Equation
    
    E2 = E1; E1 = Ey(Nz);
    Hx1 = Hx;
    Hx2 = [H2,Hx(1:Nz-1)];
    
    Ey = Ey + mEy.*(Hx1 - Hx2)/dz;    % Electric field Update Equation
    
    Ey(10) = Source_E(T) + Ey(10); % Add a E source at position z(20)
    
    % draw FDTD simulation frame
    if mod(T,10) == 0
        h3 = figure(3);set(h3,'Name','FDTD Monitor','color','w','units','normalized','outerposition',[0.215 0.2 0.78 0.4]);clf;
        subplot(1,3,[1,2])
        % draw a semi-transparent green rectangle
        rectangle('Position',[0.5*FDTD_space-0.5*Slab_THK -3 Slab_THK 6],...
            'FaceColor',[0.5 1 0.5 0.5],'EdgeColor','k');hold on;
        plot(z,Ey,'.-','color','b') % draw electric field
        plot(z,Hx,'.-','color','r') % draw magnetic field
        plot([z(10) z(10)],[-1.5 1.5],'color','k','linewidth',3)%draw source position
        hold off;box on;xlabel("\itz (\mum)");ylabel("\itE_y & H_x");title("Steps = "+T);
        xlim([0,FDTD_space]);ylim([-1.5,1.5]);
        set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',14)
        subplot(1,3,3)
        plot(t,Source_E,'color','b','linewidth',2)
        xlabel("\itt (fs)"),title('$\mathbf{\mathit{e^{-\left ( \frac{t-t_{0}}{\tau} \right )^{2}}}}$','Interpreter','latex')
        set(gca,'Fontname','times new roman');set(gca,'fontsize',14);hold on
        scatter(T*dt,Source_E(T),'o','markeredgecolor','k','markerfacecolor','r');hold off
        drawnow
    end
end
