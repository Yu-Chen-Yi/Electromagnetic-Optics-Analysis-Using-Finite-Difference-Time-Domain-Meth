%% Simulation problem setting
clear; clc; close all;
c0 = 299792458*1e6;     % Free-space phase velocity (\mum/s)
lambda = 0.632;         % Wavelength (\mum)
STEPS = 2000;           % Numbor of iterations
FDTD_space = 1.5;       % Problem space (\mum)
Slab_THK = 0.25;        % Dummy dielectric slab (\mum)
d_min = Slab_THK;       % Critical dimension
er_slab = 1;            % Relative permittivity of slab
ur_slab = 1;            % Relative permeability of slab
nbc = 1;                % Refraction index at the boundary
%% Initial Grid Resolution :
dz_Wavelength = lambda/50;           % Δ = λ/N; N>=10
dz_Structure = d_min/1;               % Δ = d_min/N; N>=1
dz = min(dz_Wavelength,dz_Structure); % The initial Grid size
%% "Snap" Grid to Critical Dimensions
Mz = ceil(Slab_THK/dz);               
dz = Slab_THK/Mz;                     % Adjust the grid size
%% The Courant Stability Condition
dt = nbc*dz/2/c0*1e15;                % Δt < n*Δ/2/c0 (fs)
%% Gaussian Pulse Source
t = (1:STEPS)*dt;                     % Simulation time (fs)
f_max = c0/lambda*1e-15;              % Maximum frequency
tao = 0.5/f_max;                      % tao = 0.5/f_max;
t0 = tao*4;                           % t0 > 3*tao
Source = exp(-((t-t0)/tao).^2);       % Gaussian Pulse
% plot(t,Source)
%% Build permittivity & permeability matrix of simulation space
Nz = FDTD_space/dz +1;
z = (0:Nz-1)*dz;
slab = z>=0.5*FDTD_space-0.5*Slab_THK & z<=0.5*FDTD_space+0.5*Slab_THK;
er_space = ones(1,Nz)+(er_slab-nbc)*ones(1,Nz).*slab;
ur_space = ones(1,Nz)+(ur_slab-nbc)*ones(1,Nz).*slab;
h1 = figure(1);
set(h1,'Name','Permittivity & Permeability','color','w'...
    ,'units','normalized','outerposition',[0.2 0.6 0.5 0.4])
subplot(2,1,1)
plot(z,er_space,'color','b','linewidth',2)
xlabel("\itz (\mum)"),ylabel("\it\epsilon_r")
set(gca,'Fontname','times new roman')
set(gca,'Fontweight','bold')
set(gca,'fontsize',14)
subplot(2,1,2)
plot(z,ur_space,'color','r','linewidth',2)
xlabel("\itz (\mum)"),ylabel("\it\mu_r")
set(gca,'Fontname','times new roman')
set(gca,'Fontweight','bold')
set(gca,'fontsize',14)
%%
mHx = c0*dt*1e-15./ur_space;
mEy = c0*dt*1e-15./er_space;
Hx = zeros(1,Nz);
Ey = zeros(1,Nz);
H1 = Hx(1);
E1 = Ey(Nz);
E2 = 0;
H2 = 0;
writerObj = VideoWriter('HW1_Video','MPEG-4');
open(writerObj);

for T = 1:1000
    
    H2 = H1; H1 = Hx(1);
    Ey2 = [Ey(2:Nz),E2];
    Ey1 = Ey;

    Hx = Hx + mHx.*(Ey2 - Ey1)/dz;    % Magnetic field Update Equation
    
    E2 = E1; E1 = Ey(Nz);
    Hx1 = Hx;
    Hx2 = [H2,Hx(1:Nz-1)];
    
    Ey = Ey + mEy.*(Hx1 - Hx2)/dz;    % Electric field Update Equation

    % draw FDTD simulation frame
    h2 = figure(2);
    set(h2,'Name','FDTD Monitor','color','w'...
        ,'units','normalized','outerposition',[0.2 0.2 0.5 0.4])
    clf
    rectangle('Position',...          % draw a semi-transparent green rectangle
        [0.5*FDTD_space-0.5*Slab_THK -3 Slab_THK 6],...
        'FaceColor',[0.5 1 0.5 0.5],'EdgeColor','k')
    hold on
    plot(z,real(Ey),'.-','color','b') % draw electric field
    
    plot(z,real(Hx),'.-','color','r') % draw magnetic field
    
%     plot([z(20) z(20)],[-1.5 1.5],... %draw source position
%         'color','k','linewidth',3) 
    hold off
    box on
    xlim([0,FDTD_space]), ylim([-1.5,1.5]), title("Steps = "+T);
    set(gca,'Fontname','times new roman'),set(gca,'Fontweight','bold'),set(gca,'fontsize',14)
%     drawnow
    F = getframe(h2)
    writeVideo(writerObj,F)
end
close(writerObj);