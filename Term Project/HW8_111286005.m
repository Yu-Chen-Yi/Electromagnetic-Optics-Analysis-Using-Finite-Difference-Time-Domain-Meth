%% Simulation problem setting
clear; clc; close all;
e0 = 8.854*10^-12;      % Permittivity of vacuum [farad/meter]
u0 = 4*pi*10^-7;        % Permeability of vacuum [henry/meter]
c0 = 1/(e0*u0)^0.5;     % Speed of light         [meter/second]
lambda =  1500e-9;      % Wavelength             [meter]
f_max = c0/lambda;      % Highest frequncy       [Hz]
Slab_THK = 1000e-9;     % Dummy dielectric slab  [meter]
er_slab = 1;            % Relative permittivity of slab
ur_slab = 1;            % Relative permeability of slab
d_min = min(Slab_THK);  % Critical dimension
nbc = 1;                % Refraction index at the boundary
N_Wavelength = 20;      % Wavelength sampling
N_Structure = 10;       % Feature of structure sampling
NPML = [0,0,20,20];     % Number of PML cells
Nyscr = 51;             % Source y position
NFREQ = 201;            % Frequency point
lambda_min = 1500*1e-9;  % Spectrum of short wavelength
lambda_max = 1600*1e-9; % Spectrum of long wavelength
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
% PML + space(3λ) + Structure + space(3λ) + PML
Nx = NPML(1) + round(3*lambda/dx) + Slab_THK/dx + round(3*lambda/dx) + NPML(2);
Ny = NPML(3) + round(3*lambda/dy) + Slab_THK/dy + round(3*lambda/dy) + NPML(3);
x = linspace(-(dx*Nx)/2,(dx*Nx)/2,Nx);
y = linspace(-(dy*Ny)/2,(dy*Ny)/2,Ny);
%% 2x grid
Nx2 = 2*Nx;
Ny2 = 2*Ny;
dx2 = dx/2;
dy2 = dy/2;
ER2 = zeros(Nx2,Ny2)+er_slab;
UR2 = zeros(Nx2,Ny2)+ur_slab;

ERzz = ER2(1:2:Nx2,1:2:Ny2);
URxx = UR2(1:2:Nx2,2:2:Ny2);
URyy = UR2(2:2:Nx2,1:2:Ny2);


h1 = figure(1);
set(h1,'Name','Permittivity & Permeability','color','w','units','normalized','outerposition',[0 0.55 0.5 0.4])
subplot(1,2,1)
pcolor(x*1e6,y*1e6,ERzz');axis equal;axis tight;
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");title("\it\epsilon_r");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);
colormap jet;colorbar;shading flat;
hold on
rectangle('Position',[min(x)*1e6 min(y)*1e6 NPML(1)*dx*1e6 Ny*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[max(x)*1e6-NPML(2)*dx*1e6 min(y)*1e6 NPML(2)*dx*1e6 Ny*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[min(x)*1e6 min(y)*1e6 Nx*dx*1e6+NPML(3)*dx NPML(3)*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[min(x)*1e6 max(y)*1e6-NPML(4)*dx*1e6 Nx*dx*1e6 NPML(4)*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
hold off
subplot(1,2,2)
pcolor(x*1e6,y*1e6,URxx');axis equal;axis tight;
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");title("\it\mu_r");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);
colormap jet;colorbar;shading flat;
hold on
rectangle('Position',[min(x)*1e6 min(y)*1e6 NPML(1)*dx*1e6 Ny*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[max(x)*1e6-NPML(2)*dx*1e6 min(y)*1e6 NPML(2)*dx*1e6 Ny*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[min(x)*1e6 min(y)*1e6 Nx*dx*1e6+NPML(3)*dx NPML(3)*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[min(x)*1e6 max(y)*1e6-NPML(4)*dx*1e6 Nx*dx*1e6 NPML(4)*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
hold off
%% Compute the Number of Steps
n_max = max([sqrt(er_slab),nbc]);     % Maximum of refractive index in space
t_prop = n_max*sqrt((Nx*dx)^2+(Ny*dy)^2)/c0;              % Time it takes a wave to propagate across the grid
tau = 1/f_max/pi;                     % tau ~= 0.5/f_max;
T = 12*tau + 1*t_prop;                % Rule of thumb : T = 12tao + 5t_prop
STEPS = round(T/dt);                  % Number of iterations
%% Gaussian Pulse Source
nsrc = 1;                             % material refractive index where source is injectd
t0 = tau*6;                           % t0 > 3*tau
t = (0:STEPS-1)*dt;                   % Simulation time (s)
Source_E = exp(-((t-t0)/tau).^2);     % E field Gaussian Pulse
Source_H = -sqrt(1/1)*exp(-((t-t0+1*dy/2/c0+0.5*dt)/tau).^2);      % H field Gaussian Pulse
h2 = figure(2);
set(h2,'Name','Gaussian Pulse Source','color','w'...
    ,'units','normalized','outerposition',[0.5 0.025 0.2 0.4])
plot(t*1e9,Source_E,'color','b','linewidth',2),hold on;
plot(t*1e9,Source_E.*sin(2*pi*f_max*t),'color','r','linewidth',2),hold off;axis tight
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12)
xlabel("\itt (ns)"),title('$\mathbf{\mathit{e^{-\frac{(t-t_{0})^{2}}{\tau^{2}}}}}$','Interpreter','latex','fontsize',20)
%% Initizlize Fourier transforms
FREQ = linspace(c0/(lambda_max),c0/(lambda_min),NFREQ);
K = exp(-1i*2*pi*dt.*FREQ);
EzR = zeros(Nx,NFREQ); 
EzT = zeros(Nx,NFREQ); 
SRC = zeros(Nx,NFREQ);  
%% PML sigma
sigx = zeros(Nx2,Ny2);
% xmin PML
nx = 2*NPML(1):-1:1;
sigx(1:2*NPML(1),:) = (0.5*e0/dt)*(nx'/2/NPML(1)).^3 .* ones(2*NPML(1),Ny2);
% xmax PML
nx = 1:2*NPML(1);
sigx((Nx2 - 2*NPML(2) +1):Nx2,:) = (0.5*e0/dt)*(nx'/2/NPML(2)).^3 .* ones(2*NPML(1),Ny2);

sigy = zeros(Nx2,Ny2);
% ymin PML
ny = 2*NPML(3):-1:1;
sigy(:,1:2*NPML(3)) = (0.5*e0/dt)*(ny/2/NPML(3)).^3 .* ones(Nx2,2*NPML(3));
% ymax PML
ny = 1:2*NPML(4);
sigy(:,(Ny2 - 2*NPML(4) +1):Ny2) = (0.5*e0/dt)*(ny/2/NPML(4)).^3 .* ones(Nx2,2*NPML(4));

h3 = figure(3);
set(h3,'Name','PML of Sigma','color','w','units','normalized','outerposition',[0.5 0.55 0.5 0.4])
subplot(1,2,1)
imagesc(x*1e6,y*1e6,sigx');xlabel("\itx axis (\mum)");ylabel("\ity axis (\mum)");axis equal;axis tight;
title("\sigma_x");colorbar;colormap jet;
subplot(1,2,2)
imagesc(x*1e6,y*1e6,sigy');xlabel("\itx axis (\mum)");ylabel("\ity axis (\mum)");axis equal;axis tight;
title("\sigma_y");colorbar;colormap jet;
%% Initialized Fields & Curl Arrays & Integration Arrays
 Ez  = zeros(Nx,Ny); Hx  = zeros(Nx,Ny); Hy  = zeros(Nx,Ny);
 CEx = zeros(Nx,Ny); CEy = zeros(Nx,Ny); CHz = zeros(Nx,Ny);
ICEx = zeros(Nx,Ny);ICEy = zeros(Nx,Ny); IDz = zeros(Nx,Ny);

%% Compute Update Coefficients
%Coefficients for Updata equation for Hx
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0  =  (1/dt) + sigHy/(2*e0);
mHx1  = ((1/dt) - sigHy/(2*e0)) ./ mHx0;
mHx2  = - c0./URxx ./ mHx0;
mHx3  = -(c0*dt/e0) * sigHx./URxx ./ mHx0;
%Coefficients for Updata equation for Hy
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0  =  (1/dt) + sigHx/(2*e0);
mHy1  = ((1/dt) - sigHx/(2*e0)) ./ mHy0;
mHy2  = - c0./URyy ./ mHy0;
mHy3  = -(c0*dt/e0) * sigHy./URyy ./ mHy0;
%Coefficients for Updata equation for Dz
sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0  = (1/dt) + (sigDx + sigDy)/(2*e0)...
               + sigDx.*sigDy*(dt/4/e0^2);
mDz1  = (1/dt) - (sigDx + sigDy)/(2*e0)...
               - sigDx.*sigDy*(dt/4/e0^2);
mDz1  = mDz1./mDz0;
mDz2  = c0./mDz0;
mDz4  = - (dt/e0^2)*sigDx.*sigDy./mDz0;
mEz1 = 1./ERzz; Dz = Ez./mEz1;
%% Visualize Fields
h4 = figure(4);
set(h4,'Name','FDTD Monitor','color','w'...
    ,'units','normalized','outerposition',[0 0.025 0.5 0.5]);
subplot(1,2,1)
f1 = imagesc(x*1e6,y*1e6,Ez');
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',18)
colormap jet;
caxis([-1 1]);
colorbar;
axis equal;axis tight
hold on
%PML : xmin、xmax、ymin、ymax
rectangle('Position',[min(x)*1e6 min(y)*1e6 NPML(1)*dx*1e6 Ny*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[max(x)*1e6-NPML(2)*dx*1e6 min(y)*1e6 NPML(2)*dx*1e6 Ny*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[min(x)*1e6 min(y)*1e6 Nx*dx*1e6+NPML(3)*dx NPML(3)*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[min(x)*1e6 max(y)*1e6-NPML(4)*dx*1e6 Nx*dx*1e6 NPML(4)*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
hold off
subplot(1,2,2)
plot((1:STEPS)*dt*1e15,Source_E,'color','b','linewidth',2);hold on;
xlabel("\itt (fs)"),title('$\mathbf{\mathit{e^{-\left ( \frac{t-t_{0}}{\tau} \right )^{2}}}}$','Interpreter','latex','fontsize',36)
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',18)
T = 1;
f2 = plot(T*dt*1e15,Source_E(T),'o','markeredgecolor','k','markerfacecolor','r');

%% Time steps loop
for T = 1:STEPS
    % Compute Curl of E 
    % (BC = Periodic)
    CEx =  diff([Ez,Ez(:,1)],1,2)/dy; 
    CEy = -diff([Ez;Ez(1,:)],1,1)/dx;
    % (BC = Dirichlet)
%     CEx =  diff([Ez,zeros(Nx,1)],1,2)/dy; 
%     CEy = -diff([Ez;zeors(1,Ny)],1,1)/dx;
    CEx(:,Nyscr-1) = diff([Ez(:,Nyscr-1),Ez(:,Nyscr)],1,2)/dy - 1/dy * Source_E(T);
    % Update H integrations
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    % Update H Field
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;
    
%     Hx(:,Nyscr) = Hx(:,Nyscr) + Source_H(T);
%     Hy(:,Nyscr) = Hy(:,Nyscr) + Source_H(T);
    % Complute Curl of H
    % BC : Periodic
    CHz = diff([Hy(end,:);Hy],1,1)/dx - diff([Hx(:,end),Hx],1,2)/dy;
    CHz(:,Nyscr) = diff([Hy(end,Nyscr);Hy(:,Nyscr)],1,1)/dx- diff([Hx(:,Nyscr-1),Hx(:,Nyscr)],1,2)/dy - 1/dy*Source_H(T);
    % BC : Dirichlet
%     CHz = diff([zeros(1,Ny);Hy],1,1)/dx - diff([zeros(Nx,1),Hx],1,2)/dy;
    

    IDz = IDz + Dz;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    % Inject Source
%     Dz(:,Nyscr) = Dz(:,Nyscr) + Source_D(T);
    % Update Ez
    Ez = mEz1.*Dz;
    
    for nfreq = 1 : NFREQ
        EzR(:,nfreq) = EzR(:,nfreq) + (K(nfreq)^T)*Ez(:,30);
        EzT(:,nfreq) = EzR(:,nfreq) + (K(nfreq)^T)*Ez(:,150);
        SRC(:,nfreq) = SRC(:,nfreq) + (K(nfreq)^T)*Ez(:,Nyscr);
    end
    REF = sum(abs(EzR),1)./sum(abs(SRC),1);
    TRN = sum(abs(EzT),1)./sum(abs(SRC),1);
    % Visualize
    subplot(1,2,1)
    set(f1,'CData',Ez')  % Reset the data.
    title("Field at steps "+T+" of "+STEPS);
    set(f2,'XData',T*dt*1e15)
    set(f2,'YData',Source_E(T))
    figure(5)
    plot(FREQ,REF,FREQ,TRN)
    drawnow
end