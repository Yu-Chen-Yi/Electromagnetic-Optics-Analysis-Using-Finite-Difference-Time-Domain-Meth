%% Simulation problem setting
clear; clc; close all;
e0 = 8.854*10^-12;      % Permittivity of vacuum [farad/meter]
u0 = 4*pi*10^-7;        % Permeability of vacuum [henry/meter]
c0 = 1/(e0*u0)^0.5;     % Speed of light         [meter/second]
lambda =  1500e-9;      % Wavelength             [meter]
f_max = c0/500e-9;     % Highest frequncy       [Hz]
Period = 500e-9;      % Period of grating        [meter]
Fill_Factor = 0.7746;   % Filling factor of grating[a.u.]
phi = 54.746;           % Slanted angle of grating [Degree]
Slab_THK = 691e-9;      % Depth of Blazed grating[meter]
er_sub = 3.476^2;       % Relative permittivity of substrate [a.u.]
ur_sub = 1;             % Relative permeability of substrate [a.u.]
d_min = min(Slab_THK);  % Critical dimension
nbc = 1;                % Refraction index at the boundary
N_Wavelength = 20;      % Wavelength sampling
N_Structure = 7;        % Feature of structure sampling
NPML = [0,0,30,30];     % Number of PML cells
Nyscr = 10 + NPML(3);   % Source y position
NFREQ = 201;            % Frequency point
lambda_min = 1500*1e-9; % Spectrum of short wavelength
lambda_max = 1600*1e-9; % Spectrum of long wavelength
Nanopost_r = 0.2e-6;
% Initial Grid Resolution :
dx_Wavelength = c0/f_max/N_Wavelength;         % Δ = λ/N; N>=10
dx_Structure = d_min/N_Structure;              % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dx = min(dx_Wavelength,dx_Structure);          % The initial Grid size
% "Snap" Grid to Critical Dimensions
Mx = ceil(Slab_THK/dx);                        
dx = Slab_THK/Mx;                              % Adjust the grid size
dy = dx;
% The Courant Stability Condition
dt = nbc/2/c0/sqrt((1/dx)^2 + (1/dy)^2);       % Δt < n/2/c0/sqrt(Δx^-2 + Δy^-2 + Δz^-2) [second]
% Initial Grid Resolution :
grating_min = Slab_THK/N_Structure*cotd(phi);
d_min = min(Slab_THK,grating_min);  % Critical dimension     [meter]
dx_Wavelength = lambda_min/N_Wavelength;       % Δ = λ/N; N>=10
dx_Structure = d_min/N_Structure;              % Δ = d_min/N; N>=1, resolve the smallest feature with at 1 to 4 cells.
dx = min(dx_Wavelength,dx_Structure);          % The initial Grid size
% "Snap" Grid to Critical Dimensions
Mx = ceil(grating_min/dx);                        
dx = grating_min/Mx;                           % Adjust the grid size
dy = dx;
% The Courant Stability Condition
dt = nbc/2/c0/sqrt((1/dx)^2 + (1/dy)^2);       % Δt < n/2/c0/sqrt(Δx^-2 + Δy^-2 + Δz^-2) [second]
% Build permittivity & permeability matrix of simulation space
Dummy_Number = round(0.5*lambda_max/dy); %space(1λ)
Period_Number = round(Period/dx);
Slab_Number = round(Slab_THK/dy);
% 0 + Period + 0
Nx = NPML(1) + Period_Number + NPML(2);
% PML + space(1λ) + Structure + space(1λ) + PML
Ny = NPML(3) + Dummy_Number*2 + Slab_Number + Dummy_Number*2 + NPML(4);
x = linspace(-(dx*Nx)/2,(dx*Nx)/2,Nx);
y = linspace(-(dy*Ny)/2,(dy*Ny)/2,Ny);
% 2x grid
Nx2 = 2*Nx;
Ny2 = 2*Ny;
ER2 = ones(Nx2,Ny2);
UR2 = ones(Nx2,Ny2);
% nanopost
ymin = -1.5e-6;
[~,P1_xidx] = min(abs(x-(0-Nanopost_r)));
[~,P1_yidx] = min(abs(y-ymin));
[~,P2_xidx] = min(abs(x-(0-Nanopost_r)));
[~,P2_yidx] = min(abs(y-(Slab_THK+ymin)));
[~,P3_xidx] = min(abs(x-(0+Nanopost_r)));
[~,P3_yidx] = min(abs(y-(Slab_THK+ymin)));
[~,P4_xidx] = min(abs(x-(0+Nanopost_r)));
[~,P4_yidx] = min(abs(y-ymin));
xidx1 = [P1_xidx P2_xidx P3_xidx P4_xidx P1_xidx]*2;
yidx1 = [P1_yidx P2_yidx P3_yidx P4_yidx P1_yidx]*2;
%
ysub = -1.5e-6;
P0_xidx = 0;
[~,P0_yidx] = min(abs(y-ysub));
[~,P5_xidx] = min(abs(x-Period/2));
[~,P5_yidx] = min(abs(y-ysub));
[~,P6_xidx] = min(abs(x-Period/2));
P6_yidx = 0;
P7_xidx = 0;
P7_yidx = 0;
%
xidx2 = [P0_xidx P5_xidx P6_xidx P7_xidx P0_xidx]*2;
yidx2 = [P0_yidx P5_yidx P6_yidx P7_yidx P0_yidx]*2;
ER2 = ER2 + double(poly2mask(yidx1,xidx1,Nx*2,Ny*2))*(er_sub-1);
UR2 = UR2 + double(poly2mask(yidx1,xidx1,Nx*2,Ny*2))*(ur_sub-1);
ER2 = ER2 + double(poly2mask(yidx2,xidx2,Nx*2,Ny*2))*(er_sub-1);
UR2 = UR2 + double(poly2mask(yidx2,xidx2,Nx*2,Ny*2))*(ur_sub-1);
ERzz = ER2(1:2:Nx2,1:2:Ny2);
URxx = UR2(1:2:Nx2,2:2:Ny2);
URyy = UR2(2:2:Nx2,1:2:Ny2);

h1 = figure(1);
set(h1,'numberTitle','off','Name','Permittivity & Permeability','color','w','units','normalized','outerposition',[0 0.225 0.25 0.6], 'Menu', 'none');
subplot(1,2,1)
imagesc(x*1e6,y*1e6,ERzz');axis equal;axis tight;
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");title("\it\epsilon_r");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);
colormap jet;colorbar;
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
imagesc(x*1e6,y*1e6,URxx');axis equal;axis tight;
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");title("\it\mu_r");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);
colormap jet;colorbar;
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
n_max = max([sqrt(er_sub),nbc]);     % Maximum of refractive index in space
t_prop = n_max*sqrt((Nx*dx)^2+(Ny*dy)^2)/c0;              % Time it takes a wave to propagate across the grid
tau = 1/f_max/pi;                     % tau ~= 0.5/f_max;
T = 12*tau + 6*t_prop;                % Rule of thumb : T = 12tao + 5t_prop
STEPS = round(T/dt);                  % Number of iterations
%% Gaussian Pulse Source
nsrc = 1;                             % material refractive index where source is injectd
t0 = tau*6;                           % t0 > 3*tau
t = (0:STEPS-1)*dt;                   % Simulation time (s)
Source_E = sind(-((t-t0)/tau).^2);     % E field Gaussian Pulse
Source_H = -sqrt(1/1)*sind(-((t-t0+1*dy/2/c0+0.5*dt)/tau).^2);      % H field Gaussian Pulse
%% Initizlize Fourier transforms
FREQ = linspace(c0/(lambda_max),c0/(lambda_min),NFREQ);
K = exp(-1i*2*pi*dt.*FREQ);
KF = kron(K,ones(Nx,1));
EzR = zeros(Nx,NFREQ);
EzT = zeros(Nx,NFREQ); 
EzS = zeros(1,NFREQ); 
Eref = zeros(Nx,NFREQ);
Etrn = zeros(Nx,NFREQ);
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
set(h3,'numberTitle','off','Name','PML of Sigma','color','w','units','normalized','outerposition',[0.24 0.225 0.25 0.6], 'Menu', 'none');
subplot(1,2,1)
imagesc(x*1e6,y*1e6,sigx');xlabel("\itx axis (\mum)");ylabel("\ity axis (\mum)");axis equal;axis tight;
title("\sigma_x");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);colormap jet;colorbar;
subplot(1,2,2)
imagesc(x*1e6,y*1e6,sigy');xlabel("\itx axis (\mum)");ylabel("\ity axis (\mum)");axis equal;axis tight;
title("\sigma_y");set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12);colormap jet;colorbar;
%% Initialized Fields & Curl Arrays & Integration Arrays
 Ez  = zeros(Nx,Ny); Hx  = zeros(Nx,Ny); Hy  = zeros(Nx,Ny);
 CEx = zeros(Nx,Ny); CEy = zeros(Nx,Ny); CHz = zeros(Nx,Ny);
ICEx = zeros(Nx,Ny);ICEy = zeros(Nx,Ny); IDz = zeros(Nx,Ny);

%% Compute Update Coefficients
% Coefficients for Updata equation for Hx
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0  =  (1/dt) + sigHy/(2*e0);
mHx1  = ((1/dt) - sigHy/(2*e0)) ./ mHx0;
mHx2  = - c0./URxx ./ mHx0;
mHx3  = -(c0*dt/e0) * sigHx./URxx ./ mHx0;
% Coefficients for Updata equation for Hy
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0  =  (1/dt) + sigHx/(2*e0);
mHy1  = ((1/dt) - sigHx/(2*e0)) ./ mHy0;
mHy2  = - c0./URyy ./ mHy0;
mHy3  = -(c0*dt/e0) * sigHy./URyy ./ mHy0;
% Coefficients for Updata equation for Dz
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
set(h4,'numberTitle','off','Name','FDTD Monitor','color','w','units','normalized','outerposition',[0.483 0.225 0.115 0.6], 'Menu', 'none');
f1 = imagesc(x*1e6,y*1e6,Ez');
xlabel("\it x axis (\mum)"); ylabel("\it y axis (\mum)");
set(gca,'Fontname','times new roman');set(gca,'Fontweight','bold');set(gca,'fontsize',12)
set(gca,'ydir','normal')
colormap jet;
caxis([-1 1]);
colorbar;
axis equal;axis tight;
hold on
% PML : xmin、xmax、ymin、ymax
rectangle('Position',[min(x)*1e6 min(y)*1e6 NPML(1)*dx*1e6 Ny*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[max(x)*1e6-NPML(2)*dx*1e6 min(y)*1e6 NPML(2)*dx*1e6 Ny*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[min(x)*1e6 min(y)*1e6 Nx*dx*1e6+NPML(3)*dx NPML(3)*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
rectangle('Position',[min(x)*1e6 max(y)*1e6-NPML(4)*dx*1e6 Nx*dx*1e6 NPML(4)*dy*1e6],...
            'FaceColor',[0.1 0.1 0.1 0.5],'EdgeColor','k');
NyTRN = Ny - (NPML(4)+4);
NyREF = NPML(3)+4;
yline(y(Nyscr) * 1e6,'-','color','k','linewidth',2);
yline(y(NyREF) * 1e6,'--','color','k','linewidth',2);
yline(y(NyTRN) * 1e6,'--','color','k','linewidth',2);
hold off
%% Time steps loop
for T = 1:4000
    % Compute Curl of E 
    % (BC = Periodic)
    CEx =  diff([Ez,Ez(:,1)],1,2)/dy; 
    CEy = -diff([Ez;Ez(1,:)],1,1)/dx;
%     % (BC = Dirichlet)
%     CEx =  diff([Ez,zeros(Nx,1)],1,2)/dy; 
%     CEy = -diff([Ez;zeors(1,Ny)],1,1)/dx;
    CEx(:,Nyscr-1) = diff([Ez(:,Nyscr-1),Ez(:,Nyscr)],1,2)/dy - 1/dy * Source_E(T);
    % Update H integrations
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    % Update H Field
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;
    
    % Complute Curl of H
    % BC : Periodic
    CHz = diff([Hy(end,:);Hy],1,1)/dx - diff([Hx(:,end),Hx],1,2)/dy;
    CHz(:,Nyscr) = diff([Hy(end,Nyscr);Hy(:,Nyscr)],1,1)/dx- diff([Hx(:,Nyscr-1),Hx(:,Nyscr)],1,2)/dy - 1/dy*Source_H(T);
    % BC : Dirichlet
%     CHz = diff([zeros(1,Ny);Hy],1,1)/dx - diff([zeros(Nx,1),Hx],1,2)/dy;
    IDz = IDz + Dz;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    Ez = mEz1.*Dz;

    % Visualize
    if mod(T,2)==0
    figure(4)
    f1.CData = Ez';  % Reset the data.
    title({"Field at steps";T}); 
    drawnow
    end
end