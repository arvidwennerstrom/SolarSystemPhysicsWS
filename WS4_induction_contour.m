clear all; clc; close all;
global R_call r_layertop ind_period mu0

% set the parameters
mu0 = 4e-7*pi; % [kgm/s^2A^2]
R_call = 2410*1e3; % [m] Callisto's radius
ind_period=10.2*3600;% [s] synodic period for Callisto
cond_layer = 0.01:0.1:100; % [S/m] range of the conductivivity of the conducting shell

% When calculating ocean
r_layertop = R_call-140*1e3; % [m] r0, outer radius of conducting sphere 
r_layerbot=linspace(1900*1e3,r_layertop-1e2,10000); % [m] r1, range for the inner radius of the conducting shell (see Figure 1)
 
% When calculating core
% r_layertop = linspace(500*1e3,700*1e3,10); % [m] outer radius of core
% r_layerbot=zeros(1,10); % [m] r1, range for the inner radius of the conducting shell (see Figure 1) 

%create a grid to calculate the amplitude and phase in dependence of the 
%conductivity and thickness
[cond_layer1,r_layerbot1]= ndgrid(cond_layer,r_layerbot);
thickness1 = (r_layertop-r_layerbot1); % [m]

%calculate the skin depth
skin=1./sqrt(mu0*pi/ind_period.*cond_layer1);
d_skin=(r_layertop-r_layerbot1)./skin;

%use the function conducting_sphere to calculate the amplitude and phase of
%the induced field
[A,phase]=conducting_sphere(cond_layer1,r_layerbot1);

%% plot the isocontours
figure(1)
subplot(1,2,1)
contourf(cond_layer1,thickness1./1000,A,(0.1:0.1:0.9),'ShowText','on')
ax1 = gca;
set(ax1,'Ydir','reverse')
set(ax1,'YScale','log')
set(ax1,'XScale','log')
grid on; 
title('Normalized amplitude: Ocean')
xlabel('Conductivity [S/m]')
ylabel('Thickness [km]')

subplot(1,2,2)
contourf(cond_layer1,thickness1./1000,phase,(-10:-10:-80),'ShowText','on')
ax2 = gca;
set(ax2,'Ydir','reverse')
set(ax2,'YScale','log')
set(ax2,'XScale','log')
grid on; 
title('Phase shift [deg]: Ocean')
xlabel('Conductivity [S/m]')
ylabel('Thickness [km]')

figure(2)
contourf(cond_layer1,thickness1./1000,d_skin,'ShowText','on')

