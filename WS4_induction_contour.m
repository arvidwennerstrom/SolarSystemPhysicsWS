clear all; clc; close all;
global R_call r_layertop ind_period mu0

% set the parameters
 mu0 = 4e-7*pi; % [kgm/s^2A^2]
 R_call = 2410*1e3; %Callisto's radius
 

 cond_layer=[0.01 1:100]; % [S/m] range of the conductivivity of the conducting shell (S/m)
 r_layertop=R_call-150*1e3; % [m] r0, outer radius of the conducting shell (see Figure 1)
 r_layerbot=linspace(1800*1e3,r_layertop-1e3,10); % [m] r1, range for the inner radius of the conducting shell (see Figure 1)
 ind_period=10.2*3600;% [s] synodic period for Callisto

  
 %create a grid to calculate the amplitude and phase in dependence of the
 %conductivity and thickness
[cond_layer1,r_layerbot1]= ndgrid(cond_layer,r_layerbot);
thickness = r_layertop-r_layerbot1;

%calculate the skin depth
skin=1./sqrt(mu0*pi/ind_period.*cond_layer1);
d_skin=(r_layertop-r_layerbot1)./skin;

%use the function conducting_sphere to calculate the amplitude and phase of
%the induced field
[A,phase]=conducting_sphere(cond_layer1,r_layerbot1);

%% plot the isocontours
figure(1)
subplot(1,2,1)
contour(log10(cond_layer1),log10(thickness),A,9,'ShowText','on')
xlabel('Conductivity [S/m]')
ylabel('Thickness [m]')
title('Normalized amplitude')
grid on

subplot(1,2,2)
contour(log10(cond_layer1),log10(thickness),phase,9,'ShowText','on')
xlabel('Conductivity [S/m]')
ylabel('Thickness [m]')
title('Phase shift [deg]')
grid on

