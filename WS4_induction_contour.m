clear all; clc; close all;
global R_call ind_period mu0

OCEAN_OR_CORE = 'o';

% set the parameters
mu0 = 4e-7*pi; % [kgm/s^2A^2]
R_call = 2410*1e3; % [m] Callisto's radius
ind_period=10.2*3600;% [s] synodic period for Callisto
cond_layer = 0.01:0.01:100; % [S/m] Conductivity range for layers

if OCEAN_OR_CORE == 'o'
    % When calculating ocean
    % Crust assumed to be 200 km deep. Ocean between 0 and 1000 km deep
    r_layertop = R_call-200*1e3; % [m] r0, outer radius of conducting sphere 
    r_layerbot = linspace(r_layertop-1000*1e3,r_layertop-100,1e3); % [m] r1, range for the inner radius of the conducting shell (see Figure 1)

elseif OCEAN_OR_CORE == 'c'
    % When calculating core
    % Starting at bottom = 0 (not exactly due to math), and core thickness
    % between 1 and 1000 km
    r_layertop = linspace(600*1e3,1500*1e3,1e3); % [m] outer radius of core
    r_layerbot = 1e-2; % [m] r1, inner radius of the conducting shell. Should be close to 0, but not equal 
end

thickness = r_layertop-r_layerbot; % [m]

%create a grid to calculate the amplitude and phase in dependence of the 
%conductivity and thickness
[cond_layer1,thickness1]= ndgrid(cond_layer,thickness);

%calculate the skin depth
skin=1./sqrt(mu0*pi/ind_period.*cond_layer1); % [m]
d_skin=(r_layertop-r_layerbot)./skin;

%use the function conducting_sphere to calculate the amplitude and phase of
%the induced field
[A,phase]=conducting_sphere(cond_layer1,r_layerbot,r_layertop);


%% plot the isocontours
figure(1)
subplot(1,2,1)
if OCEAN_OR_CORE == 'o'
    contourf(cond_layer1,thickness1./1000,A,(0.1:0.1:0.9),'ShowText','on')
    title('Normalized amplitude: Ocean')
elseif OCEAN_OR_CORE == 'c'
    contourf(cond_layer1,thickness1./1000,A,9,'ShowText','on')
    title('Normalized amplitude: Core')
end
ax1 = gca;
set(ax1,'Ydir','reverse')
set(ax1,'YScale','log')
set(ax1,'XScale','log')
grid on; 
xlabel('Conductivity [S/m]')
ylabel('Thickness [km]')

subplot(1,2,2)
contourf(cond_layer1,thickness1./1000,phase,(-10:-10:-80),'ShowText','on')
ax2 = gca;
set(ax2,'Ydir','reverse')
set(ax2,'YScale','log')
set(ax2,'XScale','log')
grid on; 
xlabel('Conductivity [S/m]')
ylabel('Thickness [km]')

if OCEAN_OR_CORE == 'o'
    title('Phase shift [deg]: Ocean')
elseif OCEAN_OR_CORE == 'c'
    title('Phase shift [deg]: Core')
end

figure(2)
semilogx(cond_layer1,skin(:,1)./1000)
grid on;
xlabel('Conductivity [S/m]')
ylabel('Skin depth [km]')
title('Skin depth')




