%% FUNCTION ws3_fitsReading_comet
%  version/date : version 00, 20221013
%  author       : Lorenz Roth - lorenzr@kth.se
%  for KTH course EF2243 - Workshop 3              
%   modified by XXX

%% DESCRIPTION
%
%  INPUT:   DATADIR         Path
%           OBSERVATION     Name fits-file.
% 
%  OUTPUT:  FITSDATASET     FITS-Image + additional Information
%  (header) about dispersion, central wavelength, etc.
%

clear
clc
close all

dir_spice = '/Users/arvwe/Documents/SPICE';
% dir_spice = 'C:/Users/tomgi/Documents/DD KTH/AA Cours/Solar System/WS1';
% dir_spice = '\Users\dansf\OneDrive\Documents\KTH\Solar_System_Physics';
% dir_spice = '/Users/irene/Documents/KTH/year 2/Solar System Physics/';

addpath( fullfile(dir_spice, 'mice/') )
addpath( fullfile(dir_spice, 'mice/lib') )
addpath( fullfile(dir_spice, 'mice/doc/html') )
addpath( fullfile(dir_spice, 'mice/src') )
addpath( fullfile(dir_spice, 'mice/src/mice') )
kernels = ["LSK/naif0010.tls","PCK/pck00010.tpc","SPK/de440s.bsp","SPK/jup230l.bsp"];
for kernel_no = 1:length(kernels)
    cspice_furnsh( fullfile(dir_spice, ['kernels/' char(kernels(kernel_no))]))
end

%  Spectral files
DATADIR     = './spectral_files/';
FITSFILE{1} = 'KeckKc_15jan12_kc1_GW';           % root name of file
FITSFILE{2} = 'KeckBra_15jan12_bra10_GW';
FITSFILE{3} = 'KeckBrac_15jan12_brac13_GW';
FITSFILE{4} = 'KeckPAH_15jan12_pah19_GW';
FITSFILE{5} = 'KeckMs_15jan12_ms7_GW';
FITSFILE{6} = 'KeckLp_15jan12_lp4_GW';
FITSFILE{7} = 'KeckH2O_15jan12_hho16_GW';

cord = zeros(7,2);    
lonlat = zeros(7,2);
Q1 = zeros(7,3);

% Earth and observatory parameters
pos_telescope = [19.826011357437878, -155.47468209999997]; % [lat,lon]
axialTilt_E = 23.4392811; % [deg] To be NPangle in xy_to_latlon-function
platescale = 0.009942; % [arcsec/pixel] Observatory resoultion

for i=1:7           
    FITSDF     = [ DATADIR FITSFILE{i} '.fits' ] ;
     
     % Read FITS-file:
    FITSINFO     = fitsinfo(FITSDF);   % read header information into structure    
    FITSDATASET   = fitsread(FITSDF);  % spectral irradiance per pixel I(x, y) [counts]
    [row,col]       = size(FITSDATASET); % Number of rows and column / image dimension


    %% Q1
    date = "2015-01-12";
    time = char(FITSINFO.PrimaryData.Keywords(15,2));
    et_time = cspice_str2et(char([date + "T" + string(time(1:5))]));
    
    [R_Io,d_angular,lon_subObs,lat_subObs] = epInfo(et_time,pos_telescope,axialTilt_E);
    R_pixels = d_angular/(2*platescale);

    Q1(i,1) = d_angular;
    Q1(i,2) = lon_subObs;
    Q1(i,3) = lat_subObs;
   
    %% Find the center and brightest volcano
    val = max(FITSDATASET, [], 'all');
    [y_volcano,x_volcano] = find(FITSDATASET==val);
    
    j = 1;
    afram = true;
    while afram
        for k=1:150
            if FITSDATASET(j,k) ~= 0
                afram = false;
                break
            end
        end
        j = j+1;
    end
    yup = j;
    xup = k;
    j = 150;
    afram = true;
    while afram
        for k=150:-1:1
            if FITSDATASET(j,k) ~= 0
                afram = false;
                break
            end
        end
        j = j-1;
    end
    ydown = j;
    xdown = k;
    
    x_cent = (150-xup-(150-xdown))/2+xup;
    y_cent = (150-yup-(150-ydown))/2+yup;
    
    dx = x_volcano-x_cent;
    dy = y_volcano-y_cent;
    
    cord(i,1) = dx;
    cord(i,2) = dy;
    
    [lon_Diff,lat_Diff] = xy_to_latlon(cord(i,1),cord(i,2),R_pixels,pos_telescope(1),pos_telescope(2),axialTilt_E);
    lonlat(i,1) = lon_Diff+Q1(i,2);
    lonlat(i,2) = lat_Diff+Q1(i,3);


%% --- Plot raw image ----

    img_plot=FITSDATASET; % multiply with 1e15 to be close to 1

    figure(i);
    
    imagesc(img_plot); 
    hold on
    plot(x_volcano,y_volcano,'o','Color', 'r');
    plot(y_cent,x_cent,'o','Color','m')
      
    set(gca,'ydir','nor') ;  % Define y axis to increase from bottom to top
    axis equal
    xlim([1 col]);
    ylim([1 row]); 
    xlabel('pixel');
    ylabel('pixel'); 
    title(i)

    % display colorbar
    c = colorbar;
    c.Label.String = 'Intensity [GW]';
    c.Label.FontSize = 11;
    clear c
    caxis; % set colormap axis limits
end


     
%% ------------------------- Functions ------------------------------------
function [R_Io,d_angular,lon,lat] = epInfo(et_time,pos_telescope,axialTilt_E) 
    [Io_spos, Io_ltime] = cspice_spkpos('501', et_time, 'ECLIPJ2000', 'LT', 'EARTH');
    d_IoToEarth = sqrt(sum(Io_spos.^2)); % [km] Distance from Earth to Io
    
    spoint = cspice_subsol('NEARPOINT','501', et_time, 'LT', 'EARTH');
    R_Io = sqrt(sum(spoint.^2)); % Radius of Io at latitude of observation
   
    d_angular = 3600*2*atand(R_Io/d_IoToEarth); % [Arcseconds]
    [lon, lat] = xy_to_latlon(spoint(1),spoint(2),R_Io,pos_telescope(1),pos_telescope(2),axialTilt_E);
end



 
