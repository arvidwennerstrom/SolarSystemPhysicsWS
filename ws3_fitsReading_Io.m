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
% dir_spice = '\Users\dansf\OneDrive\Documents\KTH\Solar System Physics';
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

%  Example:
           DATADIR     = './spectral_files/';
           FITSFILE = 'KeckKc_15jan12_kc1_GW';           % root name of file
 
    FITSDF     = [ DATADIR FITSFILE '.fits' ] ;
     
     % Read FITS-file:
     
    FITSINFO     = fitsinfo(FITSDF);   % read header information into structure
        
    FITSDATASET   = fitsread(FITSDF);  % spectral irradiance per pixel I(x, y) [counts]
    
    [row,col]       = size(FITSDATASET); % Number of rows and column / image dimension


%% --- Plot raw image ----
    img_plot=FITSDATASET; % multiply with 1e15 to be close to 1

    figure(1);
    
    imagesc(img_plot); 
      
    set(gca,'ydir','nor') ;  % Define y axis to increase from bottom to top
    axis equal
    xlim([1 col]);
    ylim([1 row]); 
    xlabel('pixel');
    ylabel('pixel'); 

    % display colorbar
    c = colorbar;
    c.Label.String = 'Intensity [GW]';
    c.Label.FontSize = 11;
    clear c
    caxis; % set colormap axis limits

    
    
%% Q1
R_Io = 1821.6; % [km] 
date = "2015-01-12";
time = char(FITSINFO.PrimaryData.Keywords(15,2));
et_time = cspice_str2et(char([date + "T" + string(time(1:5))]));
% et_arr = et_time:1:et_time+40;

[d_angular,spoint] = epInfo(et_time,R_Io);
     
    

function [d_angular,spoint] = epInfo(et_time,R_Io)
    [Io_spos, Io_ltime] = cspice_spkpos('501', et_time, 'ECLIPJ2000', 'LT', 'EARTH');
    d_IoEarth = sqrt(sum(Io_spos.^2));
    d_angular = 2*atand(R_Io/d_IoEarth)*3600; % [Arcseconds]

    spoint = cspice_subsol('NEARPOINT','501', et_time, 'LT', 'EARTH');
    pos_telescope = [19.826011357437878, -155.47468209999997];
    axialTilt_E = 23.4392811;
    [lon, lat] = xy_to_latlon(spoint(1),spoint(2),R_Io,pos_telescope(1),pos_telescope(2),axialTilt_E)


end

 