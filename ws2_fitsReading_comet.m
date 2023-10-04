%% FUNCTION ws2_fitsReading_comet
%  version/date : version 02, 220926
%  author       : Lorenz Roth - lorenzr@kth.se
%  for KTH course EF2243 - Workshop 2              
%   modified by Arvid Wennerström, Irene Nocerino, Ragnheiður
%   Tryggvadóttir, Thomas Giudicissi

%% DESCRIPTION
%
%  INPUT:   DATADIR         Path
%           OBSERVATION     Name fits-file.
% 
%  OUTPUT:  FITSDATASET     FITS-Image + additional Information
%  (header) about dispersion, central wavelength, etc.
%
clear; clc; close all

%  Example:
           DATADIR     = './HST_Data_46P/';
           OBSERVATION = 'odx605010';           % root name of file

%% main funtion
    dir_spice = '/Users/arvwe/Documents/SPICE';
    % dir_spice = 'C:/Users/tomgi/Documents/DD KTH/AA Cours/Solar System/WS1';
    % dir_spice = '\Users\dansf\OneDrive\Documents\KTH\Solar System Physics';
    % dir_spice = '/Users/irene/Documents/KTH/year 2/Solar System Physics/';


    FITSFILE    = [ OBSERVATION '_flt'];
    FITSDIR     = [ DATADIR FITSFILE '.fits' ] ;
     
     % Read FITS-file:
     
    FITSINFO     = fitsinfo(FITSDIR);   % read header information into structure
        
    FITSDATASET.RAWCOUNTS   = fitsread(FITSDIR,'image',1);  % spectral irradiance per pixel I(x, y) [counts]
    FITSDATASET.ERROR       = fitsread(FITSDIR,'image',2);  % statistical error per pixel sigma(y, x) [counts]
    
    [row,col]       = size(FITSDATASET.RAWCOUNTS); % Number of rows and column / image dimensions
           
    % -- Read specific parameters from header information --
     
    % - EITHER by directly pointing to the sub-indices of the FITSINFO structure

    dateobs = FITSINFO.PrimaryData.Keywords{37,2};  % date start exposure
    timeobs = FITSINFO.PrimaryData.Keywords{38,2};  % time start exposure
    
    % - OR by using the "KeyFinder" routine to search for your keyword by name

    % Read general parameters 
    platesc   = KeyFinder(FITSINFO,'PLATESC'); % plate scale of a detector pixel (arcseconds)
    cenwave  = KeyFinder(FITSINFO,'CENWAVE'); % wavelength (Å) at reference x pixel
    
    % Read parameters specific for image 
    exptime   = KeyFinder(FITSINFO,'EXPTIME','image');  % exposure time (s)
    crpx1    = KeyFinder(FITSINFO,'CRPIX1','image'); % reference x pixel index of target location
    crpx2    = KeyFinder(FITSINFO,'CRPIX2','image'); % reference y pixel index of target location
    cd11     = KeyFinder(FITSINFO,'CD1_1','image');  % dispersion Delta lambda (Å / pixel)
         

    % -- Define image axes 
    xaxis_pix        = ([1:col]) ; % x axis as pixel (spatial) axis
    xaxis_lambda = cenwave+([1:col]-crpx1)*cd11; % x axis as spectral axis
    
    yaxis_pix        = ([1:row]) ; % spatial y axis 
           
    xaxis_arcsec = xaxis_pix/platesc;   % axes as spatial (angular) axes in arcseconds
    yaxis_arcsec = yaxis_pix/platesc;

%% --- Plot raw image ----

    img_plot=FITSDATASET.RAWCOUNTS*1e15; % multiply with 1e15 to be close to 1
    fig_caxis = [0 3];     % Set arbitrary bounds to displayed image
    
    figure(1);
    
    imagesc(img_plot, fig_caxis); 
      
    set(gca,'ydir','nor') ;  % Define y axis to increase from bottom to top
    axis equal
    xlim([1 col]);
    ylim([1 row]); 
    xlabel('pixel');
    ylabel('pixel'); 

    % display colorbar
    c = colorbar;
    c.Label.String = 'Intensity [counts]';
    c.Label.FontSize = 11;
    clear c
    caxis; % set colormap axis limits
     
    
 %% --- Plot spectrum by summing along slit --- 
    figure(2);
 
    for jj=1:col
        total_along_the_slit(jj)  = sum(FITSDATASET.RAWCOUNTS(:,jj)) ;
        error_total_along_the_slit(jj)  = sqrt(sum( (FITSDATASET.ERROR(:,jj)).^2 )) ;
    end

    plot(xaxis_lambda,total_along_the_slit); % spectrum (flux vs column)
       
        ylabel('Intensity (counts)')  ;
        ylim([min(total_along_the_slit) max(total_along_the_slit)]);
        xlim([1700 3150]);
        tit=sprintf('Summed counts along the slit');
        title(tit,'Interpreter','none');
        grid on
   
    a1Pos = get(gca,'Position');

    hold on
    plot(xaxis_lambda,error_total_along_the_slit); % propagated error 
    hold off

    %// Place axis 2 below the 1st.
    ax2 = axes('Position',[a1Pos(1) a1Pos(2)-.05 a1Pos(3) a1Pos(4)],'Color','none','YTick',[],'YTickLabel',[]);
    xlim([min(xaxis_pix(:)) max(xaxis_pix(:))])  

    xlabel('Wavelength (Å)  -   Pixel') ;

 %% --
 %-- Helpful factors for unit conversion from [counts] to photon flux

    platesc_rad = deg2rad(platesc / 3600);   % pixel platescale in radians
    omega_pix = platesc_rad^2;       % solid angle of one pixels in steradians (rad^2)

 %% Q1
 % Offset from comet nucleus towards the Sun, in km
 offset_km = offsetDistance(dir_spice);


 %% -------------------------- Functions -------------------------------
function [dist_offset] = offsetDistance(dir_spice) 
    addpath( fullfile(dir_spice, 'mice/') )
    addpath( fullfile(dir_spice, 'mice/lib') )
    addpath( fullfile(dir_spice, 'mice/doc/html') )
    addpath( fullfile(dir_spice, 'mice/src') )
    addpath( fullfile(dir_spice, 'mice/src/mice') )
    
    kernels = ["LSK/naif0010.tls","SPK/de440s.bsp","SPK/hst.bsp","SPK/1000109.bsp"];
    for kernel_no = 1:length(kernels)
        cspice_furnsh( fullfile(dir_spice, ['kernels/' char(kernels(kernel_no))]))
    end
    
    dist_offset = []; % [km]
    angle_offset = [0,2.5,8]./3600; % [deg] Offset towards the Sun
    measurements = ["2019-01-11T09:51","2019-01-11T10:23";"2019-01-11T14:34","2019-01-11T15:06";"2019-01-14T12:32","2019-01-14T13:03"];
    for measure_no = 1:length(measurements)
        et_interval = cspice_str2et([char(measurements(measure_no,1));char(measurements(measure_no,2))]);
        et_step  = 1;
        et_arr = et_interval(1):et_step:et_interval(2);
        [hub_spos, hub_ltime] = cspice_spkpos('-48', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
        [wir_spos, wir_ltime] = cspice_spkpos('1000109', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
        dist = mean(cspice_vdist(hub_spos,wir_spos));
        dist_offset = [dist_offset dist*tand(angle_offset(measure_no))];
    end
end





