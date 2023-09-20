
%% ------------------------------------------------------------------------
%
% Lorenz Roth / 1 Sep 2022 
% Example program for WS1 on SPICE / mice for EF2243 HT2022
%
%% ------------------------------------------------------------------------

% Your directory where SPICE toolkit and kernels are located
% ! Change to the correct directory on your computer ! 
dir_spice = '/Users/arvwe/Documents/SolarSystemPhysicsWS1';
% dir_spice = 'C:/Users/tomgi/Documents/DD KTH/AA Cours/Solar System/WS1';
% dir_spice = '\Users\dansf\OneDrive\Documents\KTH\Solar System Physics';
% dir_spice = '/Users/irene/Documents/KTH/year 2/Solar System Physics/';

% Add paths to routines of toolkit so that MATLAB finds routines
% ! more paths need to be linked than mentioend in description !
addpath( fullfile(dir_spice, 'mice/') )
addpath( fullfile(dir_spice, 'mice/lib') )
addpath( fullfile(dir_spice, 'mice/doc/html') )
addpath( fullfile(dir_spice, 'mice/src') )
addpath( fullfile(dir_spice, 'mice/src/mice') )

%% -- Load kernels that are needed 
kernels = ["LSK/naif0010.tls","SPK/de440s.bsp","SPK/mar097.bsp","SPK/jup230l.bsp","SPK/sat415.bsp","SPK/ura111.bsp","SPK/nep095.bsp","SPK/plu058.bsp","PCK/pck00010.tpc","SPK/Voyager_2.m05016u.merged.bsp"];
for kernel_no = 1:length(kernels)
    cspice_furnsh( fullfile(dir_spice, ['kernels/' char(kernels(kernel_no))]))
end

% Define a start and end date for a time interval you want to use 
% and convert the date to ephemeris time (Unit: seconds)
et_mission  = cspice_str2et(['1977-08-20T15:45';'2023-09-30T12:00']); % From start of the observation.
et_step  = 4*3600 ; % time step of 1 hour in ephemeris (in seconds)
et_arr   = et_mission(1):et_step:et_mission(2); % Create a time array from start to end

cspice_str2et('2000-01-01T11:58:56');

% Calculate the positions of planets
% - uses Ecliptic Coordinate System, with third coord pointing to Ecliptic North
[mar_spos, mar_ltime] = cspice_spkpos('499', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
[jup_spos, jup_ltime] = cspice_spkpos('599', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
[sat_spos, sat_ltime] = cspice_spkpos('699', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
[ura_spos, ura_ltime] = cspice_spkpos('799', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
[nep_spos, nep_ltime] = cspice_spkpos('899', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
[plu_spos, plu_ltime] = cspice_spkpos('999', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
[voyager_spos, ltime_end] = cspice_spkpos('-32', et_arr, 'ECLIPJ2000', 'LT', 'SUN');
[triton_spos, triton_ltime] = cspice_spkpos('801', et_arr, 'ECLIPJ2000', 'LT', 'SUN');


% Convert coordinates from [km] to [AU]
AU_km = 1.495979e+8;    % 1 AU in [km]
mar_spos = mar_spos/AU_km;
jup_spos = jup_spos/AU_km;
sat_spos = sat_spos/AU_km;
ura_spos = ura_spos/AU_km;
nep_spos = nep_spos/AU_km;
plu_spos = plu_spos/AU_km;
voyager_spos=voyager_spos/AU_km;
triton_spos = triton_spos/AU_km;

% Plot the orbits of Solar system
figure(1); 
plot(0,0,'-ko','MarkerSize',10,'Color','y',LineWidth=5); % Plot sun in center
hold on;
plot(voyager_spos(1,:),voyager_spos(2,:),'r');
plot(triton_spos(1,:),triton_spos(2,:),'g');
plot(mar_spos(1,:),mar_spos(2,:),'b');
plot(jup_spos(1,:),jup_spos(2,:),'b');
plot(sat_spos(1,:),sat_spos(2,:),'b');
plot(ura_spos(1,:),ura_spos(2,:),'b');
plot(nep_spos(1,:),nep_spos(2,:),'b');
plot(plu_spos(1,:),plu_spos(2,:),'b');
%plotFullPlanetOrbits(AU_km);
% plot(jup_spos(1,i_perih),jup_spos(2,i_perih),'-ko','MarkerSize',10,'Color','b'); % Plot Jup position at perhelion
hold off;

axis([-20 40 -120 20]);  % Set axes ranges
xlabel('AU')
ylabel('AU')
legend(["Sun" "Voyager" "Triton" "Planetary orbits"],'location','southwest')

[distanceTritonFlyby,dateTritonFlyby] = TritonFlyby(voyager_spos,triton_spos,et_arr);
% Question4(et_mission)

%% ------------------------- FUNCTIONS ------------------------------------
function [distance_flyby_precise,date_flyby] = TritonFlyby(voyager,triton,et_arr)
% Calculate distance and time of closest approach to Triton for Voyager

% OUTPUT:
% distance_flyby_precise = [km] Closest distance between Voyager and Triton
% date_flyby = Date of closest approach

    distance = sqrt(sum((voyager-triton).^2)); % [AU]
    [distance_flyby,index_flyby] = min(distance);

    % Re-do calculations with better precision
    time_step_precise = 60; % time step size [s]
    et_arr_precise = et_arr(index_flyby-1):time_step_precise:et_arr(index_flyby+1);
    
    % Precise positions for Voyager and Triton
    voyager_precise = cspice_spkpos('-32', et_arr_precise, 'ECLIPJ2000', 'LT', 'SUN');
    triton_precise = cspice_spkpos('801', et_arr_precise, 'ECLIPJ2000', 'LT', 'SUN');
   
    distance_precise = sqrt(sum((voyager_precise-triton_precise).^2));
    [distance_flyby_precise,index_flyby_precise] = min(distance_precise);
    time_flyby_precise_et = et_arr_precise(index_flyby_precise);
    
    % Convert et-time to date
    date_flyby = cspice_et2utc(time_flyby_precise_et,'C',0);
end


function [] = Question4(et)
% 'VOYAGER 2': code -32

    cnfine = cspice_wninsd( et(1), et(2) );
    % Set search parameters. Select a 3-minute step.
    MAXWIN  = 1000; %dimension of workspace array

    % %MERCURY
    % 
    %       occtyp  = 'ANY'; %any kind of ellipse, total, partial etc
    %       front   = 'MERCURY';
    %       fshape  = 'ELLIPSOID';
    %       fframe  = 'IAU_MERCURY'; %frame of body in the front
    %       back    = 'SUN';
    %       bshape  = 'ELLIPSOID';
    %       bframe  = 'IAU_SUN'; %frame of body in the back
    %       obsrvr  = 'VOYAGER 2'; %observer
    %       abcorr  = 'LT+S'; %type of correction, chosen from the slides
    % 
    %       step    = 180.; %time step of 3 s
    % 
    %        result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
    % 			      back, bshape, bframe,          ...
    % 			      abcorr, obsrvr, step, cnfine,  ...
    % 			      MAXWIN); %on the columns we read the time interval in which occultation occurred?
    
    
    % %VENUS
    % 
    %       occtyp  = 'ANY'; %any kind of ellipse, total, partial etc
    %       front   = 'VENUS';
    %       fshape  = 'ELLIPSOID';
    %       fframe  = 'IAU_VENUS'; %frame of body in the front
    %       back    = 'SUN';
    %       bshape  = 'ELLIPSOID';
    %       bframe  = 'IAU_SUN'; %frame of body in the back
    %       obsrvr  = 'VOYAGER 2'; %observer
    %       abcorr  = 'LT+S'; %type of correction
    % 
    %       step    = 180.; %time step of 3 s
    % 
    %        result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
    % 			      back, bshape, bframe,          ...
    % 			      abcorr, obsrvr, step, cnfine,  ...
    %  MAXWIN);
    
    
    %%EARTH
    %
    %       occtyp  = 'ANY'; %any kind of ellipse, total, partial etc
    %       front   = 'EARTH';
    %       fshape  = 'ELLIPSOID';
    %       fframe  = 'IAU_EARTH'; %frame of body in the front
    %       back    = 'SUN';
    %       bshape  = 'ELLIPSOID';
    %       bframe  = 'IAU_SUN'; %frame of body in the back
    %       obsrvr  = 'VOYAGER 2'; %observer
    %       abcorr  = 'LT+S'; %type of correction
    % 
    %       step    = 180.; %time step of 3 s
    % 
    %        result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
    % 			      back, bshape, bframe,          ...
    % 			      abcorr, obsrvr, step, cnfine,  ...
    %  MAXWIN);
    
    
    %MARS
    
          occtyp  = 'ANY'; %any kind of ellipse, total, partial etc
          front   = 'MARS';
          fshape  = 'ELLIPSOID';
          fframe  = 'IAU_MARS'; %frame of body in the front
          back    = 'SUN';
          bshape  = 'ELLIPSOID';
          bframe  = 'IAU_SUN'; %frame of body in the back
          obsrvr  = 'VOYAGER 2'; %observer
          abcorr  = 'LT+S'; %type of correction
    
          step    = 180.; %time step of 3 s
    
           result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
			          back, bshape, bframe,          ...
			          abcorr, obsrvr, step, cnfine,  ...
     MAXWIN);
    
       
    % %JUPITER
    % 
    %       occtyp  = 'ANY';
    %       front   = 'JUPITER';
    %       fshape  = 'ELLIPSOID';
    %       fframe  = 'IAU_JUPITER';
    %       back    = 'SUN';
    %       bshape  = 'ELLIPSOID';
    %       bframe  = 'IAU_SUN';
    %       obsrvr  = 'VOYAGER 2';
    %       abcorr  = 'LT+S';
    % 
    %       step    = 180.; %s
    % 
    %        result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
    % 			      back, bshape, bframe,          ...
    % 			      abcorr, obsrvr, step, cnfine,  ...
    % 			      MAXWIN);
    
    
    %%SATURN;
    %
    % 
    %       occtyp  = 'ANY';
    %       front   = 'SATURN';
    %       fshape  = 'ELLIPSOID';
    %       fframe  = 'IAU_SATURN';
    %       back    = 'SUN';
    %       bshape  = 'ELLIPSOID';
    %       bframe  = 'IAU_SUN';
    %       obsrvr  = 'VOYAGER 2';
    %       abcorr  = 'LT+S';
    % 
    %       step    = 180.;
    % 
    %        result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
    % 			      back, bshape, bframe,          ...
    % 			      abcorr, obsrvr, step, cnfine,  ...
    % 			      MAXWIN);
    
    
    %%URANUS
    %
    %       occtyp  = 'ANY'; %any kind of ellipse, total, partial etc
    %       front   = 'URANUS';
    %       fshape  = 'ELLIPSOID';
    %       fframe  = 'IAU_URANUS'; %frame of body in the front
    %       back    = 'SUN';
    %       bshape  = 'ELLIPSOID';
    %       bframe  = 'IAU_SUN'; %frame of body in the back
    %       obsrvr  = 'VOYAGER 2'; %observer
    %       abcorr  = 'LT+S'; %type of correction
    % 
    %       step    = 180.; %time step of 3 s
    % 
    %        result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
    % 			      back, bshape, bframe,          ...
    % 			      abcorr, obsrvr, step, cnfine,  ...
    %  MAXWIN);
    
    
    %NEPTUNE
    % 
    %       occtyp  = 'ANY'; %any kind of ellipse, total, partial etc
    %       front   = 'NEPTUNE';
    %       fshape  = 'ELLIPSOID';
    %       fframe  = 'IAU_NEPTUNE'; %frame of body in the front
    %       back    = 'SUN';
    %       bshape  = 'ELLIPSOID';
    %       bframe  = 'IAU_SUN'; %frame of body in the back
    %       obsrvr  = 'VOYAGER 2'; %observer
    %       abcorr  = 'LT+S'; %type of correction
    % 
    %       step    = 180.; %time step of 3 s
    % 
    %        result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
    % 			      back, bshape, bframe,          ...
    % 			      abcorr, obsrvr, step, cnfine,  ...
    % 			      MAXWIN); %on the columns we read the time interval in which occultation occurred?
    
    
    %time interval (to run separately for each planet)
    
    deltat_vec = 0;
    
     for i = 1: length(result(1,:))
        deltat = result(2,i)-result(1,i);
        deltat_vec=[deltat_vec, deltat];
        deltat_tot= sum(deltat_vec);
     end
    
      deltat_tot=deltat_tot/3600
    
     t = hours(deltat_tot)
    
    % disp('total eclipse time is:')
    % DON'T NAME THIS 'hours', IT OVERWRITES THE
    % MATLAB FUNCTION 
    % [hours,minutes,seconds]=hms(t) 
    
    %time conversion for dates of eclipses
    TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
    for i=1:numel(result)/2
     
         [left, right] = cspice_wnfetd( result, i );
    
         output = cspice_timout( [left,right], TIMFMT );
    
         if( isequal( left, right) )
    
            disp( ['Event time: ' output(1,:)] )
    
         else
    
            disp( ['From : ' output(1,:)] )
            disp( ['To   : ' output(2,:)] )
            disp( ' ')
    
         end
    end
end


function [] = plotFullPlanetOrbits(AU_km)
    % Function to plot full planet orbits fo all planets outwards 
    % from Mars, including Pluto

    planets = ['499'; '599'; '699'; '799'; '899'; '999'];
    orbitalPeriods = [687, 12*365, 30*365, 84*365, 165*365, 248*365]*86400; %[seconds]

    for planetNumber = 1:6
        if planets(planetNumber,:) == '899'
            start_time = -90*365*86400;
        elseif planets(planetNumber,:) == '999'
            start_time = -100*365*86400;
        else
            start_time = 0;
        end
        et_time = linspace(start_time,orbitalPeriods(planetNumber)+start_time,100);
        orbit = cspice_spkpos(planets(planetNumber), et_time, 'ECLIPJ2000', 'LT', 'SUN');
        orbit = orbit/AU_km;
        plot(orbit(1,:),orbit(2,:),'b')
    end
end



