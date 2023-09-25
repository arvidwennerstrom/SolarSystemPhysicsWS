close all
%% ------------------------------------------------------------------------
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

%% General Calculations

planetaryNames = ["Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"];
planetaryIDs = ['499'; '599'; '699'; '799'; '899'; '999'];

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

all_planet_spos = [mar_spos; jup_spos; sat_spos; ura_spos; nep_spos; plu_spos];

%% Question 1 & 2
figure(); 
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
hold off;
axis([-20 40 -120 20]);  % Set axes ranges
xlabel('AU')
ylabel('AU')
legend(["Sun" "Voyager" "Triton" "Planetary orbits"],'location','southwest')
title('Voyager 2 and planetary positions projected on the ecliptic plane')


%% Question 3
disp(" "); disp("Question 3")
disp("Voyager 2 made its closest approach to the following planets: ")
for planet_no = 1:length(planetaryIDs)
    planet_spos = all_planet_spos(3*planet_no-2:3*planet_no,:);
    [distanceFlyby,dateFlyby] = ClosestApproach(voyager_spos,planet_spos,planetaryIDs(planet_no),et_arr);
    disp([planetaryNames(planet_no)] + ": At a distance of " + num2str(round(distanceFlyby)) + " km, on " + dateFlyby)
end


%% Question 4
disp(" "); disp("Question 4")
Question4


%% Question 5
[distanceTritonFlyby,dateTritonFlyby] = ClosestApproach(voyager_spos,triton_spos,'801',et_arr);
disp("Question 5")
disp(['Voyager 2 approached to ', num2str(floor(distanceTritonFlyby)),' km from Triton on the following date :', dateTritonFlyby, 'by date of receipt of data on Earth' ])
disp(['Triton time is just 4 hours before beacause of the distance.'])


%% Question 6
disp(" "); disp("Question 6")
radii  = cspice_bodvrd( '801', 'RADII', 3 );
re = radii(1);
rp = radii(3);
f = ( re-rp)/re;

et_mission  = cspice_str2et(['1977-08-20T15:45';'2023-09-30T12:00']); % From start of the observation.
et_step  = 4*3600 ; % time step of 1 hour in ephemeris (in seconds)
et_arr   = et_mission(1):et_step:et_mission(2); % Create a time array from start to end

et_view_bounds = cspice_str2et(['1989-08-25T01:22';'1989-08-26T01:22']);
et_view_step  = 60 ;
et_view   = et_view_bounds(1):et_view_step:et_view_bounds(2);
[ spoint, trgepc, srfvec] = cspice_subpnt ( 'Intercept: ellipsoid', '801', et_view,        ...
                            'IAU_TRITON', 'LT','-32');
 
[spglon, spglat, spgalt] = cspice_recpgr( '801', spoint, re,f );
spglon = spglon * cspice_dpr;
spglat = spglat * cspice_dpr;


figure();
plot(spglon(:),spglat(:),'r');            

figure();
subplot(2, 1, 1);
plot(spglat, 'r-', 'LineWidth', 2);
title('Latitude');
xlabel('Temps');
ylabel('Latitude (degrés)');

subplot(2, 1, 2);
plot(spglon, 'g-', 'LineWidth', 2);
title('Longitude');
xlabel('Temps');
ylabel('Longitude (degrés)');


%% Question 7
disp(" "); disp("Question 7")
[subSolarLon,subSolarLat] = Question7(et_view,re,f);


%% Questions 8
distance_view = sqrt(sum((voyager_view-triton_view).^2));
TritonMap(distance_view',spglon,spglat,subSolarLon,subSolarLat);


%% Question 9
% Définir la date de début et de fin
dateDebut = cspice_str2et('2023-09-21T00:00:00');
dateFin = cspice_str2et('2023-09-21T23:59:59');

% Définir la position de l'observateur (Stockholm)
%[state, lighttime] = cspice_spkezr('Triton', dateDebut, 'J2000', 'LT+S', 'Earth').
%observerState = cspice_spkssb('Earth', dateDebut, 'J2000');

% Initialiser les tableaux pour stocker les résultats
timeArray = dateDebut:60:dateFin; % Une minute d'intervalle
nPoints = length(timeArray);

declination = zeros(1, nPoints);
rightAscension = zeros(1, nPoints);
altitude = zeros(1, nPoints);
azimuth = zeros(1, nPoints);

% Boucle sur chaque instant de temps
%for i = 1:nPoints
%    et = timeArray(i);
%    [state, lighttime] = cspice_spkezr('Triton', et, 'J2000', 'LT+S', 'Earth');
%    
%    % Calculer les coordonnées dans le système de coordonnées observer-centered
%    [range, ra, dec] = cspice_recrad(state - observerState);
%    [alt, az] = cspice_radrec(range, ra, dec);
%
%    % Stocker les coordonnées dans les tableaux
%    declination(i) = dec;
%    rightAscension(i) = ra;
%    altitude(i) = alt;
%    azimuth(i) = az;
%end

% Tracer les graphiques (comme précédemment)
%figure();

% Déclinaison et Ascension Droite
%subplot(2, 1, 1);
%plot(timeArray, rad2deg(declination), 'b', 'LineWidth', 1.5);
%hold on;
%plot(timeArray, rad2deg(rightAscension), 'r', 'LineWidth', 1.5);
%xlabel('Temps');
%ylabel('Degrés');
%title('Déclinaison et Ascension Droite de Triton');
%legend('Déclinaison', 'Ascension Droite');
%grid on;

% Altitude et Azimut
%subplot(2, 1, 2);
%plot(timeArray, rad2deg(altitude), 'g', 'LineWidth', 1.5);
%hold on;
%plot(timeArray, rad2deg(azimuth), 'm', 'LineWidth', 1.5);
%xlabel('Temps');
%ylabel('Degrés');
%title('Altitude et Azimut de Triton');
%legend('Altitude', 'Azimut');
%grid on;

%% ------------------------- FUNCTIONS ------------------------------------
function [distance_flyby_precise,date_flyby] = ClosestApproach(voyager,planet,planetID,et_arr)
% Calculate distance and time of closest approach to Triton for Voyager

% OUTPUT:
% distance_flyby_precise = [km] Closest distance between Voyager and Triton
% date_flyby = Date of closest approach

    distance = sqrt(sum((voyager-planet).^2)); % [AU]
    [distance_flyby,index_flyby] = min(distance);

    % Re-do calculations with better precision
    time_step_precise = 60; % time step size [s]
    et_arr_precise = et_arr(index_flyby-1):time_step_precise:et_arr(index_flyby+1);
    
    % Precise positions for Voyager and Planet
    voyager_precise = cspice_spkpos('-32', et_arr_precise, 'ECLIPJ2000', 'LT', 'SUN');
    planet_precise = cspice_spkpos(planetID, et_arr_precise, 'ECLIPJ2000', 'LT', 'SUN');
   
    distance_precise = sqrt(sum((voyager_precise-planet_precise).^2));
    [distance_flyby_precise,index_flyby_precise] = min(distance_precise);
    time_flyby_precise_et = et_arr_precise(index_flyby_precise);
    
    % Convert et-time to date
    date_flyby = cspice_et2utc(time_flyby_precise_et,'C',0);
end

function [] = Question4
% 'VOYAGER 2': code -32

    %cnfine = cspice_wninsd( et(1), et(2) );
    result = ones(2, 5);

    % GENERAL PARAMETERS FOR SUN
    occtyp  = 'ANY'; %any kind of ellipse, total, partial etc
    back    = 'SUN';
    bshape  = 'ELLIPSOID';
    bframe  = 'IAU_SUN'; %frame of body in the back
    obsrvr  = 'VOYAGER 2'; %observer
    abcorr  = 'LT+S'; %type of correction
    step    = 180.; %time step of 3 s
    MAXWIN  = 1000; %dimension of workspace array

    %MARS
    et_mars  = cspice_str2et(['1977-08-24T12:00';'1980-12-10T12:00']);
    cnfine_mars = cspice_wninsd( et_mars(1), et_mars(2) );
    front_mars   = 'MARS';
    fshape_mars  = 'ELLIPSOID';
    fframe_mars  = 'IAU_MARS'; %frame of body in the front
    
    result(:,1) = cspice_gfoclt( occtyp, front_mars, fshape_mars, fframe_mars, ...
			          back, bshape, bframe,          ...
			          abcorr, obsrvr, step, cnfine_mars, MAXWIN);
    
       
    %JUPITER
    et_jupiter  = cspice_str2et(['1978-11-24T12:00';'1981-05-10T12:00']);
    cnfine_jupiter = cspice_wninsd( et_jupiter(1), et_jupiter(2) );
    front_jupiter   = 'JUPITER';
    fshape_jupiter  = 'ELLIPSOID';
    fframe_jupiter  = 'IAU_JUPITER';

    result(:,2) = cspice_gfoclt( occtyp, front_jupiter, fshape_jupiter, fframe_jupiter, ...
     			      back, bshape, bframe,          ...
     			      abcorr, obsrvr, step, cnfine_jupiter, MAXWIN);
    
    %SATURN;
    et_saturn  = cspice_str2et(['1980-10-24T12:00';'1982-10-10T12:00']);
    cnfine_saturn = cspice_wninsd( et_saturn(1), et_saturn(2) );
    front_saturn   = 'SATURN';
    fshape_saturn  = 'ELLIPSOID';
    fframe_saturn  = 'IAU_SATURN';

    result(:,3) = cspice_gfoclt( occtyp, front_saturn, fshape_saturn, fframe_saturn, ...
     			      back, bshape, bframe,          ...
     			      abcorr, obsrvr, step, cnfine_saturn, MAXWIN);
    
    
    %URANUS
    et_uranus  = cspice_str2et(['1984-10-24T12:00';'1986-10-10T12:00']);
    cnfine_uranus = cspice_wninsd( et_uranus(1), et_uranus(2) );
    front_uranus   = 'URANUS';
    fshape_uranus  = 'ELLIPSOID';
    fframe_uranus  = 'IAU_URANUS'; %frame of body in the front

    result(:,4) = cspice_gfoclt( occtyp, front_uranus, fshape_uranus, fframe_uranus, ...
     			      back, bshape, bframe,          ...
     			      abcorr, obsrvr, step, cnfine_uranus, MAXWIN);
    
    %NEPTUNE
    et_neptune  = cspice_str2et(['1988-10-24T12:00';'1990-04-10T12:00']);
    cnfine_neptune = cspice_wninsd( et_neptune(1), et_neptune(2) );
    front_neptune   = 'NEPTUNE';
    fshape_neptune  = 'ELLIPSOID';
    fframe_neptune  = 'IAU_NEPTUNE'; %frame of body in the front

    result(:,5) = cspice_gfoclt( occtyp, front_neptune, fshape_neptune, fframe_neptune, ...
     			      back, bshape, bframe,          ...
     			      abcorr, obsrvr, step, cnfine_neptune, MAXWIN);
    
    
    %time interval (to run separately for each planet)
    duration_hours=zeros(1,5);
    TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
    planets = {'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};

    for i=1:5
        duration_hours(i) = (result(2,i)-result(1,i))/3600;   %duration in hours
        for j=1:numel(result(:,i))/2
     
         [left, right] = cspice_wnfetd( result(:,i), j );
    
         output = cspice_timout( [left,right], TIMFMT );
    
         if( isequal( left, right) )
    
            disp( ['Event time: ' output(1,:)] )
            disp( ' ')
    
         else
            disp( planets{i})
            disp( ['From : ' output(1,:)] )
            disp( ['To   : ' output(2,:)] )
            disp( ['Duration :' num2str(duration_hours(i))])
            disp( ' ')
    
         end
        end
    end
    
    %deltat_vec = 0;
    
     %for i = 1: length(result(1,:))
     %   deltat = result(2,i)-result(1,i);
     %   deltat_vec=[deltat_vec, deltat];
     %   deltat_tot= sum(deltat_vec);
     %end
    
      %deltat_tot=deltat_tot/3600
    
     %t = hours(deltat_tot)
    
    % disp('total eclipse time is:')
    % DON'T NAME THIS 'hours', IT OVERWRITES THE
    % MATLAB FUNCTION 
    % [hours,minutes,seconds]=hms(t) 
    
    %time conversion for dates of eclipses
    
    %for i=1:numel(result)/2
    % 
    %     [left, right] = cspice_wnfetd( result, i );
    %
    %     output = cspice_timout( [left,right], TIMFMT );
    %
    %     if( isequal( left, right) )
    %
    %        disp( ['Event time: ' output(1,:)] )
    %
    %     else
    %
    %        disp( ['From : ' output(1,:)] )
    %        disp( ['To   : ' output(2,:)] )
    %        disp( ' ')
    %
    %     end
    %end
end

function [] = TritonMap(altitude,lon,lat,lon2,lat2)
    triton_r = 1353.4; % [km]
    figure()    
    r = ones(length(lon),1)*triton_r;
            
    [x1,y1,z1] = sph2cart(deg2rad(lon),deg2rad(lat),r(1));
    [x2,y2,z2] = sph2cart(deg2rad(lon2),deg2rad(lat2),r(1));
    [xs,ys,zs] = sphere(30);
    xs = xs*triton_r; ys = ys*triton_r; zs = zs*triton_r;
   
    plot3(x2,y2,z2,'.r')
    hold on
    scatter3(x1,y1,z1,20,altitude,'filled')
    % textscatter(x1(1),y1(1),z1(1),'t = 0')
    surface(xs,ys,zs,'facecolor','none','edgecolor',ones(1,3)*0.3)
    
    view(45,22.5)
    c = colorbar('eastoutside');
    c.Label.String = 'Voyager 2 altitude above Triton surface [km]';
    axis equal
    xlabel('x-coordinate [km]'); ylabel('y-coordinate [km]'); zlabel('z-coordinate [km]')
    title('Map of Triton')
    legend('Sub-solar point','Location','north')
end

function [spglon,spglat] = Question7(et_view,re,f)


    [spoint] = cspice_subsol('NEARPOINT','801', et_view, 'LT', '-32');
    
    [spglon, spglat, spgalt] = cspice_recpgr( '801', spoint, re,f );
    spglon = spglon * cspice_dpr;
    spglat = spglat * cspice_dpr;
    
    figure()
    subplot(2, 1, 1);
    plot(spglat, 'r-', 'LineWidth', 2);
    title('Latitude');
    xlabel('Temps');
    ylabel('Latitude (degrés)');
    
    subplot(2, 1, 2);
    plot(spglon, 'g-', 'LineWidth', 2);
    title('Longitude');
    xlabel('Temps');
    ylabel('Longitude (degrés)');
    
    for i=1:length(spglon)
        if spglon(i)>250
            spglon(i) = spglon(i)-360;
        end
    end
    
    figure();
    plot(spglon(:),spglat(:),'r'); 
    title('Trajectory of sub-solar point on Triton during closest approach');
    xlabel('Longitude (degrees)');
    ylabel('Latitude (degrees)');
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

