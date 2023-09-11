
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
% Standard kernels:
cspice_furnsh( fullfile(dir_spice, 'kernels/LSK/naif0010.tls') ) % Leapseconds - time conversion
cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/de440s.bsp') )   % Ephemeris - planets
cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/mar097.bsp') )   % Ephemeris - Mars + moons system
cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/jup230l.bsp') )  % Ephemeris - Jupiter + moons system
cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/sat415.bsp') )   % Ephemeris - Saturn + moons system
cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/ura111.bsp') )   % Ephemeris - Uranus + moons system
cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/nep095.bsp') )   % Ephemeris - Neptune + moons system
cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/plu058.bsp') )   % Ephemeris - Pluto
cspice_furnsh( fullfile(dir_spice, 'kernels/PCK/pck00010.tpc') ) % Orientation, size, shape - planets
cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/Voyager_2.m05016u.merged.bsp') )   % Voyager full trajectory

% Extra kernels:
% -> Example: Kernel for position of the Hubble Space Telescope (HST)
%cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/hst.bsp') ) 
%% ---

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

% Find perihelion of Jupiter in that period (minimum distance to Sun)
% ! [distance in km , array index of minimum distance ]
%[jup_perih,i_perih] = min(cspice_vnorm(jup_spos)); 
%date_perih = cspice_et2utc(et_arr(i_perih),'C',0); % Convert ET of perihelion to date string
% disp(['Jupiter´s perihelion distance is ',num2str(jup_perih), ' AU, on ',date_perih(1:11)]);

% Plot the orbits of Solar system
figure(1); 
plot(0,0,'-ko','MarkerSize',10,'Color','y',LineWidth=5); % Plot sun in center
hold on;
plot(mar_spos(1,:),mar_spos(2,:),'b');
plot(jup_spos(1,:),jup_spos(2,:),'b');
plot(sat_spos(1,:),sat_spos(2,:),'b');
plot(ura_spos(1,:),ura_spos(2,:),'b');
plot(nep_spos(1,:),nep_spos(2,:),'b');
plot(plu_spos(1,:),plu_spos(2,:),'b');
plot(voyager_spos(1,:),voyager_spos(2,:),'r');
plot(triton_spos(1,:),triton_spos(2,:),'g');
%plotFullPlanetOrbits(AU_km);

% plot(jup_spos(1,i_perih),jup_spos(2,i_perih),'-ko','MarkerSize',10,'Color','b'); % Plot Jup position at perhelion
hold off;
axis([-20 40 -120 20]);  % Set axes ranges
xlabel('AU')
ylabel('AU')

%--------------------------------------------------------------------------
function [] = plotFullPlanetOrbits(AU_km)
    
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



