
%% ------------------------------------------------------------------------
%
% Lorenz Roth / 1 Sep 2022 
% Example program for WS1 on SPICE / mice for EF2243 HT2022
%
%% ------------------------------------------------------------------------

% Your directory where SPICE toolkit and kernels are located
% ! Change to the correct directory on your computer ! 
dir_spice = '/Users/arvwe/OneDrive/Dokument/MATLAB/MATLAB_SolarSystem/';
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

% Extra kernels:
% -> Example: Kernel for position of the Hubble Space Telescope (HST)
%cspice_furnsh( fullfile(dir_spice, 'kernels/SPK/hst.bsp') ) 
%% ---

% Define a start and end date for a time interval you want to use 
% and convert the date to ephemeris time (Unit: seconds)
et_mission  = cspice_str2et(['2016-09-01T12:00';'2028-09-01T12:00']); % From start of the observation.
et_step  = 3600 ; % time step of 1 hour in ephemeris (in seconds)
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


% Convert coordinates from [km] to [AU]
AU_km = 1.495979e+8;    % 1 AU in [km]
mar_spos = mar_spos/AU_km;
jup_spos = jup_spos/AU_km;
sat_spos = sat_spos/AU_km;
ura_spos = ura_spos/AU_km;
nep_spos = nep_spos/AU_km;
plu_spos = plu_spos/AU_km;

% Find perihelion of Jupiter in that period (minimum distance to Sun)
% ! [distance in km , array index of minimum distance ]
[jup_perih,i_perih] = min(cspice_vnorm(jup_spos)); 
date_perih = cspice_et2utc(et_arr(i_perih),'C',0); % Convert ET of perihelion to date string
% disp(['Jupiter´s perihelion distance is ',num2str(jup_perih), ' AU, on ',date_perih(1:11)]);

% Plot the orbits of Solar system
outerLimit = 35;
figure(1); 
plot(0,0,'-ko','MarkerSize',10,'Color','y',LineWidth=5); % Plot sun in center
hold on;
plot(mar_spos(1,:),mar_spos(2,:),'b');
plot(jup_spos(1,:),jup_spos(2,:),'b');
plot(sat_spos(1,:),sat_spos(2,:),'b');
plot(ura_spos(1,:),ura_spos(2,:),'b');
plot(nep_spos(1,:),nep_spos(2,:),'b');
plot(plu_spos(1,:),plu_spos(2,:),'b');

% plot(jup_spos(1,i_perih),jup_spos(2,i_perih),'-ko','MarkerSize',10,'Color','b'); % Plot Jup position at perhelion
hold off;
axis([-outerLimit outerLimit -outerLimit outerLimit]);  % Set axes ranges
xlabel('AU')
ylabel('AU')

%--------------------------------------------------------------------------


