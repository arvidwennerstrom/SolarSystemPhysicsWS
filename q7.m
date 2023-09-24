%% Question 7

[spoint] = cspice_subsol('NEARPOINT','801', et_view, 'LT', '-32');

[spglon, spglat, spgalt] = cspice_recpgr( '801', spoint, re,f );
spglon = spglon * cspice_dpr;
spglat = spglat * cspice_dpr;


figure(4)
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

figure(5);
plot(spglon(:),spglat(:),'r'); 
title('Trajectory of sub-solar point on Triton during closest approach');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');