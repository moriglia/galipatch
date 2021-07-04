close all;
clear all;
clc;

%% Galileo Frequencies and Bandwidths
[f0, B] = utils.galileognss('E1');
c0 = 3e8;
lambda0 = c0/f0;

LBand = [1200e6 : 500e3 : 1800e6];

% design_band = LBand;
design_band = [1300e6 : 1e6 : 1700e6];
R_0 = 50; % Ohm

%% Dielectric height
h_sub = 2e-3;


%% Patch Design
d  = dielectric('FR4');
[Lp, ~, ~] = utils.rectagularResonantPatchParams(lambda0, d.EpsilonR, h_sub);

Wp = Lp;

p = utils.cornerTruncatedPatch(Wp, Lp, Lp/8, 'RHCP');

groundPlane = antenna.Rectangle('Width', 1.3*Wp, 'Length', 1.3*Lp);

truncatedCornerPatch = pcbStack;
truncatedCornerPatch.Name = 'Galileo E1 Patch';
truncatedCornerPatch.BoardThickness = h_sub;
truncatedCornerPatch.BoardShape = groundPlane;
truncatedCornerPatch.Layers = {p, d, groundPlane};
truncatedCornerPatch.FeedLocations(3:4) = [1 3];

figure(4)
show(truncatedCornerPatch);
title("Corner Truncated Patch for Galileo E1 Reception");


Z = impedance(truncatedCornerPatch, design_band);
reflectionCoefficient=(Z-R_0)./(Z + R_0);
reflectionCoefficient_dB = 20*log10(abs(reflectionCoefficient));

figure(2)
plot(design_band/1e6, -reflectionCoefficient_dB);
hold on ;
set(gca, 'Ydir', 'reverse');
utils.returnLossMask((f0-B/2)/1e6, (f0+B/2)/1e6);
title('Return Loss in the L Band');
xlabel('Frequency f [MHz]');
ylabel('RL [dB]')
legend('RL', 'Mask E1 signal band', "Location", "SouthWest")
hold off ;

figure(3)
smithplot(design_band, reflectionCoefficient);
hold on;
fill(10^-0.5*cos((0:.01:1)*2*pi), 10^-0.5*sin((0:.01:1)*2*pi), 'g');
alpha(0.5);
title("Smith Chart");
legend('Band L Impedance/Reflection Coefficient', 'Acceptable region');
hold off;

figure(5)
plot(design_band*1e-6, real(Z), 'k--');
grid on
hold on
plot(design_band*1e-6, imag(Z), 'r');
line([f0-B/2, f0 + B/2]*1e-6, [0, 0], 'LineWidth', 3);
title('Input Impedance in the L band');
xlabel('Frequency f [MHz]');
ylabel("Impedance [\Omega]");
legend('\Re\{Z_{IN}\}', '\Im\{Z_{IN}\}', 'Galileo E1 Band', "Location", "SouthWest");
hold off;


%% Pattern
elevation_span = -180 : 1 : 180;
zenith_span = (90 - elevation_span)*pi/180;

RHCPxz = pattern(truncatedCornerPatch, f0, 0, elevation_span, 'Polarization', 'RHCP');
LHCPxz = pattern(truncatedCornerPatch, f0, 0, elevation_span, 'Polarization', 'LHCP');
RHCPyz = pattern(truncatedCornerPatch, f0, 90, elevation_span, 'Polarization', 'RHCP');
LHCPyz = pattern(truncatedCornerPatch, f0, 90, elevation_span, 'Polarization', 'LHCP');

figure(6);
polaraxes('ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top', 'rlim', [-25, 10]);
hold on;
polarplot(zenith_span, RHCPxz);
polarplot(zenith_span, LHCPxz);
title("Copolar and Crosspolar components @ E1 carrier, xz-plane (\phi=0)");
legend("RHCP", "LHCP");
hold off;

figure(7);
polaraxes('ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top', 'rlim', [-30, 10]);
hold on;
polarplot(zenith_span, RHCPyz);
polarplot(zenith_span, LHCPyz);
title("Copolar and Crosspolar components @ E1 carrier, yz-plane (\phi=90°)");
legend("RHCP", "LHCP");
hold off;


clear elevation_span zenith_span ;

%% Axial Ratio
ARxz        = axialRatio(truncatedCornerPatch, f0, 0, [0:180]);
ARyz        = axialRatio(truncatedCornerPatch, f0, 90, [0:180]);
AXboresight = axialRatio(truncatedCornerPatch, design_band, 0, 90);

figure(10);
plot((90:-1:-90), ARxz, 'k--');
hold on;
grid on;
plot((90:-1:-90), ARyz, 'k-.');
line([-45, 45], [3, 3], 'Color', 'Red', 'LineWidth', 2);
legend("xz-plane (\phi = 0) AR", "yz-plane (\phi = \pi/2) AR", "3dB AR limit");
title(['Axial ratio @ ', num2str(f0*1e-6), 'MHz']);
ylabel('AR [dB]');
xlabel('\theta [deg]');
hold off;

figure(11);
plot(design_band*1e-6, AXboresight);
hold on;
grid on;
line([design_band(1), design_band(end)]*1e-6, [3, 3], 'Color', 'Red', 'LineWidth', 2);
title("Axial Ratio over the L band @ boresight");
xlabel("Frequency f [MHz]");
ylabel("Axial Ratio [dB]");
ylim([0, 30]);
utils.axialRatioMask((f0-B/2)*1e-6, (f0+B/2)*1e-6);
legend("AR [dB] @ \theta = 0", "3dB Limit", "Galileo E1 Band");
hold off;


