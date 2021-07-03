close all;
clear all;
clc;

%% Galileo Frequencies and Bandwidths
signalNameMap = containers.Map;
signalNameMap('E1')     = 1;
signalNameMap('E5a')    = 2;
signalNameMap('E5b')    = 3;
signalNameMap('E6')     = 4;

f0 = [1575.420 1176.450 1207.150 1278.750]*1e6;
B = [ 24.552 20.460 20.460 40.920]*1e6;

c0 = 3e8;
lambda0 = c0./f0;

LBand = [1200e6 : 1e6 : 1800e6];
R_0 = 50; % Ohm

%% Dielectric
d  = dielectric('Air');
h_sub = 2.95e-3;


% %% Rectangular Resonant patch for E1 signal
% [Lp, ~, ~] = utils.rectagularResonantPatchParams(lambda0(1), d.EpsilonR, h_sub);
% 
% Wp = Lp;
% 
%
% %% Truncated Corner
% 
% p = utils.cornerTruncatedPatch(Wp, Lp, Lp/10, 'RHCP');
% 
% groundPlane = antenna.Rectangle('Width', 3*Wp, 'Length', 3*Lp);
% 
% truncatedCornerPatch = pcbStack;
% truncatedCornerPatch.Name = 'Galileo  E1 Patch';
% truncatedCornerPatch.BoardThickness = h_sub;
% truncatedCornerPatch.BoardShape = groundPlane;
% truncatedCornerPatch.Layers = {p, d, groundPlane};
% truncatedCornerPatch.FeedLocations = [-0.4*Wp/2 0 1 3];
% 
% figure(1)
% show(truncatedCornerPatch);
% 
% Z = impedance(truncatedCornerPatch, LBand);
% reflectionCoefficient=(Z-R_0)./(Z + R_0);
% reflectionCoefficient_dB = 20*log10(abs(reflectionCoefficient));
% 
% figure(2)
% plot(LBand, -reflectionCoefficient_dB);
% set(gca, 'Ydir', 'reverse');
% 
% figure(3)
% smithplot(LBand, reflectionCoefficient);

%% Adjust resonance frequency
[Lp, ~, ~] = utils.rectagularResonantPatchParams(lambda0(1)*1.065, d.EpsilonR, h_sub);

Wp = Lp;

p = utils.cornerTruncatedPatch(Wp, Lp, Lp/10, 'RHCP');

groundPlane = antenna.Rectangle('Width', 1.5*Wp, 'Length', 1.5*Lp);

truncatedCornerPatch = pcbStack;
truncatedCornerPatch.Name = 'Galileo E1 Patch';
truncatedCornerPatch.BoardThickness = h_sub;
truncatedCornerPatch.BoardShape = groundPlane;
truncatedCornerPatch.Layers = {p, d, groundPlane};
truncatedCornerPatch.FeedLocations(3:4) = [1 3];

figure(4)
show(truncatedCornerPatch);


Z = impedance(truncatedCornerPatch, LBand);
reflectionCoefficient=(Z-R_0)./(Z + R_0);
reflectionCoefficient_dB = 20*log10(abs(reflectionCoefficient));

figure(2)
plot(LBand, -reflectionCoefficient_dB);
hold on ;
set(gca, 'Ydir', 'reverse');
utils.returnLossMask(f0(1)-B(1)/2, f0(1)+B(1)/2);
hold off ;

figure(3)
hold on
smithplot(LBand, reflectionCoefficient);

figure(5)
plot(LBand, real(Z), 'k--');
grid on
hold on
plot(LBand, imag(Z), 'r');
line([f0(1)-B(1)/2, f0(1) + B(1)/2], [0, 0], 'LineWidth', 3);
hold off
