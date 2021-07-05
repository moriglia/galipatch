close all ;
clc ;

%% Load pre-evaluated slow-to-evaluate variables
if exist('slowVariables.mat', 'file')
    load('slowVariables.mat');
end

%% Galileo Frequencies and Bandwidths
[f0, B] = utils.galileognss('E1');
c0 = 3e8;
lambda0 = c0/f0;

LBand = 1200e6 : 500e3 : 1800e6;

% design_band = LBand;
design_band = 1300e6 : 1e6 : 1700e6 ;
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

figure(1)
show(truncatedCornerPatch);
title("Corner Truncated Patch for Galileo E1 Reception");

if exist('Z', 'var') == 0
    Z = impedance(truncatedCornerPatch, design_band);
end
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

RHCP_xz = pattern(truncatedCornerPatch, f0, 0, elevation_span, 'Polarization', 'RHCP');
LHCP_xz = pattern(truncatedCornerPatch, f0, 0, elevation_span, 'Polarization', 'LHCP');
RHCP_yz = pattern(truncatedCornerPatch, f0, 90, elevation_span, 'Polarization', 'RHCP');
LHCP_yz = pattern(truncatedCornerPatch, f0, 90, elevation_span, 'Polarization', 'LHCP');

figure(6);
polaraxes('ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top', 'rlim', [-25, 10]);
hold on;
polarplot(zenith_span, RHCP_xz, 'g');
polarplot(zenith_span, LHCP_xz, 'r--');
title("Copolar and Crosspolar components @ E1 carrier, xz-plane (\phi=0)");
legend("RHCP", "LHCP");
hold off;

figure(7);
polaraxes('ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top', 'rlim', [-30, 10]);
hold on;
polarplot(zenith_span, RHCP_yz, 'g');
polarplot(zenith_span, LHCP_yz, 'r--');
title("Copolar and Crosspolar components @ E1 carrier, yz-plane (\phi=90°)");
legend("RHCP", "LHCP");
hold off;


%% Axial Ratio
ARxz        = axialRatio(truncatedCornerPatch, f0, 0, [0:180]);
ARyz        = axialRatio(truncatedCornerPatch, f0, 90, [0:180]);
if exist('ARboresight', 'var') == 0
    ARboresight = axialRatio(truncatedCornerPatch, design_band, 0, 90);
end

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
plot(design_band*1e-6, ARboresight);
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

%% Field on a z plane
X_v=linspace(-0.5, 0.5, 201)*10*Lp;
Y_v=linspace(-0.5, 0.5, 201)*10*Wp;
[X,Y] = meshgrid(X_v, Y_v);

z_quote=1.5*h_sub;
field_grid = [X(:)'; Y(:)'; z_quote*ones(1,numel(X))];

if exist('E', 'var') == 0 || exist('H', 'var') == 0
    [E, H]=EHfields(truncatedCornerPatch, f0, field_grid);
end

% Only tangent component will be used in the equivalence principle
E_tangent_modulus = sqrt(abs(E(1,:)).^2+abs(E(2,:)).^2);
E_tangent_modulus_normalized = 10*log10(E_tangent_modulus / max(E_tangent_modulus, [], 'all'));
H_tangent_modulus = sqrt(abs(H(1,:)).^2+abs(H(2,:)).^2);
H_tangent_modulus_normalized = 10*log10(H_tangent_modulus / max(H_tangent_modulus, [], 'all'));

figure(12);
[c, h]=contourf(X_v*1e3, Y_v*1e3,reshape(E_tangent_modulus_normalized,length(X_v),length(Y_v)));
axis equal;
clabel(c, h);
colorbar ;
xlabel('X-axis [mm]')
ylabel('Y-axis [mm]')
title(['Tangent E field @ z=',num2str(z_quote*1e3),'mm']);
hold on;
fill([-1, 1, 1, -1]*100, [-1, -1, 1, 1]*100, 'r--');
alpha(0.05);
hold off;

figure(13);
[c, h]=contourf(X_v*1e3, Y_v*1e3,reshape(H_tangent_modulus_normalized,length(X_v),length(Y_v)));
axis equal;
clabel(c, h);
colorbar ;
xlabel('X-axis [mm]')
ylabel('Y-axis [mm]')
title(['Tangent H field @ z=',num2str(z_quote*1e3),'mm']);
hold on;
fill([-1, 1, 1, -1]*100, [-1, -1, 1, 1]*100, 'r--');
alpha(0.05);
hold off;

% The tangent component of the fields is very low out of a square of side
% length of 20 cm, hence I will neglect the fields outside that square. So
% I can redefine my problem on a restricted domain for the equivalence
% principle application.

%% Equivalence Principle Application
% Reduced domain
X_v_reduced=X_v(abs(X_v) < 10e-2);
Y_v_reduced=Y_v(abs(Y_v) < 10e-2);

% Reduce domain
% utils.fixnan cleans up some NaN's that are yielded by EHfields
E_x = reshape(E(1,:), length(Y_v), length(X_v));
E_x = utils.fixnan(E_x);
E_x = E_x(abs(Y_v) < 10e-2, abs(X_v) < 10e-2);
E_y = reshape(E(2,:), length(Y_v), length(X_v));
E_y = utils.fixnan(E_y);
E_y = E_y(abs(Y_v) < 10e-2, abs(X_v) < 10e-2);

H_x = reshape(H(1,:), length(Y_v), length(X_v));
H_x = utils.fixnan(H_x);
H_x = H_x(abs(Y_v) < 10e-2, abs(X_v) < 10e-2);
H_y = reshape(H(2,:), length(Y_v), length(X_v));
H_y = utils.fixnan(H_y);
H_y = H_y(abs(Y_v) < 10e-2, abs(X_v) < 10e-2);

dX = X_v_reduced(2) - X_v_reduced(1);
dY = Y_v_reduced(2) - Y_v_reduced(1);


zenith_angles = (-1: .025 :1)*pi/2;
phase_constant = 2*pi/lambda0;

eta_medium = 377;

% XZ plane
kernel_xz_plane = exp(1i*phase_constant*sin(zenith_angles.')*X_v_reduced);
P_x_xz = kernel_xz_plane*sum(E_x, 1).'*dX*dY;
P_y_xz = kernel_xz_plane*sum(E_y, 1).'*dX*dY;

Q_x_xz = kernel_xz_plane*sum(H_x, 1).'*dX*dY;
Q_y_xz = kernel_xz_plane*sum(H_y, 1).'*dX*dY;

E_theta_xz = eta_medium*cos(zenith_angles).*Q_y_xz.' + P_x_xz.';
E_phi_xz = -eta_medium*Q_x_xz.' + cos(zenith_angles).*P_y_xz.';

E_theta_xz = 1j*phase_constant/4/pi*E_theta_xz;
E_phi_xz = 1j*phase_constant/4/pi*E_phi_xz;

% YZ plane
kernel_yz_plane = exp(1j*phase_constant*sin(zenith_angles.')*Y_v_reduced);
P_x_yz = kernel_yz_plane*sum(E_x, 2)*dX*dY;
P_y_yz = kernel_yz_plane*sum(E_y, 2)*dX*dY;

Q_x_yz = kernel_yz_plane*sum(H_x, 2)*dX*dY;
Q_y_yz = kernel_yz_plane*sum(H_y, 2)*dX*dY;

E_theta_yz = -eta_medium*cos(zenith_angles).*Q_x_yz.' + P_y_yz.';
E_phi_yz = -eta_medium*Q_y_yz.' - cos(zenith_angles).*P_x_yz.';

E_theta_yz = 1j*phase_constant/4/pi*E_theta_yz;
E_phi_yz = 1j*phase_constant/4/pi*E_phi_yz;


%% Pattern comparison
E_pattern_xz = 10*log10(abs(E_theta_xz).^2 + abs(E_phi_xz).^2);
E_pattern_xz = E_pattern_xz - max(E_pattern_xz);  % Normalize

E_pattern_yz = 10*log10(abs(E_theta_yz).^2 + abs(E_phi_yz).^2);
E_pattern_yz = E_pattern_yz - max(E_pattern_yz);  % Normalize

Full_pattern_xz = pattern(truncatedCornerPatch, f0, 0, linspace(0, 180, numel(zenith_angles)), 'Normalize', true);
Full_pattern_yz = pattern(truncatedCornerPatch, f0, 90, linspace(0, 180, numel(zenith_angles)), 'Normalize', true);

figure(14);
polaraxes('ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top', 'rlim', [-20, 0]);
hold on;
polarplot(zenith_angles, E_pattern_xz);
polarplot(zenith_angles, Full_pattern_xz);
title("Power Pattern @ E1 carrier, xz-plane (\phi=0)");
legend("Equivalence", "Full Wave");
hold off;

figure(15);
polaraxes('ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top', 'rlim', [-20, 0]);
hold on;
polarplot(zenith_angles, E_pattern_yz);
polarplot(zenith_angles, Full_pattern_yz);
title("Power Pattern @ E1 carrier, xz-plane (\phi=0)");
legend("Equivalence", "Full Wave");
hold off;

%% Co- and Cross-polar component comparison

if exist('NormalizationOffset', 'var') == 0
    RHCP_fw_xz = pattern(truncatedCornerPatch, f0, 0, elevation_span, 'Polarization', 'RHCP', 'Type', 'powerdb');
    LHCP_fw_xz = pattern(truncatedCornerPatch, f0, 0, elevation_span, 'Polarization', 'LHCP', 'Type', 'powerdb');
    RHCP_fw_yz = pattern(truncatedCornerPatch, f0, 90, elevation_span, 'Polarization', 'RHCP', 'Type', 'powerdb');
    LHCP_fw_yz = pattern(truncatedCornerPatch, f0, 90, elevation_span, 'Polarization', 'LHCP', 'Type', 'powerdb');
    
    NormalizationOffset = max(RHCP_fw_xz);
    
    % Use the same normalization to keep rato between RHCP and LHCP. Note
    % that on the two planes the RHCP component max power must be the same,
    % since the maximum directivity is @ boresight, of which the direction
    % is included in both planes.
    RHCP_fw_xz = RHCP_fw_xz - NormalizationOffset;
    RHCP_fw_yz = RHCP_fw_yz - NormalizationOffset;
    LHCP_fw_xz = LHCP_fw_xz - NormalizationOffset;
    LHCP_fw_yz = LHCP_fw_yz - NormalizationOffset;
    
    RHCP_fw_xz(RHCP_fw_xz < -30) = -30;
    LHCP_fw_xz(LHCP_fw_xz < -30) = -30;
    RHCP_fw_yz(RHCP_fw_yz < -30) = -30;
    LHCP_fw_yz(LHCP_fw_yz < -30) = -30;
end

RHCP_eq_xz = 20*log10(abs(sqrt(0.5)*(E_theta_xz + 1j*E_phi_xz)));
RHCP_eq_xz = RHCP_eq_xz - max(RHCP_eq_xz);
LHCP_eq_xz = 20*log10(abs(sqrt(0.5)*(E_theta_xz - 1j*E_phi_xz)));
LHCP_eq_xz = LHCP_eq_xz - max(RHCP_eq_xz);

RHCP_eq_yz = 20*log10(abs(sqrt(0.5)*(E_theta_yz + 1j*E_phi_yz)));
RHCP_eq_yz = RHCP_eq_yz - max(RHCP_eq_yz);
LHCP_eq_yz = 20*log10(abs(sqrt(0.5)*(E_theta_yz - 1j*E_phi_yz)));
LHCP_eq_yz = LHCP_eq_yz - max(RHCP_eq_yz);

figure(16);
polaraxes('ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top', 'rlim', [-30, 0]);
hold on;
polarplot(zenith_span, RHCP_fw_xz, 'k-.');
polarplot(zenith_span, LHCP_fw_xz, 'k:');
polarplot(zenith_angles, RHCP_eq_xz, 'g-');
polarplot(zenith_angles, LHCP_eq_xz, 'r--');
polarplot(zenith_angles, -3*ones(1,length(zenith_angles)), 'r.');
title("Copolar and Crosspolar components @ E1 carrier, xz-plane (\phi=0)");
legend("RHCP full wave", "LHCP full wave", "RHCP equivalence", "LHCP equivalence", "-3dB Curve");
hold off;

figure(17);
polaraxes('ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top', 'rlim', [-30, 0]);
hold on;
polarplot(zenith_span, RHCP_fw_yz, 'k-.');
polarplot(zenith_span, LHCP_fw_yz, 'k:');
polarplot(zenith_angles, RHCP_eq_yz, 'g-');
polarplot(zenith_angles, LHCP_eq_yz, 'r--');
polarplot(zenith_angles, -3*ones(1,length(zenith_angles)), 'r.');
title("Copolar and Crosspolar components @ E1 carrier, yz-plane (\phi=\pi/2)");
legend("RHCP full wave", "LHCP full wave", "RHCP equivalence", "LHCP equivalence", "-3dB Curve");
hold off;
