close all ;
clc ;

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
if exist('AXboresight', 'var') == 0
    AXboresight = axialRatio(truncatedCornerPatch, design_band, 0, 90);
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

%% Field on a z plane
X_v=linspace(-0.5, 0.5, 201)*10*Lp;
Y_v=linspace(-0.5, 0.5, 201)*10*Wp;
[X,Y] = meshgrid(X_v, Y_v);

z_quote=1.5*h_sub;
field_grid = [X(:)'; Y(:)'; z_quote*ones(1,numel(X))];

if exist('E', 'var') == 0 || exist('H', 'var') == 0
    [E, H]=EHfields(truncatedCornerPatch, f0, field_grid);
end
% Efield=20*log10(sqrt(abs(E(1,:)).^2+abs(E(2,:)).^2+abs(E(3,:)).^2));
E_tangent_modulus = sqrt(abs(E(1,:)).^2+abs(E(2,:)).^2);
E_tangent_modulus_normalized = 10*log10(E_tangent_modulus / max(E_tangent_modulus, [], 'all'));
H_tangent_modulus = sqrt(abs(H(1,:)).^2+abs(H(2,:)).^2);
H_tangent_modulus_normalized = 10*log10(H_tangent_modulus / max(H_tangent_modulus, [], 'all'));

figure(12);
% cmin=10;
% cmax=45;
% v=[cmin:5:cmax];
[c, h]=contourf(X_v*1e3, Y_v*1e3,reshape(E_tangent_modulus_normalized,length(X_v),length(Y_v)));
axis equal;
clabel(c, h);
colorbar ;
% caxis([cmin,cmax]);
xlabel('X-axis [mm]')
ylabel('Y-axis [mm]')
title(['Tangent E field @ z=',num2str(z_quote*1e3),'mm']);
hold on;
fill([-1, 1, 1, -1]*100, [-1, -1, 1, 1]*100, 'r--');
alpha(0.05);
hold off;

figure(13);
% cmin=10;
% cmax=45;
% v=[cmin:5:cmax];
[c, h]=contourf(X_v*1e3, Y_v*1e3,reshape(H_tangent_modulus_normalized,length(X_v),length(Y_v)));
axis equal;
clabel(c, h);
colorbar ;
% caxis([cmin,cmax]);
xlabel('X-axis [mm]')
ylabel('Y-axis [mm]')
title(['Tangent H field @ z=',num2str(z_quote*1e3),'mm']);
hold on;
fill([-1, 1, 1, -1]*100, [-1, -1, 1, 1]*100, 'r--');
alpha(0.05);
hold off;

%% Equivalence Principle Application
% Reduced domain
X_v_reduced=X_v(abs(X_v) < 10e-2);
Y_v_reduced=Y_v(abs(Y_v) < 10e-2);

z_quote=1.5*h_sub;
field_grid = [X(:)'; Y(:)'; z_quote*ones(1,numel(X))];

E_reduced = utils.fixnan(E(:, abs(X_v) < 10e-2));
H_reduced = utils.fixnan(H(:, abs(X_v) < 10e-2));
dX = X_v_reduced(2) - X_v_reduced(1);
dY = Y_v_reduced(2) - Y_v_reduced(1);

E_x = reshape(E_reduced(1,:), length(X_v_reduced), length(Y_v_reduced));
E_y = reshape(E_reduced(2,:), length(X_v_reduced), length(Y_v_reduced));
H_x = reshape(H_reduced(1,:), length(X_v_reduced), length(Y_v_reduced));
H_y = reshape(H_reduced(2,:), length(X_v_reduced), length(Y_v_reduced));

zenith_angles = (-1: .025 :1)*pi/2;
phase_constant = 2*pi/lambda0;

eta_medium = 377;

% XZ plane
kernel_xz_plane = exp(1i*phase_constant*sin(zenith_angles.')*X_v);
P_x_xz = kernel_xz_plane*sum(E_x, 1).'*dX*dY;
P_y_xz = kernel_xz_plane*sum(E_y, 1).'*dX*dY;

Q_x_xz = kernel_xz_plane*sum(H_x, 1).'*dX*dY;
Q_y_xz = kernel_xz_plane*sum(H_y, 1).'*dX*dY;

E_theta_xz = eta_medium*cos(zenith_angles).*Q_y_xz.' + P_x_xz.';
E_phi_xz = -eta_medium*Q_x_xz.' + cos(zenith_angles).*P_y_xz.';

% YZ plane
kernel_yz_plane = exp(1j*phase_constant*sin(zenith_angles.')*Y_v);
P_x_yz = kernel_yz_plane*sum(E_x, 2)*dX*dY;
P_y_yz = kernel_yz_plane*sum(E_y, 2)*dX*dY;

Q_x_yz = kernel_yz_plane*sum(H_x, 2)*dX*dY;
Q_y_yz = kernel_yz_plane*sum(H_y, 2)*dX*dY;

E_theta_yz = -eta_medium*cos(zenith_angles).*Q_x_yz.' + P_y_yz.';
E_phi_yz = -eta_medium*Q_y_yz.' - cos(zenith_angles).*P_x_yz.';

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