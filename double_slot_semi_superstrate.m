close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');

c = physconst('LightSpeed');
wave_impedance = 376.730313668;

%% PARAMETERS
% Wave parameters
wave.f = 15 * 1e9;
% Stratification parameters
stratification.h = 10 * 1e-3;
stratification.er = 10;
% Coordinate system parameters
N = 1001;
R = 1;

%% DEPENDENT PARAMETERS
% Wave dependent parameters
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
% Double slot dependent parameters
double_slot.L = wave.wavelength / 2;
double_slot.W = wave.wavelength / 20;
% Single slot dependent parameters
single_slot.L = wave.wavelength / 2;
single_slot.W = wave.wavelength / 20;

%% TM0 PROPAGATION CONSTANT
krho_tm0 = find_krho_tm0(wave.k0, 'SemiInfiniteSuperstrate', ...
    stratification.h, stratification.er);

%% OPTIMUM DISTANCE FOR DOUBLE SLOT
double_slot.d = pi / real(krho_tm0);

%% SPHERICAL COORDINATE SYSTEM
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, 100);
phi = linspace(0, 2 * pi, 400);
sph_grid = meshgrid_comb(theta, phi);

%% ELEVATION
z = R * cos(sph_grid(:, :, 1));

%% WAVE VECTOR
[k_comp, ~, kx, ky] ...
    = wave_vector_multi_freq(stratification.er, wave.k0, sph_grid);
% [k_comp, k] = wave_vector(stratification.er, wave.k0, sph_grid);
k = wave.k0 * sqrt(stratification.er);
KRHO = sqrt(kx .^ 2 + ky .^ 2);
k_comp(:, :, 3) = - 1j * sqrt(- k ^ 2 + KRHO .^ 2);
% k_comp(:, :, 3) = 1j * sqrt(- k ^ 2 + KRHO .^ 2);
    
%% VOLTAGE AND CURRENT FIELDS OF STRATIFIED MEDIA
[vte, ite, vtm, itm] = stratified_media(wave.k0, KRHO, z, ...
    'SemiInfiniteSuperstrate', stratification.h, stratification.er);
        
%% SPECTRAL GREEN'S FUNCTIONS
SGF = spectral_gf(stratification.er, k, kx, ky, vtm, vte, itm, ite, ...
    'E', 'M');

%% SINGLE SLOT MAGNETIC CURRENT
single_slot.M = ft_current(wave.k0, k_comp, single_slot.W, ...
    single_slot.L, 1, 'dipole', 'x');

%% DOUBLE SLOT MAGNETIC CURRENT
double_slot.M = single_slot.M .* 2 .* cos(kx * double_slot.d / 2);

%% SINGLE SLOT ELECTRIC FAR-FIELD
single_slot.E = farfield(k, R, sph_grid, k_comp(:, :, 3), z, SGF, ...
    single_slot.M);
single_slot.Etotal = total_field(single_slot.E);

%% DOUBLE SLOT ELECTRIC FAR-FIELD
double_slot.E = farfield(k, R, sph_grid, k_comp(:, :, 3), z, SGF, ...
    double_slot.M);
double_slot.Etotal = total_field(double_slot.E);

%% SINGLE SLOT DIRECTIVITY
[single_slot.dir, ~, ~] = directivity(stratification.er, single_slot.E, ...
    sph_grid, R);
single_slot.dir_broadside = single_slot.dir(1, 1);

%% DOUBLE SLOT DIRECTIVITY
[double_slot.dir, ~, ~] = directivity(stratification.er, double_slot.E, ...
    sph_grid, R);
double_slot.dir_broadside = double_slot.dir(1, 1);

%% PLOT ELECTRIC FAR-FIELD
% Elevation angle
theta_plot = NaN(1, 2 * length(theta));
theta_plot(1 : length(theta)) =  - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
% Normalization
ss_Enorm = norm_magnitude(single_slot.Etotal, 'dB');
ds_Enorm = norm_magnitude(double_slot.Etotal, 'dB');
% Plane indecies
plane_0 = find(round(phi * 180 / pi, 0) == 0, 1);
plane_45 = find(round(phi * 180 / pi, 0) == 45, 1);
plane_90 = find(round(phi * 180 / pi, 0) == 90, 1);
plane_180 = find(round(phi * 180 / pi, 0) == 180, 1);
plane_225 = find(round(phi * 180 / pi, 0) == 225, 1);
plane_270 = find(round(phi * 180 / pi, 0) == 270, 1);
% Single slot antenna
ss_0_plane = NaN(1, 2 * length(theta));
ss_0_plane(1 : length(theta)) = fliplr(ss_Enorm(plane_180, :));
ss_0_plane(length(theta) + 1 : end) = ss_Enorm(plane_0, :);
ss_45_plane = NaN(1, 2 * length(theta));
ss_45_plane(1 : length(theta)) = fliplr(ss_Enorm(plane_225, :));
ss_45_plane(length(theta) + 1 : end) = ss_Enorm(plane_45, :);
ss_90_plane = NaN(1, 2 * length(theta));
ss_90_plane(1 : length(theta)) = fliplr(ss_Enorm(plane_270, :));
ss_90_plane(length(theta) + 1 : end) = ss_Enorm(plane_90, :);
% Double slot antenna
ds_0_plane = NaN(1, 2 * length(theta));
ds_0_plane(1 : length(theta)) = fliplr(ds_Enorm(plane_180, :));
ds_0_plane(length(theta) + 1 : end) = ds_Enorm(plane_0, :);
ds_45_plane = NaN(1, 2 * length(theta));
ds_45_plane(1 : length(theta)) = fliplr(ds_Enorm(plane_225, :));
ds_45_plane(length(theta) + 1 : end) = ds_Enorm(plane_45, :);
ds_90_plane = NaN(1, 2 * length(theta));
ds_90_plane(1 : length(theta)) = fliplr(ds_Enorm(plane_270, :));
ds_90_plane(length(theta) + 1 : end) = ds_Enorm(plane_90, :);
% Plot single slot antenna
figure('Position', [250 250 750 400]);
plot(theta_plot, ss_0_plane, 'LineWidth', 2.0, ...
    'DisplayName', '\phi = 0 deg');
hold on;
plot(theta_plot, ss_45_plane, 'LineWidth', 2.0, ...
    'DisplayName', 'phi = 45 deg');
hold on;
plot(theta_plot, ss_90_plane, 'LineWidth', 2.0, ...
    'DisplayName', 'phi = 90 deg');
grid on;
xticks(-30 : 5 : 30);
ylim([-40 0]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / dB');
title(['|E^{FF}| @ Semi-Infinite Superstrate, Single-Slot Antenna, ' ...
    'h = ' num2str(stratification.h * 1e3) ' mm, and ' ...
    '\epsilon_{r} = ' num2str(stratification.er) '']);
saveas(gcf, 'figures\single_slot_eff.fig');
% H-plane
% Plot double slot antenna
figure('Position', [250 250 750 400]);
plot(theta_plot, ds_0_plane, 'LineWidth', 2.0, ...
    'DisplayName', '\phi = 0 deg');
hold on;
plot(theta_plot, ds_45_plane, 'LineWidth', 2.0, ...
    'DisplayName', '\phi = 45 deg');
hold on;
plot(theta_plot, ds_90_plane, 'LineWidth', 2.0, ...
    'DisplayName', '\phi = 90 deg');
grid on;
xticks(-30 : 5 : 30);
ylim([-40 0]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / dB');
title(['|E^{FF}| @ Semi-Infinite Superstrate, Double-Slor Antenna, ' ...
    'h = ' num2str(stratification.h * 1e3) ' mm, and ' ...
    '\epsilon_{r} = ' num2str(stratification.er)]);
saveas(gcf, 'figures\double_slot_eff.fig');

%% PRINT DOUBLE-SLOT ANTENNA DIRECTIVITY
fprintf('Double-slot antenna directivity: %.2f dB\n', ...
    10 * log10(double_slot.dir_broadside));
