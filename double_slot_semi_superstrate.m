close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');

c = physconst('LightSpeed');
wave_impedance = 376.730313668;

Nf = 101;

%% PARAMETERS
wave.f = 15 * 1e9;
stratification.h = 10 * 1e-3;
stratification.er = 10;
N = 1001;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
double_slot.L = wave.wavelength / 2;
double_slot.W = wave.wavelength / 20;
single_slot.L = wave.wavelength / 2;
single_slot.W = wave.wavelength / 20;

%% TM0 PROPAGATION CONSTANT
krho_tm0 = find_krho_tm0(wave.k0, 'SemiInfiniteSuperstrate', ...
    stratification.h, stratification.er);

%% OPTIMUM DISTANCE
% kx_lw = krho_lw
% k
% wavelength_tm0 = 
double_slot.d = pi / real(krho_tm0);

%% SPHERICAL COORDINATE SYSTEM
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, 100);
phi = linspace(0, 2 * pi, 400);
sph_grid = meshgrid_comb(theta, phi);

%% ELEVATION
z = R * cos(sph_grid(:, :, 1));

%% WAVE VECTOR
[k_comp, k] = wave_vector(stratification.er, wave.k0, sph_grid);
KRHO = sqrt(k_comp(:, :, 1) .^ 2 + k_comp(:, :, 2) .^ 2);
k_comp(:, :, 3) = 1j * sqrt(- k ^ 2 + KRHO .^ 2);
    
%% VOLTAGE AND CURRENT FIELDS OF STRATIFIED MEDIA
[vte, ite, vtm, itm] = stratified_media(wave.k0, KRHO, z, ...
    'SemiInfiniteSuperstrate', stratification.h, stratification.er);
        
%% SPECTRAL GREEN'S FUNCTIONS
SGF = spectral_gf(stratification.er, k, k_comp(:, :, 1), ...
    k_comp(:, :, 2), vtm, vte, itm, ite, 'E', 'M');

%% SINGLE SLOT MAGNETIC CURRENT
single_slot.M = ft_current(wave.k0, k_comp, single_slot.W, ...
    single_slot.L, 1, 'dipole', 'x');

%% DOUBLE SLOT MAGNETIC CURRENT
double_slot.M = single_slot.M ...
    .* 2 .* cos(k_comp(:, :, 1) * double_slot.d / 2);

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
% E-plane
e_plane_idx_1 = find(round(phi * 180 / pi, 0) == 0, 1);
e_plane_idx_2 = find(round(phi * 180 / pi, 0) == 180, 1);
ss_e_plane = NaN(1, 2 * length(theta));
ss_e_plane(1 : length(theta)) = fliplr(ss_Enorm(e_plane_idx_2, :));
ss_e_plane(length(theta) + 1 : end) = ss_Enorm(e_plane_idx_1, :);
ds_e_plane = NaN(1, 2 * length(theta));
ds_e_plane(1 : length(theta)) = fliplr(ds_Enorm(e_plane_idx_2, :));
ds_e_plane(length(theta) + 1 : end) = ds_Enorm(e_plane_idx_1, :);
figure('Position', [250 100 650 650]);
subplot(2, 1, 1);
plot(theta_plot, ss_e_plane, 'LineWidth', 2.0, ...
    'DisplayName', 'single-slot')
hold on;
plot(theta_plot, ds_e_plane, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'double-slot')
grid on;
ylim([-40 0]);
% xlim([-40 40]);
% xticks(-40 : 8 : 40);
legend show;
legend('location', 'bestoutside');
ylabel('|E| / dB');
title('E-plane');
% H-plane
h_plane_idx_1 = find(round(phi * 180 / pi, 0) == 90, 1);
h_plane_idx_2 = find(round(phi * 180 / pi, 0) == 270, 1);
ss_h_plane = NaN(1, 2 * length(theta));
ss_h_plane(1 : length(theta)) = fliplr(ss_Enorm(h_plane_idx_2, :));
ss_h_plane(length(theta) + 1 : end) = ss_Enorm(h_plane_idx_1, :);
ds_h_plane = NaN(1, 2 * length(theta));
ds_h_plane(1 : length(theta)) = fliplr(ds_Enorm(h_plane_idx_2, :));
ds_h_plane(length(theta) + 1 : end) = ds_Enorm(h_plane_idx_1, :);
subplot(2, 1, 2);
plot(theta_plot, ss_h_plane, 'LineWidth', 2.0, ...
    'DisplayName', 'single-slot')
hold on;
plot(theta_plot, ds_h_plane, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'double-slot')
grid on;
ylim([-40 0]);
% xlim([-40 40]);
% xticks(-40 : 8 : 40);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / dB');
title('H-plane');
sgtitle(['|E^{FF}| @ Semi-Infinite Superstrate, \epsilon_{r} = ' ...
    num2str(stratification.er) ', and h = ' ...
    num2str(stratification.h * 1e3) ' mm'], 'FontWeight', 'bold', ...
    'FontSize', 12);
saveas(gcf, 'figures\double_slot_eff.fig')

%% PRINT DOUBLE-SLOT ANTENNA DIRECTIVITY
fprintf('Double-slot antenna directivity: %.2f dB\n', ...
    10 * log10(double_slot.dir_broadside));
