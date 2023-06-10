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
wave.f = linspace(8, 12, 201) * 1e9;
% Stratification parameters
stratification.h = 15 * 1e-3;
stratification.er = linspace(1, 25, 121);
% Coordinate system parameters
N = 1001;
R = 1;

%% DEPENDENT PARAMETERS
% Wave dependent parameters
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
% Stratification dependent parameters
stratification.hs = wave.wavelength ./ (4 * sqrt(stratification.er'));
% Dipole dependent parameters
dipole.L = wave.wavelength / 2;
dipole.W = wave.wavelength / 20;

%% SPHERICAL COORDINATE SYSTEM
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, 100);
phi = linspace(0, 2 * pi, 400);
sph_grid = meshgrid_comb(theta, phi);

%% ELEVATION
z = R * cos(sph_grid(:, :, 1));

%% WAVE VECTOR COMPONENTS
[k_comp, ~, kx, ky, kz] = wave_vector_multi_freq(1, wave.k0, sph_grid);
krho = sqrt(kx .^ 2 + ky .^ 2);
        
%% MAGNETIC CURRENT DENSITY
Mx = ft_current_multi_freq(wave.k0, k_comp, dipole.W, dipole.L, 1, ...
    'dipole', 'x');

Es = NaN( [size(sph_grid, 1, 2), 3] );
dir_broadside = NaN(length(stratification.er), length(wave.f));
fh = NaN(1, length(stratification.er));
fl = NaN(1, length(stratification.er));
for er_idx = 1 : 1 : length(stratification.er)
    %% VOLTAGE AND CURRENT FIELDS OF STRATIFIED MEDIA
    [vte, ite, vtm, itm] ...
        = stratified_media_multi_freq(wave.k0, krho, z, ...
        'Superstrate', stratification.h, stratification.hs(er_idx, :), ...
        stratification.er(er_idx));
            
    %% SPECTRAL GREEN'S FUNCTIONS
    SGF = spectral_gf_multi_freq(1, wave.k0, kx, ky, vtm, vte, itm, ite, ...
        'E', 'M');
            
    %% ELECTRIC FIELD
    E = farfield_multi_freq(wave.k0, R, sph_grid, kz, z, SGF, Mx);
    if er_idx == length(stratification.er)
        Es(:, :, 1) = sqrt( abs(E(:, :, 1, 1)) .^ 2 ...
            + abs(E(:, :, 2, 1)) .^ 2 + abs(E(:, :, 3, 1)) .^ 2);
        Es(:, :, 2) = sqrt( abs(E(:, :, 1, 101)) .^ 2 ...
            + abs(E(:, :, 2, 101)) .^ 2 + abs(E(:, :, 3, 101)) .^ 2);
        Es(:, :, 3) = sqrt( abs(E(:, :, 1, end)) .^ 2 ...
            + abs(E(:, :, 2, end)) .^ 2 + abs(E(:, :, 3, end)) .^ 2);
    end
            
    %% DIRECTIVITY
    dir_broadside(er_idx, :) = broadside_directivity(1, E, sph_grid, R);

    %% -3 dB POINTS
    [peak, peak_idx] = max(dir_broadside(er_idx, :));
    if ~isempty(peak)
        % Low frequency point
        fl_temp ...
            = find(dir_broadside(er_idx, 1 : peak_idx) < (peak / 2), ...
            1, 'last');
        if ~isempty(fl_temp)
            fl(er_idx) = wave.f(fl_temp);
        end
        % High frequency point
        fh_temp ...
            = find(dir_broadside(er_idx, peak_idx : end) < (peak / 2), ...
            1) + peak_idx - 1;
        if ~isempty(fh_temp)
            fh(er_idx) = wave.f(fh_temp);
        end
    end
end

%% BANDWIDTH
BW = 200 * (fh - fl) ./ (fh + fl);

%% PLOT DIRECTIVITY
figure('Position', [250 250 750 400]);
for idx = 1 : 20 : length(stratification.er)
    plot(wave.f * 1e-9, 10 * log10(dir_broadside(idx, :)), ...
        'LineWidth', 2.0, 'DisplayName', ...
        ['dir, \epsilon_{r} = ' num2str(stratification.er(idx))]);
    hold on;
end
hold off;
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('f / GHz');
ylabel('D(\theta=0,\phi=0) / dB');
title(['Broadside Directivity @ Superstrate, ' ...
    num2str(stratification.h * 1e3) ' mm, and ' ...
    'h_{s} = \lambda_{0} / (4 * sqrt(\epsilon_{r}))']);
saveas(gcf, 'figures\superstrate_dir.fig');

%% PLOT ELECTRIC FAR-FIELD
% Normalized electric far-field
Enorm = NaN( [size(Es, 1, 2, 3)] );
for e_idx = 1 : 1 : size(Es, 3)
    Enorm(:, :, e_idx) = norm_magnitude(Es(:, :, e_idx), 'dB');
end
% Elevation angle
theta_plot = NaN(1, 2 * length(theta));
theta_plot(1 : length(theta)) =  - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
freq = [8 10 12];
% Plot at phi in 0 deg
figure('Position', [250 250 650 650]);
subplot(2, 1, 1);
for e_idx = 1 : 1 : size(Enorm, 3)
    plane_idx_1 = find(round(phi * 180 / pi, 0) == 0, 1);
    plane_idx_2 = find(round(phi * 180 / pi, 0) == 180, 1);
    phi_plot = NaN(1, 2 * length(theta));
    phi_plot(1 : length(theta)) = fliplr(Enorm(plane_idx_1, :, e_idx));
    phi_plot(length(theta) + 1 : end) = Enorm(plane_idx_2, :, e_idx);
    plot(theta_plot, phi_plot, 'LineWidth', 2.0, 'DisplayName', ...
        ['f = ' num2str(freq(e_idx)) ' GHz'])
    hold on;
end
hold off;
grid on;
ylim([-40 0]);
xlim([-90 90]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / dB');
title('\phi = 0 deg');
% title(['E Far-Field @ Superstrate, \phi = 0 deg, \epsilon_{r} = ' ...
%     num2str(stratification.er(end))]);
% saveas(gcf, 'figures\superstrate_Eff_phi0.fig');
% Plot at phi in 90 deg
% figure('Position', [250 250 750 400]);
subplot(2, 1, 2);
for e_idx = 1 : 1 : size(Enorm, 3)
    plane_idx_1 = find(round(phi * 180 / pi, 0) == 90, 1);
    plane_idx_2 = find(round(phi * 180 / pi, 0) == 270, 1);
    phi_plot = NaN(1, 2 * length(theta));
    phi_plot(1 : length(theta)) = fliplr(Enorm(plane_idx_2, :, e_idx));
    phi_plot(length(theta) + 1 : end) = Enorm(plane_idx_1, :, e_idx);
    plot(theta_plot, phi_plot, 'LineWidth', 2.0, 'DisplayName', ...
        ['f = ' num2str(freq(e_idx)) ' GHz'])
    hold on;
end
hold off;
grid on;
ylim([-40 0]);
xlim([-90 90]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E| / dB');
title('\phi = 90 deg');
sgtitle(['E Far-Field @ Superstrate, \epsilon_{r} = ' ...
    num2str(stratification.er(end))], 'FontWeight', 'bold', ...
    'FontSize', 12);
saveas(gcf, 'figures\superstrate_Eff.fig');

%% PLOT BANDWIDTH
figure('Position', [250 250 750 400]);
plot(stratification.er, BW, 'LineWidth', 2.0, 'DisplayName', 'BW');
grid on;
ylim([0 29]);
legend show;
legend('location', 'bestoutside');
xlabel('\epsilon_{r}');
ylabel('BW / %');
title(['Bandwidth @ Superstrate, h = ' num2str(stratification.h * 1e3) ...
    ' mm']);
saveas(gcf, 'figures\superstrate_bw.fig');
