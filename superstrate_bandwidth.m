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
wave.f = linspace(8, 12, 101) * 1e9;
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
