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
wave.f = linspace(9, 11, 101) * 1e9;
stratification.h = 15 * 1e-3;
stratification.er = 1 : 0.5 : 25;
N = 1001;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
stratification.hs = wave.wavelength ./ (4 * sqrt(stratification.er'));
dipole.L = wave.wavelength / 2;
dipole.W = wave.wavelength / 20;

%% SPHERICAL COORDINATE SYSTEM
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, 100);
phi = linspace(0, 2 * pi, 400);
sph_grid = meshgrid_comb(theta, phi);

%% ELEVATION
z = R * cos(sph_grid(:, :, 1));

%% WAVE VECTOR COMPONENTS
kx = permute(wave.k0, [3 1 2]) .* sin(sph_grid(:, :, 1)) .* cos(sph_grid(:, :, 2));
ky = permute(wave.k0, [3 1 2]) .* sin(sph_grid(:, :, 1)) .* sin(sph_grid(:, :, 2));
krho = sqrt(kx .^ 2 + ky .^ 2);
kz = - 1j * sqrt(- permute(repmat(wave.k0', 1, 400, 100), [2 3 1]) .^ 2 + krho .^ 2);

%% DIRECTIVITY AND BANDWIDTH
krho_te = NaN(length(stratification.er), length(wave.f));
krho_tm = NaN(length(stratification.er), length(wave.f));
dir_broadside = NaN(length(stratification.er), length(wave.f));
fh = NaN(1, length(stratification.er));
fl = NaN(1, length(stratification.er));
for er_idx = 1 : 1 : length(stratification.er)
    for f_idx = 1 : 1 : length(wave.f)
        %% VOLTAGE AND CURRENT FIELDS OF STRATIFIED MEDIA
        [vte, ite, vtm, itm] = stratified_media(wave.k0(f_idx), krho(:, :, f_idx), z, ...
            'Superstrate', stratification.h, stratification.hs(er_idx, f_idx), stratification.er(er_idx));
        
        %% SPECTRAL GREEN'S FUNCTIONS
        SGF = spectral_gf(1, wave.k0(f_idx), kx(:, :, f_idx), ky(:, :, f_idx), vtm, vte, itm, ite, 'E', 'M');
        
        %% MAGNETIC CURRENT DENSITY
        k_comp = NaN( [size(sph_grid, 1, 2), 3] );
        k_comp(:, :, 1) = kx(:, :, f_idx);
        k_comp(:, :, 2) = ky(:, :, f_idx);
        Mx = ft_current(wave.k0(f_idx), k_comp, dipole.W(f_idx), dipole.L(f_idx), 1, 'dipole', 'x');
        
        %% ELECTRIC FIELD
        E = farfield(wave.k0(f_idx), R, sph_grid, kz(:, :, f_idx), z, SGF, Mx);
        
        %% DIRECTIVITY
        [dir, ~, ~] = directivity(1, E, sph_grid, R);
        dir_broadside(er_idx, f_idx) = dir(1, 1);
    end

    %% RESONANCE FREQUENCY
    [peak, peak_idx_temp] = findpeaks(round(dir_broadside(er_idx, :), 4));
    if ~isempty(peak)
        %% LOW AND HIGH FREQUENCY POINTS
        % Low frequency cut
        fl_temp = find(dir_broadside(er_idx, 1 : peak_idx_temp) ...
            < peak - 3, 1, 'last');
        if ~isempty(fl_temp)
            fl(er_idx) = wave.f(fl_temp);
        end
        % High frequency cut
        fh_temp = find(dir_broadside(er_idx, peak_idx_temp : end) ...
            < peak - 3, 1) + peak_idx_temp - 1;
        if ~isempty(fh_temp)
            fh(er_idx) = wave.f(fh_temp);
        end
    end
end

%% BANDWIDTH
BW = 200 * (fh - fl) ./ (fh + fl);

%% PLOT DIRECTIVITY
figure('Position', [250 250 750 400]);
for idx = 1 : 2 : length(stratification.er)
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
title(['Broadside Directivity @ Superstrate, h = ' ...
    num2str(stratification.h * 1e3) ' mm, ' ...
    'h_{s} = \lambda_{0} / (4 * sqrt(\epsilon_{r}))']);
saveas(gcf, 'figures\superstrate_dir_dedicated.fig');

%% PLOT BANDWIDTH
figure('Position', [250 250 750 400]);
plot(stratification.er, BW, 'LineWidth', 2.0, 'DisplayName', 'BW');
grid on;
% ylim([0 13.5]);
legend show;
legend('location', 'bestoutside');
xlabel('\epsilon_{r}');
ylabel('BW / %');
title('Bandwidth @ Superstrate');
saveas(gcf, 'figures\superstrate_bw_dedicated.fig');
