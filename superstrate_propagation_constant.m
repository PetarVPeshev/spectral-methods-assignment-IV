close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');

c = physconst('LightSpeed');
wave_impedance = 376.730313668;

Ner = 102;
Nf = 101;
stratification.er = NaN(1, Ner);

%% PARAMETERS
wave.f = linspace(9, 11, Nf) * 1e9;
stratification.h = 15 * 1e-3;
stratification.er(1 : end - 1) = linspace(1, 25, Ner - 1);
stratification.er(end) = 12;
N = 1001;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
% frequency - col; er - row
stratification.hs = wave.wavelength ./ (4 * sqrt(stratification.er'));
stratification.hs(Ner, :) = 2.1 * 1e-3;
dipole.L = wave.wavelength / 2;
dipole.W = wave.wavelength / 20;

%% SPHERICAL COORDINATE SYSTEM
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, 100);
phi = linspace(0, 2 * pi, 400);
sph_grid = meshgrid_comb(theta, phi);

%% ELEVATION
z = R * cos(sph_grid(:, :, 1));

krho_te = NaN(length(stratification.er), length(wave.f));
krho_tm = NaN(length(stratification.er), length(wave.f));
dir_broadside = NaN(length(stratification.er), length(wave.f));
Es = NaN( [size(sph_grid, 1, 2), 3, 3] );
peak_dir = NaN(1, Ner);
peak_idx = NaN(1, Ner);
fh = NaN(1, Ner);
fl = NaN(1, Ner);
idx = 1;
for er_idx = 1 : 1 : length(stratification.er)
    krho_norm = linspace(1, sqrt(stratification.er(er_idx)), N);

    for f_idx = 1 : 1 : length(wave.f)
        %% TE1 AND TM1 PROPAGATION CONSTANTS
        krho = krho_norm * wave.k0(f_idx);
    
        [krho_te(er_idx, f_idx), krho_tm(er_idx, f_idx)] ...
            = find_krho(wave.k0(f_idx), krho, 'Superstrate', ...
            stratification.h, stratification.hs(er_idx, f_idx), stratification.er(er_idx));
        
        %% WAVE VECTOR COMPONENTS
        [k_comp, ~] = wave_vector(1, wave.k0(f_idx), sph_grid);
        KRHO = sqrt(k_comp(:, :, 1) .^ 2 + k_comp(:, :, 2) .^ 2);

        %% VOLTAGE AND CURRENT FIELDS OF STRATIFIED MEDIA
        [vte, ite, vtm, itm] = stratified_media(wave.k0(f_idx), KRHO, z, ...
            'Superstrate', stratification.h, stratification.hs(er_idx, f_idx), stratification.er(er_idx));
        
        %% SPECTRAL GREEN'S FUNCTIONS
        SGF = spectral_gf(1, wave.k0(f_idx), k_comp(:, :, 1), k_comp(:, :, 2), vtm, vte, itm, ite, 'E', 'M');
        
        %% MAGNETIC CURRENT DENSITY
        Mx = ft_current(wave.k0(f_idx), k_comp, dipole.W(f_idx), dipole.L(f_idx), 1, 'dipole', 'x');
        
        %% ELECTRIC FIELD
        E = farfield(wave.k0(f_idx), R, sph_grid, k_comp(:, :, 3), z, SGF, Mx);
        if er_idx == Ner - 1
            if wave.f(f_idx) == 9e9 || wave.f(f_idx) == 10e9 || wave.f(f_idx) == 11e9
                Es(:, :, :, idx) = E;
                idx = idx + 1;
            end
        end
        
        %% DIRECTIVITY
        [dir, ~, ~] = directivity(1, E, sph_grid, R);
        dir_broadside(er_idx, f_idx) = dir(1, 1);
    end

    %% RESONANCE FREQUENCY
    [peak, peak_idx_temp] = findpeaks(round(dir_broadside(er_idx, :), 4));
    if ~isempty(peak)
        peak_dir(er_idx) = peak;
        peak_idx(er_idx) = peak_idx_temp;

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
    
%% PLOT PROPAGATION CONSTANT FOR CONSTANT PERMITTIVITY
figure('Position', [250 250 750 400]);
plot(stratification.h ./ wave.wavelength, real(krho_te(end, :) ./ wave.k0), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.h ./ wave.wavelength, imag(krho_te(end, :) ./ wave.k0), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.h ./ wave.wavelength, real(krho_tm(end, :) ./ wave.k0), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.h ./ wave.wavelength, imag(krho_tm(end, :) ./ wave.k0), ...
    '--', 'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{TM1\}');
grid on;
xticks(0.45 : 0.01 : 0.55);
xlim([0.45 0.55]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{0}');
ylabel('k_{\rho} / k_{0}');
title(['Normalized k_{\rho} @ Superstrate, h = ' ...
    num2str(stratification.h * 1e3) ' mm, h_{s} = ' ...
    num2str(stratification.hs(end, 1) * 1e3) ' mm, ' ...
    'and \epsilon_{r} = ' num2str(stratification.er(end))]);
saveas(gcf, ['figures\superstrate_krho_er_const_' ...
    num2str(stratification.er(end)) '.fig']);
    
%% PLOT PROPAGATION CONSTANT
figure('Position', [250 250 750 400]);
plot(stratification.er(1 : end - 1), real(krho_te(1 : end - 1, ceil(Nf / 2)) ./ wave.k0(ceil(Nf / 2))), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.er(1 : end - 1), imag(krho_te(1 : end - 1, ceil(Nf / 2)) ./ wave.k0(ceil(Nf / 2))), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.er(1 : end - 1), real(krho_tm(1 : end - 1, ceil(Nf / 2)) ./ wave.k0(ceil(Nf / 2))), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.er(1 : end - 1), imag(krho_tm(1 : end - 1, ceil(Nf / 2)) ./ wave.k0(ceil(Nf / 2))), ...
    '--', 'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{TM1\}');
grid on;
xticks(min(stratification.er(1 : end - 1)) : 2 : max(stratification.er(1 : end - 1)));
xlim([min(stratification.er(1 : end - 1)) max(stratification.er(1 : end - 1))]);
legend show;
legend('location', 'bestoutside');
xlabel('\epsilon_{r}');
ylabel('k_{\rho} / k_{0}');
title(['Normalized k_{\rho} @ Superstrate, f = ' ...
    num2str(wave.f(ceil(Nf / 2)) * 1e-9) ' GHz, h = ' ...
    num2str(stratification.h * 1e3) ' mm, ' ...
    'h_{s} = \lambda_{0} / (4 * sqrt(\epsilon_{r}))']);
saveas(gcf, ['figures\superstrate_krho_f_const_' ...
    num2str(wave.f(ceil(Nf / 2)) * 1e-9) 'GHz.fig']);

%% PLOT DIRECTIVITY
er_idx_7 = find(stratification.er == 7, 1);
er_idx_12 = find(stratification.er == 12, 1);
er_idx_25 = find(stratification.er == 25, 1);
figure('Position', [250 250 750 400]);
plot(wave.f * 1e-9, 10 * log10(dir_broadside(er_idx_7, :)), ...
    'LineWidth', 2.0, 'DisplayName', ...
    ['dir, \epsilon_{r} = ' num2str(stratification.er(er_idx_7))]);
hold on;
plot(wave.f * 1e-9, 10 * log10(dir_broadside(er_idx_12, :)), ...
    'LineWidth', 2.0, 'DisplayName', ...
    ['dir, \epsilon_{r} = ' num2str(stratification.er(er_idx_12))]);
hold on;
plot(wave.f * 1e-9, 10 * log10(dir_broadside(er_idx_25, :)), ...
    'LineWidth', 2.0, 'DisplayName', ...
    ['dir, \epsilon_{r} = ' num2str(stratification.er(er_idx_25))]);
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('f / GHz');
ylabel('D(\theta=0,\phi=0) / dB');
title(['Broadside Directivity @ Superstrate, h = ' ...
    num2str(stratification.h * 1e3) ' mm, ' ...
    'h_{s} = \lambda_{0} / (4 * sqrt(\epsilon_{r}))']);
saveas(gcf, 'figures\superstrate_dir.fig');

%% PLOT ELECTRIC FAR-FIELD
% Normalized electric far-field
Enorm = NaN( [size(sph_grid, 1, 2), size(Es, 4)] );
for E_idx = 1 : 1 : size(Es, 4)
    E_tot = total_field(Es(:, :, :, E_idx));
    Enorm(:, :, E_idx) = norm_magnitude(E_tot, 'dB');
end
% Elevation angle
theta_plot = NaN(1, 2 * length(theta));
theta_plot(1 : length(theta)) =  - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
freq = [9 10 11];
color_styles = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
figure('Position', [250 250 750 400]);
for E_idx = 1 : 1 : size(Enorm, 3)
    % E-plane
    plane_idx_1 = find(round(phi * 180 / pi, 0) == 0, 1);
    plane_idx_2 = find(round(phi * 180 / pi, 0) == 180, 1);
    phi0_plot = NaN(1, 2 * length(theta));
    phi0_plot(1 : length(theta)) = fliplr(Enorm(plane_idx_2, :, E_idx));
    phi0_plot(length(theta) + 1 : end) = Enorm(plane_idx_1, :, E_idx);
    plot(theta_plot, phi0_plot, 'Color', color_styles(E_idx, :), ...
        'LineWidth', 2.0, 'DisplayName', ['f = ' num2str(freq(E_idx)) ...
        ' GHz, E-plane'])
    hold on;
    % H-plane
    plane_idx_1 = find(round(phi * 180 / pi, 0) == 90, 1);
    plane_idx_2 = find(round(phi * 180 / pi, 0) == 270, 1);
    phi90_plot = NaN(1, 2 * length(theta));
    phi90_plot(1 : length(theta)) = fliplr(Enorm(plane_idx_2, :, E_idx));
    phi90_plot(length(theta) + 1 : end) = Enorm(plane_idx_1, :, E_idx);
    plot(theta_plot, phi90_plot, '--', 'Color', color_styles(E_idx, :), ...
        'LineWidth', 2.0, 'DisplayName', ['f = ' num2str(freq(E_idx)) ...
        ' GHz, H-plane'])
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
title(['E Far-Field @ Superstrate, \epsilon_{r} = ' ...
    num2str(stratification.er(Ner - 1))]);
saveas(gcf, 'figures\superstrate_Eff.fig');

%% PLOT BANDWIDTH
figure('Position', [250 250 750 400]);
plot(stratification.er(1 : end - 1), BW(1 : end - 1), 'LineWidth', 2.0, ...
    'DisplayName', 'BW');
grid on;
ylim([0 13.5]);
legend show;
legend('location', 'bestoutside');
xlabel('\epsilon_{r}');
ylabel('BW / %');
title('Bandwidth @ Superstrate');
saveas(gcf, 'figures\superstrate_bw.fig');
