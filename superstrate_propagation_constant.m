close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');
c = physconst('LightSpeed');
wave_impedance = 376.730313668;

Ner = 1002;
Nf = 1001;
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
E = NaN( [size(sph_grid, 1, 2), 3, 3] );
idx = 0;
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
                E(:, :, :, idx) = E;
                idx = idx + 1;
            end
        end
        
        %% DIRECTIVITY
        [dir, ~, ~] = directivity(1, E, sph_grid, R);
        dir_broadside(er_idx, f_idx) = dir(1, 1);
    end
end
    
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
    'and \epsilon_{r} = ' num2str(stratification.er(102))]);
saveas(gcf, ['figures\superstrate_krho_er_const_' ...
    num2str(stratification.er(end)) '.fig']);
    
%% PLOT PROPAGATION CONSTANT
figure('Position', [250 250 750 400]);
plot(stratification.er(1 : 101), real(krho_te(1 : end - 1, ceil(Nf / 2)) ./ wave.k0(ceil(Nf / 2))), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.er(1 : 101), imag(krho_te(1 : end - 1, ceil(Nf / 2)) ./ wave.k0(ceil(Nf / 2))), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.er(1 : 101), real(krho_tm(1 : end - 1, ceil(Nf / 2)) ./ wave.k0(ceil(Nf / 2))), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.er(1 : 101), imag(krho_tm(1 : end - 1, ceil(Nf / 2)) ./ wave.k0(ceil(Nf / 2))), ...
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
