close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');

c = physconst('LightSpeed');
wave_impedance = 376.730313668;

Ner = 101;
Nf = 101;

%% PARAMETERS
wave.f = linspace(10, 20, Nf) * 1e9;
stratification.h = 10 * 1e-3;
stratification.er = linspace(1, 25, Ner);
N = 1001;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
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
for er_idx = 1 : 1 : length(stratification.er)
    krho_norm = linspace(1, sqrt(stratification.er(er_idx)), N);

    for f_idx = 1 : 1 : length(wave.f)
        %% TE1 AND TM1 PROPAGATION CONSTANTS
        krho = krho_norm * wave.k0(f_idx);

        [krho_te(er_idx, f_idx), krho_tm(er_idx, f_idx)] ...
            = find_krho(wave.k0(f_idx), krho, 'SemiInfiniteSuperstrate', ...
            stratification.h, stratification.er(er_idx));
    end
end
    
%% PLOT PROPAGATION CONSTANT FOR CONSTANT PERMITTIVITY
figure('Position', [250 250 750 400]);
plot(stratification.h * sqrt(stratification.er(end)) ./ wave.wavelength, real(krho_te(end, :) ./ (wave.k0 * sqrt(stratification.er(end)))), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.h * sqrt(stratification.er(end)) ./ wave.wavelength, imag(krho_te(end, :) ./ (wave.k0 * sqrt(stratification.er(end)))), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.h * sqrt(stratification.er(end)) ./ wave.wavelength, real(krho_tm(end, :) ./ (wave.k0 * sqrt(stratification.er(end)))), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.h * sqrt(stratification.er(end)) ./ wave.wavelength, imag(krho_tm(end, :) ./ (wave.k0 * sqrt(stratification.er(end)))), ...
    '--', 'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{TM1\}');
grid on;
xticks(1.65 : 0.15 : 3.3);
xlim([1.65 3.3]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('k_{\rho} / k_{d}');
title(['Normalized k_{\rho} @ Semi-Infinite Superstrate, h = ' ...
    num2str(stratification.h * 1e3) ' mm, ' ...
    'and \epsilon_{r} = ' num2str(stratification.er(end))]);
saveas(gcf, ['figures\semi_superstrate_krho_er_const_' ...
    num2str(stratification.er(end)) '.fig']);
    
%% PLOT PROPAGATION CONSTANT
figure('Position', [250 250 750 400]);
plot(stratification.er, real(krho_te(:, ceil(Nf / 2)) ./ (wave.k0(ceil(Nf / 2)) .* sqrt(stratification.er'))), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.er, imag(krho_te(:, ceil(Nf / 2)) ./ (wave.k0(ceil(Nf / 2)) .* sqrt(stratification.er'))), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.er, real(krho_tm(:, ceil(Nf / 2)) ./ (wave.k0(ceil(Nf / 2)) .* sqrt(stratification.er'))), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.er, imag(krho_tm(:, ceil(Nf / 2)) ./ (wave.k0(ceil(Nf / 2)) .* sqrt(stratification.er'))), ...
    '--', 'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{TM1\}');
grid on;
xticks(min(stratification.er) : 2 : max(stratification.er));
xlim([min(stratification.er) max(stratification.er)]);
legend show;
legend('location', 'bestoutside');
xlabel('\epsilon_{r}');
ylabel('k_{\rho} / k_{0}');
title(['Normalized k_{\rho} @ Semi-Infinite Superstrate, f = ' ...
    num2str(wave.f(ceil(Nf / 2)) * 1e-9) ' GHz, and h = ' ...
    num2str(stratification.h * 1e3) ' mm']);
saveas(gcf, ['figures\semi_superstrate_krho_f_const_' ...
    num2str(wave.f(ceil(Nf / 2)) * 1e-9) 'GHz.fig']);
