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
wave.f = 10 * 1e9;
stratification.h = 15 * 1e-3;
stratification.er = linspace(1, 25, 1001);
N = 1001;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
stratification.hs = wave.wavelength ./ (4 * sqrt(stratification.er));
dipole.L = wave.wavelength / 2;
dipole.W = wave.wavelength / 20;

krho_te = NaN(1, length(stratification.er));
krho_tm = NaN(1, length(stratification.er));
for er_idx = 1 : 1 : length(stratification.er)
    krho = linspace(1, sqrt(stratification.er(er_idx)), N) * wave.k0;
    
    [krho_te(er_idx), krho_tm(er_idx)] = find_krho(wave.k0, krho, 'Superstrate', ...
        stratification.h, stratification.hs(er_idx), stratification.er(er_idx));
end
    
%% PLOT PROPAGATION CONSTANT
figure('Position', [250 250 750 400]);
plot(stratification.er, real(krho_te / wave.k0), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.er, imag(krho_te / wave.k0), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.er, real(krho_tm / wave.k0), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.er, imag(krho_tm / wave.k0), ...
    '--', 'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{TM1\}');
grid on;
xlim([min(stratification.er) max(stratification.er)]);
legend show;
xticks(min(stratification.er) : 2 : max(stratification.er));
legend('location', 'bestoutside');
xlabel('\epsilon_{r}');
ylabel('k_{\rho} / k_{0}');
title(['Normalized k_{\rho} @ Superstrate, h = ' ...
    num2str(stratification.h * 1e3) ' mm, ' ...
    'h_{s} = \lambda_{0} / (4 * sqrt(\epsilon_{r}))']);
saveas(gcf, ['figures\propagation_constant_f_' num2str(wave.f * 1e-9) ...
    'GHz.fig']);
