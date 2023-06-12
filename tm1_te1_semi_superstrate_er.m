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
wave.f = 15 * 1e9;
stratification.h = 10 * 1e-3;
stratification.er = linspace(1, 25, 1001);
N = 1001;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;

krho_te = NaN(1, length(stratification.er));
krho_tm = NaN(1, length(stratification.er));
for er_idx = 1 : 1 : length(stratification.er)
    krho = linspace(1, sqrt(stratification.er(er_idx)), N) * wave.k0;

    [krho_te(er_idx), krho_tm(er_idx)] = find_krho(wave.k0, krho, ...
        'SemiInfiniteSuperstrate', stratification.h, stratification.er(er_idx));
end
    
%% PLOT PROPAGATION CONSTANT
krho_te_norm = krho_te ./ (wave.k0 * sqrt(stratification.er));
krho_tm_norm = krho_tm ./ (wave.k0 * sqrt(stratification.er));
figure('Position', [250 250 750 400]);
plot(stratification.er, real(krho_te_norm), 'LineWidth', 2.0, ...
    'Color', [0 0.4470 0.7410], 'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.er, imag(krho_te_norm), '--', 'LineWidth', 2.0, ...
    'Color', [0 0.4470 0.7410], 'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.er, real(krho_tm_norm), 'LineWidth', 2.0, ...
    'Color', [0.8500 0.3250 0.0980], 'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.er, imag(krho_tm_norm), '--', 'LineWidth', 2.0, ...
    'Color', [0.8500 0.3250 0.0980], 'DisplayName', '\Im\{TM1\}');
grid on;
xticks(min(stratification.er) : 2 : max(stratification.er));
xlim([min(stratification.er) max(stratification.er)]);
legend show;
legend('location', 'bestoutside');
xlabel('\epsilon_{r}');
ylabel('k_{\rho} / k_{d}');
title(['Normalized k_{\rho} @ Semi-Infinite Superstrate, f = ' ...
    num2str(wave.f * 1e-9) ' GHz, and h = ' ...
    num2str(stratification.h * 1e3) ' mm']);
saveas(gcf, ['figures\semi_superstrate_te1_tm1_f_' ...
    num2str(wave.f * 1e-9) 'GHz.fig']);
