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
wave.f = linspace(10, 20, 1001) * 1e9;
stratification.h = 10 * 1e-3;
stratification.er = 25;
N = 1001;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;

krho_norm = linspace(1, sqrt(stratification.er), N);

krho_te = NaN(1, length(wave.f));
krho_tm = NaN(1, length(wave.f));
for f_idx = 1 : 1 : length(wave.f)
    krho = krho_norm * wave.k0(f_idx);

    [krho_te(f_idx), krho_tm(f_idx)] = find_krho(wave.k0(f_idx), krho, ...
        'SemiInfiniteSuperstrate', stratification.h, stratification.er);
end
    
%% PLOT PROPAGATION CONSTANT FOR CONSTANT PERMITTIVITY
kd = wave.k0 * sqrt(stratification.er);
figure('Position', [250 250 750 400]);
plot(wave.f * 1e-9, real(krho_te ./ kd), 'LineWidth', 2.0, ...
    'Color', [0 0.4470 0.7410], 'DisplayName', '\Re\{TE1\}');
hold on;
plot(wave.f * 1e-9, imag(krho_te ./ kd), '--', 'LineWidth', 2.0, ...
    'Color', [0 0.4470 0.7410], 'DisplayName', '\Im\{TE1\}');
hold on;
plot(wave.f * 1e-9, real(krho_tm ./ kd), 'LineWidth', 2.0, ...
    'Color', [0.8500 0.3250 0.0980], 'DisplayName', '\Re\{TM1\}');
hold on;
plot(wave.f * 1e-9, imag(krho_tm ./ kd), '--', 'LineWidth', 2.0, ...
    'Color', [0.8500 0.3250 0.0980], 'DisplayName', '\Im\{TM1\}');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('f / GHz');
ylabel('k_{\rho} / k_{d}');
title(['Normalized k_{\rho} @ Semi-Infinite Superstrate, h = ' ...
    num2str(stratification.h * 1e3) ' mm, ' ...
    'and \epsilon_{r} = ' num2str(stratification.er(end))]);
saveas(gcf, ['figures\semi_superstrate_te1_tm1_er_' ...
    num2str(stratification.er) '.fig']);
