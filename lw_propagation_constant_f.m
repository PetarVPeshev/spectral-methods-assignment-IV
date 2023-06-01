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
wave.f = linspace(9, 11, 1001) * 1e9;
stratification.h = 15 * 1e-3;
stratification.hs = 2.1 * 1e-3;
stratification.er = 12;
N = 1001;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
dipole.L = wave.wavelength / 2;
dipole.W = wave.wavelength / 20;

krho_norm = linspace(1, sqrt(stratification.er), N);

krho_te = NaN(1, length(wave.f));
krho_tm = NaN(1, length(wave.f));
for f_idx = 1 : 1 : length(wave.f)
    krho = krho_norm * wave.k0(f_idx);

    [krho_te(f_idx), krho_tm(f_idx)] = find_krho(wave.k0(f_idx), krho, ...
        'Superstrate', stratification.h, stratification.hs, stratification.er);
end
    
%% PLOT PROPAGATION CONSTANT
figure('Position', [250 250 750 400]);
plot(stratification.h ./ wave.wavelength, real(krho_te ./ wave.k0), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.h ./ wave.wavelength, imag(krho_te ./ wave.k0), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.h ./ wave.wavelength, real(krho_tm ./ wave.k0), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.h ./ wave.wavelength, imag(krho_tm ./ wave.k0), ...
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
    num2str(stratification.hs * 1e3) ' mm, and \epsilon_{r} = ' ...
    num2str(stratification.er)]);
saveas(gcf, ['figures\propagation_constant_er_' ...
    num2str(stratification.er) '.fig']);
