close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');
c = physconst('LightSpeed');
wave_impedance = 376.730313668;

stratification.er = NaN(1, 1002);

%% PARAMETERS
wave.f = linspace(9, 11, 1001) * 1e9;
stratification.h = 15 * 1e-3;
stratification.er(1 : 1001) = linspace(1, 25, 1001);
stratification.er(1002) = 12;
N = 1001;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
% frequency - col; er - row
stratification.hs = wave.wavelength ./ (4 * sqrt(stratification.er'));
stratification.hs(1002, :) = 2.1 * 1e-3;
dipole.L = wave.wavelength / 2;
dipole.W = wave.wavelength / 20;

krho_te = NaN(length(stratification.er), length(wave.f));
krho_tm = NaN(length(stratification.er), length(wave.f));
for er_idx = 1 : 1 : length(stratification.er)
    krho_norm = linspace(1, sqrt(stratification.er(er_idx)), N);

    for f_idx = 1 : 1 : length(wave.f)
        krho = krho_norm * wave.k0(f_idx);
    
        [krho_te(er_idx, f_idx), krho_tm(er_idx, f_idx)] ...
            = find_krho(wave.k0(f_idx), krho, 'Superstrate', ...
            stratification.h, stratification.hs(er_idx, f_idx), stratification.er(er_idx));
    end
end
    
%% PLOT PROPAGATION CONSTANT FOR CONSTANT PERMITTIVITY
figure('Position', [250 250 750 400]);
plot(stratification.h ./ wave.wavelength, real(krho_te(1002, :) ./ wave.k0), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.h ./ wave.wavelength, imag(krho_te(1002, :) ./ wave.k0), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.h ./ wave.wavelength, real(krho_tm(1002, :) ./ wave.k0), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.h ./ wave.wavelength, imag(krho_tm(1002, :) ./ wave.k0), ...
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
    num2str(stratification.hs(1002, 1) * 1e3) ' mm, ' ...
    'and \epsilon_{r} = ' num2str(stratification.er(1002))]);
saveas(gcf, ['figures\superstrate_krho_er_const_' ...
    num2str(stratification.er(1002)) '.fig']);
    
%% PLOT PROPAGATION CONSTANT
figure('Position', [250 250 750 400]);
plot(stratification.er(1 : 1001), real(krho_te(1 : 1001, 501) ./ wave.k0(501)), ...
    'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{TE1\}');
hold on;
plot(stratification.er(1 : 1001), imag(krho_te(1 : 1001, 501) ./ wave.k0(501)), ...
    '--', 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{TE1\}');
hold on;
plot(stratification.er(1 : 1001), real(krho_tm(1 : 1001, 501) ./ wave.k0(501)), ...
    'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{TM1\}');
hold on;
plot(stratification.er(1 : 1001), imag(krho_tm(1 : 1001, 501) ./ wave.k0(501)), ...
    '--', 'LineWidth', 2.0, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{TM1\}');
grid on;
xlim([min(stratification.er(1 : 1001)) max(stratification.er(1 : 1001))]);
legend show;
xticks(min(stratification.er(1 : 1001)) : 2 : max(stratification.er(1 : 1001)));
legend('location', 'bestoutside');
xlabel('\epsilon_{r}');
ylabel('k_{\rho} / k_{0}');
title(['Normalized k_{\rho} @ Superstrate, f = ' ...
    num2str(wave.f(501) * 1e-9) ' GHz, h = ' ...
    num2str(stratification.h * 1e3) ' mm, ' ...
    'h_{s} = \lambda_{0} / (4 * sqrt(\epsilon_{r}))']);
saveas(gcf, ['figures\superstrate_krho_f_const_' ...
    num2str(wave.f(501) * 1e-9) 'GHz.fig']);
