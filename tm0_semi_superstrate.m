close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');

c = physconst('LightSpeed');
wave_impedance = 376.730313668;

Nf = 101;

%% PARAMETERS
wave.f = linspace(13, 17, Nf) * 1e9;
stratification.h = 10 * 1e-3;
stratification.er = 10;
N = 1001;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;

%% SPHERICAL COORDINATE SYSTEM
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, 100);
phi = linspace(0, 2 * pi, 400);
sph_grid = meshgrid_comb(theta, phi);

krho_tm0 = NaN(1, length(wave.f));
for f_idx = 1 : 1 : length(wave.f)
    krho_tm0(f_idx) = find_krho_tm0(wave.k0(f_idx), ...
        'SemiInfiniteSuperstrate', stratification.h, stratification.er);
end
    
%% PLOT PROPAGATION CONSTANT FOR CONSTANT PERMITTIVITY
krho_tm0_norm = krho_tm0 ./ (wave.k0 * sqrt(stratification.er));
figure('Position', [250 250 750 400]);
plot(wave.f * 1e-9, real(krho_tm0_norm), 'LineWidth', 2.0, ...
    'Color', [0 0.4470 0.7410], 'DisplayName', '\Re\{TM0\}');
hold on;
plot(wave.f * 1e-9, imag(krho_tm0_norm), '--', 'LineWidth', 2.0, ...
    'Color', [0 0.4470 0.7410], 'DisplayName', '\Im\{TM0\}');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('f / GHz');
ylabel('k_{\rho} / k_{d}');
title(['Normalized TM0 k_{\rho} @ Semi-Infinite Superstrate, h = ' ...
    num2str(stratification.h * 1e3) ' mm, ' ...
    'and \epsilon_{r} = ' num2str(stratification.er)]);
saveas(gcf, ['figures\semi_superstrate_tm0_er_' ...
    num2str(stratification.er) '.fig']);
