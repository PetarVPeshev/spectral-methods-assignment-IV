function [k_vector, k, kx, ky, kz] ...
    = wave_vector_multi_freq(relat_permit, k0, sph_grid)
%WAVE_VECTOR This function calculates the wave vector and medium
%propagation constant
%   Detailed explanation goes here
% FIXME: taking k0 as input and then calculating k is confusing, function
% should take f or k (in case of k, do not return k)
    theta = repmat(sph_grid(:, :, 1), [1 1 length(k0)]);
    phi = repmat(sph_grid(:, :, 2), [1 1 length(k0)]);

    k = repmat(k0, [1 1 size(sph_grid, 1, 2)]) * sqrt(relat_permit);
    k = permute(k, [3 4 2 1]);

    kx = k .* sin(theta) .* cos(phi);
    ky = k .* sin(theta) .* sin(phi);
    kz = - 1j * sqrt(- k .^ 2 + kx .^ 2 + ky .^ 2);
    
    k_vector = NaN( [size(k, 1, 2, 3), 3] );
    k_vector(:, :, :, 1) = kx;
    k_vector(:, :, :, 2) = ky;
    k_vector(:, :, :, 3) = kz;
    k_vector = permute(k_vector, [1 2 4 3]);
end