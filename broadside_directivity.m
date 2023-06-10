function dir_br = broadside_directivity(relat_permit, field, sph_grid, r)
%DIRECTIVITY This function calculates the directivity, radiation intensity,
% and radiated power
%   Detailed explanation goes here
    wave_impedance = 376.730313668 / sqrt(relat_permit);

    dth = sph_grid(1, 2, 1) - sph_grid(1, 1, 1);
    dph = sph_grid(2, 1, 2) - sph_grid(1, 1, 2);
    theta = repmat(sph_grid(:, :, 1), [1 1 size(field, 4)]);

    field_total = sqrt( abs(field(:, :, 1, :)) .^ 2 ...
        + abs(field(:, :, 2, :)) .^ 2 + abs(field(:, :, 3, :)) .^ 2);
    field_total = permute(field_total, [1 2 4 3]);
    
    rad_intensity = (field_total .^ 2) * (r ^ 2) / (2 * wave_impedance);
    [row, col] = find( isnan(rad_intensity) );
    if ~isempty(col)
        warning(['Found singularities (NaN) in radiation intensity; ' ...
            'set NaN cells to zero in radiation intensity calculation.']);
        rad_intensity(row, col) = 0;
    end

    rad_power = sum(rad_intensity .* sin(theta), [1 2]) * dth * dph;
    rad_power = permute(repmat(permute(rad_power, [3 1 2]), ...
        [1 1 size(theta, 1, 2)]), [3 4 1 2]);
    
    dir = 4 * pi * rad_intensity ./ rad_power;
    dir_br = permute(dir(1, 1, :), [2 3 1]);
end
