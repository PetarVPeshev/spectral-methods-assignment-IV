function field = farfield_multi_freq(k, r, sph_grid, kz, z, sgf, ...
    current_ft, varargin)
%FARFIELD This function evaluates the electric or magnetic far-field in 
%spherical coordinates for a given Fourier Transform of electric or 
%magnetic current density
%   Detailed explanation goes here
    z_source = 0;
    if ~isempty(varargin)
        z_source = varargin{1};
    end

    k = permute(repmat(k, [1 1 size(kz, 1, 2)]), [3 4 2 1]);
    z = repmat(z, [1 1 size(kz, 3)]);

    theta = repmat(sph_grid(:, :, 1), [1 1 size(kz, 3)]);
    phi = repmat(sph_grid(:, :, 2), [1 1 size(kz, 3)]);

    % z is the observation point, z_source is the location of 
    % the source on the z-axis
    field_const = 1j * kz .* exp(1j * kz .* abs(z - z_source)) ...
        .* exp(- 1j * k * r) / (2 * pi * r);

    % both sgf and current_ft must be 2D
    sgf_xx = permute(sgf(:, :, 1, 1, :), [1 2 5 3 4]);
    sgf_xy = permute(sgf(:, :, 1, 2, :), [1 2 5 3 4]);
    sgf_yx = permute(sgf(:, :, 2, 1, :), [1 2 5 3 4]);
    sgf_yy = permute(sgf(:, :, 2, 2, :), [1 2 5 3 4]);
    sgf_zx = permute(sgf(:, :, 3, 1, :), [1 2 5 3 4]);
    sgf_zy = permute(sgf(:, :, 3, 2, :), [1 2 5 3 4]);
    current_ft_x = permute(current_ft(:, :, 1, :), [1 2 4 3]);
    current_ft_y = permute(current_ft(:, :, 2, :), [1 2 4 3]);

    field_x = field_const ...
        .* (sgf_xx .* current_ft_x + sgf_xy .* current_ft_y);
    field_y = field_const ...
        .* (sgf_yx .* current_ft_x + sgf_yy .* current_ft_y);
    field_z = field_const ...
        .* (sgf_zx .* current_ft_x + sgf_zy .* current_ft_y);

    field = NaN( [size(kz, 1, 2, 3), 3] );
    field(:, :, :, 1) = field_x .* sin(theta) .* cos(phi) ...
        + field_y .* sin(theta) .* sin(phi) + field_z .* cos(theta);
    field(:, :, :, 2) = field_x .* cos(theta) .* cos(phi) ...
        + field_y .* cos(theta) .* sin(phi) - field_z .* sin(theta);
    field(:, :, :, 3) = - field_x .* sin(phi) + field_y .* cos(phi);

    field = permute(field, [1 2 4 3]);
end

