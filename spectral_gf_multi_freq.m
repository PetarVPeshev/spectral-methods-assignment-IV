function sgf = spectral_gf_multi_freq(relat_permit, k, kx, ky, ...
    v_tm, v_te, i_tm, i_te, field, current)
%SPECTRAL_GF Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668 / sqrt(relat_permit);

    k = permute(repmat(k, [1 1 size(kx, 1, 2)]), [3 4 2 1]);
    krho = sqrt(kx .^ 2 + ky .^ 2);

    sgf = NaN( [size(k, 1, 2, 3), 3, 2] );
    if strcmp(field, 'E')
        if strcmp(current, 'J')
            sgf(:, :, :, 1, 1) = - (v_tm .* (kx .^ 2) + v_te .* (ky .^ 2)) ./ (krho .^ 2);
            sgf(:, :, :, 1, 2) = (v_te - v_tm) .* kx .* ky ./ (krho .^ 2);
            sgf(:, :, :, 2, 1) = sgf(:, :, 1, 2);
            sgf(:, :, :, 2, 2) = - (v_te .* (kx .^ 2) + v_tm .* (ky .^ 2)) ./ (krho .^ 2);
            sgf(:, :, :, 3, 1) = wave_impedance * kx .* i_tm ./ k;
            sgf(:, :, :, 3, 2) = wave_impedance * ky .* i_tm ./ k;
        elseif strcmp(current, 'M')
            kz = - 1j * sqrt( - k .^ 2 + krho .^ 2 );
            Ztm = wave_impedance * kz ./ k;
            sgf(:, :, :, 1, 1) = (v_tm - v_te) .* kx .* ky ./ (krho .^ 2);
            sgf(:, :, :, 1, 2) = - (v_tm .* (kx .^ 2) + v_te .* (ky .^ 2)) ./ (krho .^ 2);
            sgf(:, :, :, 2, 1) = (v_tm .* (ky .^ 2) + v_te .* (kx .^ 2)) ./ (krho .^ 2);
            sgf(:, :, :, 2, 2) = (v_te - v_tm) .* kx .* ky ./ (krho .^ 2);
            % TODO: check if correct
            sgf(:, :, :, 3, 1) = - Ztm .* i_tm .* ky ./ kz;
            sgf(:, :, :, 3, 2) = Ztm .* i_tm .* kx ./ kz;
        else
            error('Error. Invalid argument.');
        end
    else
        error('Error. Invalid argument.');
    end

    sgf = permute(sgf, [1 2 4 5 3]);
end

