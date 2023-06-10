function FT = ft_current_multi_freq(k0, varargin)
%FT_CURRENT This function calculates the current's Fourier Transform
%   Detailed explanation goes here
% , k_comp, width, length, antenna, orientation
    antenna = varargin{end - 1};
    orientation = varargin{end};
    
    if strcmp(antenna, 'dipole')
        k_comp = varargin{1};
        width = varargin{2};
        length = varargin{3};
        dielectric_er = varargin{4};

        length ...
            = permute(repmat(length, [1 1 size(k_comp, 1, 2)]), [3 4 2 1]);
        width ...
            = permute(repmat(width, [1 1 size(k_comp, 1, 2)]), [3 4 2 1]);

        kx = permute(k_comp(:, :, 1, :), [1 2 4 3]);
        ky = permute(k_comp(:, :, 2, :), [1 2 4 3]);
        k0 = permute(repmat(k0, [1 1 size(k_comp, 1, 2)]), [3 4 2 1]);
        dielectric_k = k0 * sqrt(dielectric_er);
        keq = (k0 + dielectric_k) / 2;

        T = sinc( ky .* width / (2 * pi) );
        F = 2 * keq .* ( cos(kx .* length / 2) - cos(keq .* length / 2) )...
            ./ ( (keq .^ 2 - kx .^ 2) .* sin(keq .* length / 2) );
        ft = F .* T;

        FT = zeros( [size(k_comp, 1, 2), size(k0, 3), 2] );
        if strcmp(orientation, 'x')
            FT(:, :, :, 1) = ft;
        elseif strcmp(orientation, 'y')
            FT(:, :, :, 2) = ft;
        else
            error('Not implemented.');
        end

        FT = permute(FT, [1 2 4 3]);
    else
        error('Not implemented.');
    end
end

