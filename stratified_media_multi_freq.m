function varargout = stratified_media_multi_freq(k0, krho, z, varargin)
%STRATIFIED_MEDIA Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668;
    if strcmp(varargin{1}, 'Superstrate')
        air_length = varargin{2};
        slab_length = varargin{3};
        dielectric_er = varargin{4};

        k0 = permute(repmat(k0, [1 1 size(krho, 1, 2)]), [3 4 2 1]);
        z = repmat(z, 1, 1, size(k0, 3));
        slab_length = permute(repmat(slab_length, ...
            [1 1 size(krho, 1, 2)]), [3 4 2 1]);

        % Air propagation vector
        air_kz = - 1j * sqrt( - k0 .^ 2 + krho .^ 2 );

        % Dielectric propagation vector
        dielectric_k = k0 * sqrt(dielectric_er);
        dielectric_kz = - 1j * sqrt( - dielectric_k .^ 2 + krho .^ 2 );

        % Air impedance
        Zair_te = wave_impedance * k0 ./ air_kz;
        Zair_tm = wave_impedance * air_kz ./ k0;

        % Substrate impedance
        Zs = wave_impedance / sqrt(dielectric_er);
        Zs_te = Zs * dielectric_k ./ dielectric_kz;
        Zs_tm = Zs * dielectric_kz ./ dielectric_k;

        % Input impedance seen at interface 1
        Zin_te = Zs_te ...
            .* (Zair_te + 1j * Zs_te .* tan(dielectric_k .* slab_length)) ...
            ./ (Zs_te + 1j * Zair_te .* tan(dielectric_k .* slab_length));
        Zin_tm = Zs_tm ...
            .* (Zair_tm + 1j * Zs_tm .* tan(dielectric_k .* slab_length)) ...
            ./ (Zs_tm + 1j * Zair_tm .* tan(dielectric_k .* slab_length));

        % Reflection coefficients at interface 1
        gamma1_te = (Zin_te - Zair_te) ./ (Zin_te + Zair_te);
        gamma1_tm = (Zin_tm - Zair_tm) ./ (Zin_tm + Zair_tm);

        % Reflection coefficients at interface 2
        gamma2_te = (Zair_te - Zs_te) ./ (Zair_te + Zs_te);
        gamma2_tm = (Zair_tm - Zs_tm) ./ (Zair_tm + Zs_tm);

        % Air medium 2 voltage and current
        V3p_te = (1 + gamma1_te) .* (1 + gamma2_te) ...
            .* exp(1j * air_kz .* (2 * air_length + slab_length)) ...
            .* exp(-1j * dielectric_kz .* slab_length) ...
            ./ ( (gamma1_te + exp(2j * air_kz * air_length)) ...
            .* (1 + gamma2_te .* exp(-2j * dielectric_kz .* slab_length)) );
        v_te = V3p_te .* exp(-1j * air_kz .* z);
        i_te = v_te ./ Zair_te;
        V3p_tm = (1 + gamma1_tm) .* (1 + gamma2_tm) ...
            .* exp(1j * air_kz .* (2 * air_length + slab_length)) ...
            .* exp(-1j * dielectric_kz .* slab_length) ...
            ./ ( (gamma1_tm + exp(2j * air_kz * air_length)) ...
            .* (1 + gamma2_tm .* exp(-2j * dielectric_kz .* slab_length)) );
        v_tm = V3p_tm .* exp(-1j * air_kz .* z);
        i_tm = v_tm ./ Zair_tm;

        varargout{1} = v_te;
        varargout{2} = i_te;
        varargout{3} = v_tm;
        varargout{4} = i_tm;
    elseif strcmp(varargin{1}, 'SemiInfiniteSuperstrate')
        air_length = varargin{2};
        dielectric_er = varargin{3};

        k0 = permute(repmat(k0, [1 1 size(krho, 1, 2)]), [3 4 2 1]);
        z = repmat(z, 1, 1, size(k0, 3));
        air_length = repmat(air_length, [size(krho, 1, 2, 3)]);

        % Air propagation vector
        air_kz = - 1j * sqrt( - k0 .^ 2 + krho .^ 2 );

        % Dielectric propagation vector
        dielectric_k = k0 * sqrt(dielectric_er);
        dielectric_kz = - 1j * sqrt( - dielectric_k .^ 2 + krho .^ 2 );

        % Air impedance
        Zair_te = wave_impedance * k0 ./ air_kz;
        Zair_tm = wave_impedance * air_kz ./ k0;

        % Dielectric impedance
        Zs = wave_impedance / sqrt(dielectric_er);
        Zs_te = Zs * dielectric_k ./ dielectric_kz;
        Zs_tm = Zs * dielectric_kz ./ dielectric_k;

        % Reflection coefficients at interface 1
        gamma_te = (Zs_te - Zair_te) ./ (Zs_te + Zair_te);
        gamma_tm = (Zs_tm - Zair_tm) ./ (Zs_tm + Zair_tm);

        % Dielectric medium voltage and current
        v_te = exp(-1j * air_kz .* air_length) ...
            .* exp(1j * dielectric_kz .* air_length) ...
            .* (1 + gamma_te) .* exp(-1j * dielectric_kz .* z) ...
            ./ (1 + gamma_te .* exp(-2j * air_kz .* air_length));
        i_te = v_te ./ Zs_te;
        v_tm = exp(-1j * air_kz .* air_length) ...
            .* exp(1j * dielectric_kz .* air_length) ...
            .* (1 + gamma_tm) .* exp(-1j * dielectric_kz .* z) ...
            ./ (1 + gamma_tm .* exp(-2j * air_kz .* air_length));
        i_tm = v_tm ./ Zs_tm;

        varargout{1} = v_te;
        varargout{2} = i_te;
        varargout{3} = v_tm;
        varargout{4} = i_tm;
    else
        error('Error. Invalid argument.');
    end
end

