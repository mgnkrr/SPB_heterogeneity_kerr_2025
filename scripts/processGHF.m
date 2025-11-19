function [HF_interp, Gmin, Gdif] = processGHF(HF_input, Ts, Mb, H)

    T0 = 273.15;
    Cp = 2009; % heat capacity of ice
    rho = 910; % ice density
    n = 3; % Glen's exponent
    gamma = 8.7e-4;
    secperyear = 31556926;
    K = 2.1 * secperyear;  % W/yr/m/K
    dzeta = 0.05;
    zeta = (1:-dzeta:0)'; % scaled coordinate system, 0 = surface of the ice sheet

    lambda = (zeta.^(n+3) - 1) / (n + 1) / (n + 3) - ...
             (n + 2) * (zeta.^2 - 1) / 2 / (n + 1) + zeta - 1;

    W = zeros(size(H,1), size(H,2), numel(zeta));
    
    for i = 2:numel(zeta)
        W(:,:,i) = W(:,:,i-1) - ...
            exp(0.5 * (lambda(i) + lambda(i-1)) .* H .* Mb * ...
            rho * Cp / K) * dzeta;
    end
    
    denom = H .* (W(:,:,1) - W(:,:,end));
    denom(abs(denom) < 1e-8) = NaN;

    Gmin = K * (T0 - gamma .* H - Ts) ./ denom;
    Gmin = Gmin / secperyear * 1000;

    HF_interp = HF_input;  % in case it’s already interpolated
    Gdif = HF_interp - Gmin;
end

