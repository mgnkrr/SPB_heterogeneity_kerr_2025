%% spb_rms_from_specularity_optionB.m
% Reproduce Schroeder et al. (2015) Section III-style forward constraint:
%   Specularity content Sc as a function of ice/water interface RMS roughness,
%   using:
%     - YOUR survey geometry (h, d, eps_r, fc)
%     - Schroeder apertures (Option B): L1=700 m, L2=2000 m
%     - Approximate aperture angle spans Phi1, Phi2 from L1, L2 and your geometry
%   Then "read off" the RMS roughness required for Sc >= thresholds.
%
% Key outputs:
%   1) Plot Sc vs RMS (your Figure-5-like curve)
%   2) RMS limits for Sc thresholds (e.g., 0.10, 0.20, 0.30)
%   3) Helpful sanity prints: Phi1, Phi2, phiF, and ratios Phi/phiF
%
% IMPORTANT:
% - This is an order-of-magnitude physical constraint, not an inversion.
% - Phi_i approximation is the biggest assumption; if you later extract
%   the true angle spans from your focuser, swap them in directly.

clear; close all;

%% ---------------- USER INPUTS ----------------
% Radar / medium
fc_Hz  = 60e6;     % MARFA center frequency [Hz]
eps_r  = 3.15;     % relative permittivity of ice
c0     = 299792458;

% Your representative geometry (set from SPB domain)
h_m    = 500;      % survey height above surface [m]
d_m    = 3000;     % ice thickness [m]  (use median/typical in your region)

% Option B apertures: match Schroeder et al. (2015)
L1_m   = 700;      % [m]
L2_m   = 2000;     % [m]

% RMS sweep (cm)
rms_cm = linspace(0.1, 50, 700); % 0.1–50 cm
rms_m  = rms_cm/100;

% Specularity thresholds of interest
Sc_targets = [0.10, 0.20, 0.30];

% Angle mapping choice for Phi_i:
%   'optical' : R_eff = h + d*sqrt(eps_r)  (conservative, includes air+ice path scale)
%   'iceonly' : R_eff = d*sqrt(eps_r)      (bed-path scale)
phi_range_model = 'iceonly';

%% ---------------- DERIVED GEOMETRY ----------------
sqrt_eps = sqrt(eps_r);

% (1) Fresnel-zone angular width phiF (Schroeder Eq. 13–14)
D1   = sqrt( (2*c0/fc_Hz) * (h_m + d_m*sqrt_eps) ); % first Fresnel diameter-like scale
phiF = 2*atan( D1/(2*d_m) );                        % radians

% (2) Aperture angle spans Phi1, Phi2 from Schroeder L1/L2 and your geometry (Option B)
switch lower(phi_range_model)
    case 'optical'
        R_eff = h_m + d_m*sqrt_eps;
    case 'iceonly'
        R_eff = d_m*sqrt_eps;
    otherwise
        error('phi_range_model must be ''optical'' or ''iceonly''.');
end

Phi1 = 2*atan( (L1_m/2) / R_eff ); % radians
Phi2 = 2*atan( (L2_m/2) / R_eff ); % radians

if Phi2 <= Phi1
    error('Phi2 <= Phi1. Check h_m, d_m, eps_r or phi_range_model.');
end

%% ---------------- SCHROEDER-STYLE SCATTERING MODEL ----------------
% g(RMS) (Schroeder Eq. 10)
g = (4*pi .* rms_m .* fc_Hz .* sqrt_eps) ./ c0;

% erf-model focused energies (Schroeder Eq. 12, using their erf approximation)
% Ei = erf( (Phi_i/phiF)/(2g) ) * erf( 1/(2g) )
arg_common = 1 ./ (2*g);

E1 = erf( (Phi1./phiF) .* arg_common ) .* erf( arg_common );
E2 = erf( (Phi2./phiF) .* arg_common ) .* erf( arg_common );

% Convert E1, E2 to specularity content Sc using the S + D*(Phi/pi) mixing model:
%   E1 = S + D*(Phi1/pi)
%   E2 = S + D*(Phi2/pi)
% Solve for S, D then Sc = S/(S+D)
D = (E2 - E1) ./ ((Phi2 - Phi1)/pi);
S = E1 - D .* (Phi1/pi);
Sc = S ./ (S + D);

% Clamp numerical noise
Sc = max(0, min(1, Sc));

%% ---------------- REPORT / SANITY PRINTS ----------------
fprintf('\n--- Geometry / params ---\n');
fprintf('h = %.0f m, d = %.0f m, fc = %.1f MHz, eps_r = %.2f\n', h_m, d_m, fc_Hz/1e6, eps_r);
fprintf('Apertures (Option B): L1 = %.0f m, L2 = %.0f m\n', L1_m, L2_m);
fprintf('Phi mapping model: %s (R_eff = %.0f m)\n', phi_range_model, R_eff);
fprintf('Phi1 = %.3f deg, Phi2 = %.3f deg\n', rad2deg(Phi1), rad2deg(Phi2));
fprintf('phiF = %.3f deg\n', rad2deg(phiF));
fprintf('Phi1/phiF = %.3f, Phi2/phiF = %.3f\n', Phi1/phiF, Phi2/phiF);

fprintf('\n--- RMS roughness required for given S_c thresholds ---\n');
for k = 1:numel(Sc_targets)
    target = Sc_targets(k);
    idx = find(Sc >= target, 1, 'last'); % Sc decreases with RMS
    if isempty(idx)
        fprintf('S_c >= %.2f: not achieved within RMS sweep (%.1f–%.1f cm).\n', target, rms_cm(1), rms_cm(end));
    else
        fprintf('S_c >= %.2f: RMS <= %.2f cm (approx; within sweep grid)\n', target, rms_cm(idx));
    end
end

%% ---------------- PLOT ----------------
figure('Color','w');
plot(rms_cm, Sc, 'LineWidth', 2);
grid on;
xlabel('Ice/water interface RMS roughness (cm)');
ylabel('Specularity content, S_c');
title(sprintf('S_c vs RMS (Option B: L1=%.0fm, L2=%.0fm; h=%.0fm, d=%.0fm; %s)', ...
    L1_m, L2_m, h_m, d_m, phi_range_model));

hold on;
for k = 1:numel(Sc_targets)
    yline(Sc_targets(k), '--', sprintf('S_c = %.2f', Sc_targets(k)), ...
        'LabelVerticalAlignment','bottom');
end
hold off;

%% ---------------- OPTIONAL: quick sensitivity on d (uncomment) ----------------
%{
d_list = [1500, 2000, 2500, 3000]; % m
targets = [0.10 0.20];
fprintf('\n--- Sensitivity: RMS limit vs d (keeping h fixed) ---\n');
for dd = d_list
    d_m_tmp = dd;

    % recompute phiF
    D1_tmp   = sqrt( (2*c0/fc_Hz) * (h_m + d_m_tmp*sqrt_eps) );
    phiF_tmp = 2*atan( D1_tmp/(2*d_m_tmp) );

    % recompute Phi using same apertures, updated geometry
    switch lower(phi_range_model)
        case 'optical'
            R_eff_tmp = h_m + d_m_tmp*sqrt_eps;
        case 'iceonly'
            R_eff_tmp = d_m_tmp*sqrt_eps;
    end
    Phi1_tmp = 2*atan( (L1_m/2) / R_eff_tmp );
    Phi2_tmp = 2*atan( (L2_m/2) / R_eff_tmp );

    % energies
    E1_tmp = erf( (Phi1_tmp./phiF_tmp) .* arg_common ) .* erf( arg_common );
    E2_tmp = erf( (Phi2_tmp./phiF_tmp) .* arg_common ) .* erf( arg_common );

    D_tmp = (E2_tmp - E1_tmp) ./ ((Phi2_tmp - Phi1_tmp)/pi);
    S_tmp = E1_tmp - D_tmp .* (Phi1_tmp/pi);
    Sc_tmp = max(0, min(1, S_tmp ./ (S_tmp + D_tmp)));

    fprintf('d = %.0f m: ', d_m_tmp);
    for t = targets
        idx = find(Sc_tmp >= t, 1, 'last');
        if isempty(idx)
            fprintf('Sc>=%.2f: NA  ', t);
        else
            fprintf('Sc>=%.2f: RMS<=%.2f cm  ', t, rms_cm(idx));
        end
    end
    fprintf('\n');
end
%}
