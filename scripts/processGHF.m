% =========================
% processGHF.m
% =========================
function [HF_interp, Gmin, Gdiff, Gmin_adv, Gdiff_adv, Delta_q_adv] = ...
    processGHF(HF_input, Ts, Mb, H, Delta_q_adv, Delta_topo)
% PROCESSGHF  Compute Gmin and model differences on the H grid.
% Robust to transpose mismatches. Optionally applies:
%   - Delta_q_adv (mW/m^2) to compute advected fields.
%   - Colgan unitless topographic correction Delta_topo = dG/G:
%       HF_corr = HF*(1 + Delta_topo)
%
% Usage:
%   [HF_interp,Gmin,Gdiff] = processGHF(HF, Ts, Mb, H)
%   [HF_interp,Gmin,Gdiff,Gmin_adv,Gdiff_adv] = processGHF(HF, Ts, Mb, H, Delta_q_adv)
%   [HF_interp,Gmin,Gdiff,~,~] = processGHF(HF, Ts, Mb, H, [], Delta_topo)

if nargin < 5
    Delta_q_adv = [];
end
if nargin < 6
    Delta_topo = [];
end

% ---- enforce inputs match H grid (transpose-safe) ----
[Ts, ~] = align_to_H_local(Ts, H, 'Ts');
[Mb, ~] = align_to_H_local(Mb, H, 'Mb');
[HF_interp, ~] = align_to_H_local(HF_input, H, 'HF_input');

if ~isempty(Delta_q_adv)
    [Delta_q_adv, ~] = align_to_H_local(Delta_q_adv, H, 'Delta_q_adv');
end

if ~isempty(Delta_topo)
    [Delta_topo, ~] = align_to_H_local(Delta_topo, H, 'Delta_topo');
end

% ---- Colgan topographic correction (unitless) ----
if ~isempty(Delta_topo)
    HF_interp = HF_interp .* (1 + Delta_topo);
end

T0 = 273.15;
Cp = 2009;        % heat capacity of ice
rho = 910;        % ice density
n = 3;            % Glen's exponent
gamma = 8.7e-4;
secperyear = 31556926;
K = 2.1 * secperyear;  % W/yr/m/K
dzeta = 0.05;
zeta = (1:-dzeta:0)';  % 0 = surface

lambda = (zeta.^(n+3) - 1) / (n + 1) / (n + 3) - ...
         (n + 2) * (zeta.^2 - 1) / 2 / (n + 1) + zeta - 1;

W = zeros(size(H,1), size(H,2), numel(zeta), 'like', H);

% Accumulate W
for i = 2:numel(zeta)
    expo = 0.5 * (lambda(i) + lambda(i-1)) .* H .* Mb .* (rho * Cp / K);
    W(:,:,i) = W(:,:,i-1) - exp(expo) .* dzeta;
end

denom = H .* (W(:,:,1) - W(:,:,end));
denom(abs(denom) < 1e-8) = NaN;

Gmin = K * (T0 - gamma .* H - Ts) ./ denom;
Gmin = Gmin / secperyear * 1000;  % -> mW/m^2

% Differences
Gdiff = HF_interp - Gmin;

% Optional advection-adjusted fields
Gmin_adv  = [];
Gdiff_adv = [];
if ~isempty(Delta_q_adv)
    Gmin_adv  = Gmin - Delta_q_adv;
    Gdiff_adv = HF_interp - Gmin_adv;
end

end

% ---- local helper ----
function [Aout, didTranspose] = align_to_H_local(Ain, H, name)
didTranspose = false;

if isequal(size(Ain), size(H))
    Aout = Ain;
    return;
end

if isequal(size(Ain), fliplr(size(H)))
    Aout = Ain.';
    didTranspose = true;
    return;
end

error('processGHF:SizeMismatch', ...
    '%s is %s but H is %s. Provide data on H grid or transpose/regrid.', ...
    name, mat2str(size(Ain)), mat2str(size(H)));
end
