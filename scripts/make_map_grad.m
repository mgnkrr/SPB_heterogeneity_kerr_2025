function H = make_map_grad(~)
% MAKE_MAP_GRAD  Build a gradient GHF map using μ, A, θ⊥ and X/Y loaded from a MAT file.
%   H = make_map_grad(cfg)
%
% Required:
%    mu                 : target mean (mW m^-2)
%    A                  : half-range (mW m^-2)  [sign flips direction]
%    theta_perp_deg     : angle PERPENDICULAR to gradient (deg CCW from +x)
%    mat_path           : path to MAT with grid coordinates (e.g., coldex_icethk.mat)
%
% Optional:
%    x_name /  y_name: variable names in MAT (default: 'Xgrid','Ygrid')
%    Gmin,  Gmax     : physics bounds (defaults 20..150)
%    exact_mean         : true => enforce mean after clamping (default true)
%    title              : custom figure title
%
% Returns: struct with fig, ax, cb, G, X, Y, phi

% --------- defaults ---------
%must = @(f) (isfield(cfg,f) && ~isempty( (f)));
% assert(must('mu') && must('A') && must('theta_perp_deg') && must('mat_path'), ...
%     'Provide  mu,  A,  theta_perp_deg,  mat_path');

mu = 42.87;
A  = 104.88;               % negative is ok; direction auto-flips
theta = -25.29;   % angle ⟂ to gradient
mat_path = 'datasets_for_gmin/coldex_icethk.mat';  % path to your ice-thickness MAT
% Optional: if your file uses different names for grids
x_name = 'Xgrid';  y_name = 'Ygrid';
exact_mean = true;
title('GHF gradient from ice-thickness grid');

%mu    = double( mu);
%A     = double( A);
%theta = double( theta_perp_deg);
Gmin  = double(20);
Gmax  = double(150);

% --------- load X/Y from MAT ---------
S = load(mat_path);
[X, Y] = extract_XY_from_mat(S,  x_name,  y_name);

% sanity on units: expect meters; if looks like km, convert to m
if max(abs(X(:))) < 1e4 && max(abs(Y(:))) < 1e4
    % looks like kilometers; convert to meters for math
    X = X * 1000; Y = Y * 1000;
end

% --------- gradient direction (φ) from perpendicular θ⊥ ---------
phi = theta - 90;                 % convert to gradient azimuth
if A < 0                          % absorb sign into direction so A is half-range
    phi = phi + 180;
    A   = abs(A);
end

% --------- build unit gradient U ∈ [-1, +1] along φ ---------
t  = deg2rad(phi); nx = cos(t); ny = sin(t);
xc = X - mean(X(:),'omitnan');
yc = Y - mean(Y(:),'omitnan');
P  = xc*nx + yc*ny;
pmin = min(P(:)); pmax = max(P(:));
if ~isfinite(pmin) || ~isfinite(pmax) || pmax==pmin
    U = zeros(size(P));
else
    U = 2*(P - (pmin+pmax)/2) / (pmax - pmin);
end

% --------- raw field and physics gating ---------
G_raw = mu + A * U;

G = gate_with_exact_mean(G_raw, mu, Gmin, Gmax, 1e-6);

% --------- persistent, interactive figure ---------
x_km = X(1,:)/1000;  y_km = Y(:,1)/1000;
fig = figure('Color','w','Name','GHF map (interactive)','NumberTitle','off','Visible','on');
ax  = axes(fig);
imagesc(ax, x_km, y_km, G);
axis(ax,'image'); set(ax,'YDir','normal'); box(ax,'on');
xlabel(ax,'x (km)'); ylabel(ax,'y (km)');

%settitle(sprintf('GHF (\\mu=%.2f, A=%.2f, \\phi=%.1f°)', mu, A, phi));

colormap(ax, parula);
cb = colorbar(ax); cb.Label.String = 'mW m^{-2}';
clim(ax, [Gmin Gmax]);

fprintf('Grid %dx%d | G[min,max]=[%.2f, %.2f] | φ=%.2f°\n', size(G,1), size(G,2), min(G(:)), max(G(:)), phi);

% return
H = struct('fig',fig,'ax',ax,'cb',cb,'G',G,'X',X,'Y',Y,'phi',phi);
end

% ================== helpers ==================

function [X, Y] = extract_XY_from_mat(S, xname, yname)
% Try exact names, then common fallbacks, then meshgrid if vectors.
candidatesX = unique([string(xname), "Xgrid","X","x","X_km","Xkm"]);
candidatesY = unique([string(yname), "Ygrid","Y","y","Y_km","Ykm"]);

X = []; Y = [];

for nm = candidatesX
    if isfield(S, nm); X = S.(nm); break; end
end
for nm = candidatesY
    if isfield(S, nm); Y = S.(nm); break; end
end

% Some files store inside a struct (e.g., D.Xgrid)
if isempty(X) || isempty(Y)
    fns = fieldnames(S);
    for k = 1:numel(fns)
        v = S.(fns{k});
        if isstruct(v)
            for nm = candidatesX
                if isfield(v, nm); X = v.(nm); end
            end
            for nm = candidatesY
                if isfield(v, nm); Y = v.(nm); end
            end
        end
    end
end

assert(~isempty(X) && ~isempty(Y), 'Could not find X/Y in MAT. Provide  x_name/ y_name.');

% If vectors, build grid
if isvector(X) && isvector(Y)
    [X, Y] = meshgrid(X(:).', Y(:));
end

% Ensure X and Y match sizes
assert(isequal(size(X), size(Y)), 'X and Y must have same size.');
end

function Gc = gate_with_exact_mean(G0, mu_target, L, U, tol)
% Clamp to [L,U] and enforce mean(Gc) ≈ mu_target after clamping via bisection.
if nargin<5 || isempty(tol), tol = 1e-6; end
assert(isfinite(L) && isfinite(U) && L<=U, 'Invalid bounds');

if mu_target <= L
    Gc = min(max(G0, L), U); Gc(:) = L; return
elseif mu_target >= U
    Gc = min(max(G0, L), U); Gc(:) = U; return
end

% Scalar objective: mean(clamp(G0+b)) - mu_target
f = @(b) mean(min(max(G0 + b, L), U), 'all', 'omitnan') - mu_target;

gfin   = G0(isfinite(G0));
spread = max(U - L, max(abs(gfin - mean(gfin,'omitnan'))) * 4 + 1);
blo = -spread; bhi =  spread;
flo = f(blo);   fhi = f(bhi);

k = 0;
while (~(flo <= 0 && fhi >= 0)) && k < 20
    blo = blo - spread; bhi = bhi + spread;
    flo = f(blo);       fhi = f(bhi);
    k = k + 1;
end

if ~(flo <= 0 && fhi >= 0)
    % Rare fallback
    H  = min(max(G0, L), U);
    mH = mean(H(:),'omitnan');
    Gc = min(max(H - (mH - mu_target), L), U);
    return
end

for k = 1:60
    bmid = 0.5*(blo + bhi);
    fmid = f(bmid);
    if abs(fmid) <= tol, break; end
    if fmid > 0, bhi = bmid; else, blo = bmid; end
end

bsol = 0.5*(blo + bhi);
Gc = min(max(G0 + bsol, L), U);
end
