%=======================================================
function [peakL0, npeakL0, C] = deco_eilersL0(y,kappa0)
%=======================================================
% Written by S. Krishnan, 09/06/2009
% for TNO deco, peak identification
%
% input(s):
% y: regressand data, e.g. Chromatogram data
% kappa0: kappa1 is the coefficient for L0 penalty
%
% output(s):
% peakL0: Estimated unknowns with Eiler's L0 penalty,
% e.g. identified peaks
%=======================================================
% y0 = sum(y,2);
% y = y0/max(y0); % Normalisation
kappa1 = 0.01;  % kappa1 is the coefficient for L1 penalty
[peakL1, peakL0, C]=deco_f_fit_L1L0(y,kappa1,kappa0);

% Number of peaks
dpp = diff(peakL0);
dppo = padarray(dpp, 1, 'post');
dppr = padarray(dpp, 1, 'pre');
dps = dppo - dppr;
dps = padarray(dps, 1, 1);
pid0 = find(dps < 0);
cond0 = dps(pid0) ./ dps(pid0+1);
cond1 = dps(pid0) ./ dps(pid0-1);
pid1 = find(cond0 < -1 & cond1 < -1);
npeakL0 = size(pid1,1);

end

%=======================================================
function [a1,a0,C]=deco_f_fit_L1L0(y,kappa1,kappa0)
%=======================================================
% Written originally by Paul. H. Eilers,
% Deconvolutie van pulsen
% met ZEN (L_0) penalty
% Modified by S. Krishnan, 09/06/2009
% added non-negativity constrain + some additional
% modification to use the function for peak 
% identification for TNO Deco.
%=======================================================
m=length(y);

% Simuleer pulsvorm
nc = 41;
u = (1:nc)' / nc;
sigma = 0.1;
uu = 2 * (u - 0.5) / sigma;
c = exp(-uu .^ 2 / 2);
c = c / max(c);

% Build convolution matrix
C = zeros(m + nc - 1, m);
for k = 1:m;
    C((1:nc) + k - 1, k) = c;
end
C = C(floor(nc / 2) + (1:m), :);
n = m;

% Pulses on random positions
% q = zeros(n, 1);
% p = cumsum(floor(exp(rand(m, 1))* 3 + 4));
% p = p(p < n);
% q(p) = rand(length(p), 1);
% y = C * q + randn(m, 1) * 0.03;

% Fit met L1 penalty
w = ones(n, 1);
%kappa1 = 0.01;
s0 = 0;
beta = 0.0001;
for it = 1:50
    W = spdiags(w, 0, n, n);
    a1 = (W + C' * C) \ (C' * y);
    w = kappa1 ./ sqrt(beta + a1 .^ 2);
    z = C * a1;
    r = y - z;
    s = r' * r + kappa1 * a1' * a1;
    disp(s - s0);
    s0 = s;
end

% Fit met L0 penalty
w = ones(n, 1);
%kappa0 = 0.003;
s0 = 0;
beta = 0.0003;
for it = 1:50
    W = spdiags(w, 0, n, n);
    a0 = (W + C' * C) \ (C' * y);
    idx0 = find(a0<0); 
    a0(idx0)=0; % Non-negativity constrain
    w = kappa0 ./ (beta + a0 .^ 2);
    z = C * a0;
    r = y - z;
    s = r' * r + kappa0 * a0' * a0;
    disp(s - s0);
    s0 = s;
end

end