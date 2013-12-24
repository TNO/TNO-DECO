function c = deco_peakpos(x, npeaks, c)

CORCOEFFL = 0.95;
pos = zeros(npeaks, 1);
val = zeros(npeaks, 1);
mc = [];
for i = 1:npeaks
    [r, a] = max(max(x,[], 2));
    [m, b] = max(max(x,[], 1));
    c(a, i) = r;
    mc = [mc, x(:, b)];
    x(:, b) = 0;
end

end