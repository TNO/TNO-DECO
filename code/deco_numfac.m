function [nsig,blk,sim] = deco_numfac(data,maxn)

if nargin < 2
    maxn=10;
end

disp('deco numfac')

factor = 2;
ssq = sum(sum(deco_mncn(data) .^2));
simblock = randn(size(data)) * (ssq/numel(data)) .^ (1/2);
ssq1 = sum(sum(simblock .^2));

blk = svd(data,0);    blk = blk .^ 2;
sim = svd(simblock,0); sim = sim .^2;

nsig = 1;
while(nsig <= maxn && blk(nsig)>sim(nsig)* factor) 
    nsig = nsig+1;
end

blk=blk(1:10)';
sim=sim(1:10)';
nsig
