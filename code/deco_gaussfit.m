function [errVal,symmval,A,yh,sigma]=deco_gaussfit(x,y)
%----------------------------------------   
% deco gaussfit(xaxis,peakprofile)
% returns fit-error,peak-symmetry,area,height
% written by J.T.W.E. Vogels
% adapted by F. v/d Kloet 
% TNO Zeist 2008
%-----------------------------------------
 y = y';
 [m,i] = max(y);
 [Ar,m,s] = deco_moments(y);
 [bestcoeffs]=fminsearch(@gaussf,[double(s) i double(sum(y))],optimset('display','off'),x,y);
 sigma = bestcoeffs(1);
 mu = bestcoeffs(2);
 A = bestcoeffs(3);

 yh = A/(sigma*sqrt(2*pi))*exp(-(x-mu) .^ 2 / (2*sigma ^2));
 
 errVal = double(sqrt(sum((y-yh).^2)));
 yhm = yh;
 yhm(yhm>0.01*m)=1;

 symmval = double(deco_peaksymmetry(y.*yhm));
 

%% gauss function
function out=gaussf(coeff,x,Y)
 sigma = coeff(1);
 mu = coeff(2);
 A = coeff(3);
 y = A/(sigma*sqrt(2*pi))*exp(-(x-mu) .^ 2 / (2*sigma ^2));
 DIFF = y - Y; 
 SQ_DIFF = DIFF.^2; 
 out = sum(SQ_DIFF);

%% Deco moments 
function [a,m,s] = deco_moments(g)
x=1:length(g);
m0 = sum(g);  % area 
m1 = sum( x .* g) / m0; % average rt 
m2 = sqrt(sum( (x - m1) .^ 2 .* g) / m0); % calc std deviation
m3 = sum( (x - m1) .^ 3 .* g) / m0; % calc asymmetry of zone profile
%sprintf('Area=%f Rt=%f  m2:%f m3:%f m4:%f\n',m0,m1,m2,m3,m4)
a = m0; % area
m = m1; % pos
s = m2; % std dev 