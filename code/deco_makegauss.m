function [respons] = deco_makegauss(x,mu,sigma);
%------------------------------------------------
% this function calculates a gaussian respons
% based upon the following input:
% x     : x axis scale
% mu    : mean
% sigma : standard deviation
% F. v/d Kloet juni 2008
%-----------------------------------------------
% function [respons] = makegauss(x,mu,sigma);
% example:
% x = 1:100; mu = 50; sigma = 10;
% y = makegauss(x,mu,sigma)
% plot(x,y);
% xlabel('x');
% ylabel('y');
% title('gaussian peak with mean = 50 and sigma = 10');

if nargin==0,
     help makegauss;
    return
end

respons = 1/(sigma*sqrt(2*pi)) * exp(-0.5*((x-mu)/sigma).^2);
