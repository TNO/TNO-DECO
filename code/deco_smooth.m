function sm = deco_smooth(data,window)
%######################################################################
% Function sm = smooth(data,window);
% Purpose: quick and dirty method to smooth spectra using window function
% input paramters:
% X  input data vectro
% M  smoothing window

%output Parameters
% sm  smoothed vector
%#####################################################################
%# written by J.T.W.E. Vogels 15/7/2004
%#####################################################################

if (nargin<2)
    window = 25; 
end

data = real(data); % get real part of the spectrum

[n,m] = size(data);

if (n>m)
    data= data';
    [n,m] = size(data);
end

sm = zeros(1,m);
for i=1:m   % for all points
     lmin = max(1,i-window); % lower limit
     lmax = min(m,i+window); % upper limit
     diff = lmax-lmin+1;
        sm(i) = sum(data(1,lmin:lmax))/diff; % smooth 
    end
end
        