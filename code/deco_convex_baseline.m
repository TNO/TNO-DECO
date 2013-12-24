function [v,c,l,m] = deco_convex_baseline(data,window)
%----------------------------------------%
% use convexhull to calculate baseline
% usage [v,c] = convex_baseline(data,16)
% input : data (data vector)
%       :  window = smoothing window
% output : v new baseline
%        : c corrected spectrum
%------------------------------------------
%J. Vogels 6-11-2007 (c) tno KvL Zeist
% added final point to list 27-05-2008
% be carefull not to set the window too high
%----------------------------------------%
data = data(:);
if nargin < 2
    window = 16;
end

specdiff = max(data) - min(data);
if (specdiff<1.0E-14) 
    %disp('no signal')
    v = ones(length(data),1);
    c = zeros(length(data),1);
    return;
end
     
p1 = deco_smooth(double(data),window);
x=1:length(data);

k=convhull(x,p1,{'Qt','Pp'}); % convex hull on smoothed data
lenx = max(size(x));
[a,b]=min(k);
l=circshift(k,-(b-1));
idx=diff(l)>0;
l=l(idx);
l = unique([l; lenx]); % add the final point to the list 
% l has at least two points
v=interp1(x(l),data(l),x); % compute baseline
v=v(:);
idx = find(isnan(v));
v(idx)=0; % set nan points to zero
c = data - v;

