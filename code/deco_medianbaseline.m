function [baseline] = deco_medianbaseline(LC_MS, slicesize,slicemove,threshold)
% median baseline uisng sliding windows
% P. Horvatovic, 2006.01.16
% 
% parameters
%input
% LC_MS         lc-ms data in matrix form (column retention time, rows m/z)
% slicesize     length of windows to calculate the median (gen. value 1200)
% slicemove     displacement of sliding windows (resolution of baseline) (gen. value 10)
% threshold     smallest value to avoid unprobable baseline (gen. value 200)
%
%outpout
% baseline      median baseline of the same size and format than the LC_MS
%-----------------------------------------

[a,b]=size(LC_MS);
middlepoint=round((slicesize/2)-(slicemove/2));
nslice=round((a-(slicesize-slicemove/2))/slicemove)-1;
baseline=zeros(a,b);
for k=1:b,
    ttemp=median(LC_MS(1:slicesize,k));
    if ttemp<=threshold, ttemp=threshold; end
    baseline(1:middlepoint,k)=ttemp;
    for i=0:nslice,
        ttemp=median(LC_MS((i*slicemove+1):(i*slicemove+slicesize),k));
        if ttemp<=threshold, ttemp=threshold; end
        baseline((i*slicemove+middlepoint+1):(i*slicemove+middlepoint+slicemove),k)=ttemp;
    end
    ttemp=median(LC_MS((nslice*slicemove+1):end,k));
    if ttemp<=threshold, ttemp=threshold; end
    baseline(((nslice+1)*slicemove+middlepoint+1):end,k)=ttemp;
end