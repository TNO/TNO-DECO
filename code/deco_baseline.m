function [matrix]=deco_baseline(matrix,bl_options,name,method,current,all)
%
% I/O: [matrix]=rem_bl_gui(matrix,options,name);
% The function removes baselines per mass trace.
% Input:
%   matrix:                matrix with time (rows) * mass (columns) Intensities
%   options_bl:            structure
%   eilers
%   bl_options.lambda:     smoothing parameter (see asysm.m)
%   bl_options.p:          asymmetrie parameter (see asysm.m)
%   bl_options.d:          order of differences (see asysm.m)
% groningen
%   bl_options.gsize       slize size
%   bl_options.gmove       slizemove
%   bl_options.gth         threshold
%   name                   is name of the fieles
% Output:
%   matrix:                baseline corrected matrix
if  nargin <= 2 
    name=' ';
end

[dir,fname] = fileparts(name);
str = ['Baseline removal:' fname '( ' num2str(current) ' of ' num2str(all) ') in progress, please wait'];
h=waitbar(0,str);

[traces,masses]=size(matrix);
if ~isfield(bl_options,'TH');
    bl_options.TH=0;
end

for i=1:masses,
    % only use those data points for which signal > 0;
    [idx]=find(matrix(:,i)~=0);
    if length(idx)>(traces/100) && max(matrix(idx,i))>bl_options.TH, % minimum amount of data points at mass #i > 1% of total traces;
        if method==1 % Paul Eiler method
            if bl_options.lambda ~= 0,
                matrix(idx,i)=matrix(idx,i)-single(deco_asysm(double(matrix(idx,i)),bl_options.lambda,bl_options.asym,bl_options.order));
            end
        elseif method==2 % groningen method
            matrix(idx,i) = matrix(idx,i)-deco_medianbaseline(matrix(idx,i),bl_options.gsize,bl_options.gmove,bl_options.gthresh);
        elseif method==3 % Frans method
            matrix(:,i) = matrix(:,i)-deco_blcor(matrix(:,i),bl_options.basewin);
        elseif method==4 % convex hull (not yet tested)
            matrix(:,i) = matrix(:,i) - deco_convex_baseline(matrix(:,i),16);
        else
          matrix(idx,i)=0;
        end
    end
    waitbar(i/masses,h);
    set(h,'Name',[num2str(deco_fixfloat(100*i/masses,2)),'% data processed (mass:', num2str(i),')']);
end
close(h);
