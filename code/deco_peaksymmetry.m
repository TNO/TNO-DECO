function [y] = deco_peaksymmetry(data,level)

%---------------------------------------------------
%calculate the symmmetry of a peak 
%at level ( default at 10% of max peak intensity)
% written by J.Vogels 27/03/2008
%--------------------------------------------------

if nargin < 2
    level =3;
end

if level <= 0 
    level =1;
end

if level >= 100
    level =99;
end

[mn,pos] = max(data);
idx = find(data < (mn/level));
pval = abs(min(idx(find(idx>pos))) - pos); %#ok<FNDSB>
mval = abs(max(idx(find(idx<pos))) - pos); %#ok<FNDSB>


% altered by frans
if isempty(pval) || isempty(mval)
    if(isempty(pval))
        [pval, mp] = min(data(pos:end));
        if (pos-mp+1>0) 
            mval = data(pos-mp+1);
        else
            mval=0.0; %#ok<NASGU>
            y=0.0;
            return;
        end
    else
        [mval, mp] = min(data(1:pos));
        if pos+mp-1<length(data)
            pval = data(pos+mp-1);
        else
            pval=0.0; %#ok<NASGU>
            y=0.0;
            return;
        end
    end
end    
    
y = min(pval/mval,mval/pval);





