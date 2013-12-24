function [mass,peakh] = deco_unique_mass(x,val)

%-------------------------
% find unique masses
% deco_unique  fixes the problem of multiple peaks per mass
% 3/8/2006 J.V.
%-------------------------
mass = unique(x);
peakh = ones(1,length(mass));
for i=1:length(mass)
    idx=find(x==mass(i));
    peakh(i) = sum(val(idx));
end


