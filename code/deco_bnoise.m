function bas_noise = deco_bnoise(spc,NN)
%---------------------------------------------------------------------
% BNOISE - Calculate spectral noise.
%          Required in   mean_bspts.
%          As described in
%              Sergey Golotvin, Antony Williams, Journal of Magnetic
%              Resonance, Vol. 146, No. 1, Sep 2000, pp. 122-125
%--------------------------------------------------------------------
npts = length(spc);
if ~exist('NN')
  NN   = 16;           % sections in the spectrum to calc noise
end
% For spectra with less than 64 points
while NN>=npts
  NN=ceil(NN/2);
end

for k=1:NN
    Knoise(k) = std(spc(round ((npts/NN)*(k-1) + (1:npts/NN))));
end

bas_noise = min(Knoise);
return
