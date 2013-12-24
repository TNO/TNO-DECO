function [d,p] = deco_corr(x,y,maxlag)

%x=x(:)
%y=y(:y)
M = length(x);

if (maxlag>M) maxlag = M; end;
    
X = fft(x,2^nextpow2(2*M-1));
Y = fft(y,2^nextpow2(2*M-1));
c = ifft(X.*conj(Y));
p = -maxlag:maxlag;
d = [c(end - maxlag+1:end); c(1:maxlag+1)];