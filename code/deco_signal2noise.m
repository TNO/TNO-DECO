function y = deco_signal2noise(x)

n=length(x);
mx= max(x);
mi= min(x);

warning off MATLAB:rankDeficientMatrix
% query by entering  [a,b] = lastwarn
dx = sqrt(sum(diff(x) .^ 2)/(n-2));
if (dx ~=0) 
    y = (mx-mi)/dx;
else
    y =0.0;
end


warning on MATLAB:rankDeficientMatrix



