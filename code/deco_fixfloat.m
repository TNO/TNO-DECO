function y=deco_fixfloat(x,n)
% Rounds off the number x to n numbers after decimal point.
%
% Input
% ------
% x    the number.
% n    the numbers after decimal point.
%
% Output
% ------
% y    the result
%
% I/O: y=fixfloat(x,n);
%
% See also built-in functions ROUND, FLOOR, CEIL, FIX.
%
% Written by Florian Wulfert.

y=round(x*10^n)/10^n;