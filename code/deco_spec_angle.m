function y = deco_spec_angle(r1,r2)
%------------------------------------------------
% similarity between vectors = 
% cos (angle) between them 
% (a.b) = |a||b|.cos(angle)
% cos(angle) = (a.b) / (|a||b|)
% written by J.T.W.E. Vogels to replace calcsimm
% 28-08-2008
%------------------------------------------------
ss1 = sum(r1 .* r1); 
ss2 = sum(r2 .* r2);
ss3 = sum(r1 .* r2);
y = ss3 / (sqrt(ss1)*sqrt(ss2));