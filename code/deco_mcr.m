function [c,s] = deco_mcr(x,c0,options)
%disp('don mymcr')

%MCR Multivariate curve resolution with constraints
%  INPUTS:
%        x = the matrix to be decomposed as X = CS, and
%       c0 = the initial guess for (c) or (s) depending on its size.
%            For X (M by N) then C is (M by K) and S is (K by N) where K is the
%            number of factors that is determined from the size of input (c0).
%            If (c0) is size (M by K) it is the initial guess for (c) (also if M==N).
%            If (c0) is size (K by N) it is the initial guess for (s).
%  OPTIONAL INPUT:
%   options = structure with the following fields:
%     display: [ 'off' | {'on'} ]      governs level of display to command window,
%       plots: [ 'none' | {'final'} ]  governs level of plotting,
%        ccon: [ 'none' | {'nonneg'} ] non-negativity on concentrations,
%        scon: [ 'none' | {'nonneg'} ] non-negativity on spectra,
%          cc: [ ]  concentration equality constraints,
%              MxK matrix with NaN where equality contraints are not applied
%              and real value of the constraint where they are applied,
%          sc: [ ]  spectra equality constraints,
%              KxN matrix with NaN where equality contraints are not applied
%              and real value of the constraint where they are applied,
%        sclc: [ ]  concentration scale axis,
%              vector with M elements otherwise 1:M is used,
%        scls: [ ]  spectra scale axis,
%              vector with N elements otherwise 1:N is used,
%        tolc: [ {1e-5} ]  tolerance on non-negativity for concentrations,
%        tols: [ {1e-5} ]  tolerance on non-negativity for spectra,
%       ittol: [ {100} ]   convergence criteria,
%              if 0<ittol<1 then this is convergence tolerance {default = 1e-8}, and
%              if ittol>=1 and integer, this is maximum number of
%              iterations.
%  OUTPUTS:
%         c = are the estimated concentrations, and 
%         s = estimated pure component spectra.
%  Note: unconstrained factors have spectra scaled to unit length.
%  This can result in output spectra with different scales. 
%I/O: [c,s] = mcr(x,c0,options);
%I/O: mcr demo

options = [];
options.name    = 'options';
options.display = 'off';
options.plots   = 'none';
options.ccon    = 'nonneg';
options.scon    = 'nonneg';
options.cc      = [];
options.sc      = [];
options.sclc    = [];
options.scls    = [];
options.tolc    = 1e-5;  %tolerance on non-negativity for concentrations
options.tols    = 1e-5;  %tolerance on non-negativity for spectra
options.ittol   = 100;

if nargin == 0; x = 'io'; end
varargin{1} = x;

if ischar(varargin{1});
  disp('Here we go checking the options');
  options = [];
  options.name    = 'options';
  options.display = 'off';
  options.plots   = 'none';
  options.ccon    = 'nonneg';
  options.scon    = 'nonneg';
  options.cc      = [];
  options.sc      = [];
  options.sclc    = [];
  options.scls    = [];
  options.tolc    = 1e-5;  %tolerance on non-negativity for concentrations
  options.tols    = 1e-5;  %tolerance on non-negativity for spectra
  options.ittol   = 100; 
  if nargout==0;
      %disp('no out')
      evriio(mfilename,varargin{1},options); 
  else
      %disp('nargout')
      c = options; %evriio(mfilename,varargin{1},options); 
  end
  return; 
end

if nargin<3     %set default options
  disp('kleiner 3');
  options  = deco_mcr('options');
else
  %disp('gelijk 3')  
  %options
  %disp('na envirro')
  %options = reconopts(options,deco_mcr('options'));
  %options
end

warning off backtrace
switch lower(options.ccon)
case {1, 'nonneg'}
  options.ccon = 1;
case {2, 'fastnnls'}
  options.ccon = 2;
case {0, 'none'}
  options.ccon = 0;
otherwise
  warning('Option.ccon not recognized. Reset to ''none''.')
end


switch lower(options.scon)
case {1, 'nonneg'}
  options.scon = 1;
case {2, 'fastnnls'}
  options.scon = 2;
case {0, 'none'}
  options.scon = 0;
otherwise
  warning('Option.scon not recognized. Reset to ''none''.');
end


[m,n]   = size(x);
if size(c0,1)==m
  ka    = size(c0,2);  %initial guess for concentration
  s0    = zeros(ka,n);
  c0int = true;
elseif size(c0,2)==n
  ka    = size(c0,1);  %initial guess for spectra
  s0    = c0;
  c0    = zeros(m,ka);
  c0int = false;
else
  error('c0 must be size(x,1) by #factors or #factors by size(x,2)')
end
if isempty(options.cc)
  cc1 = false;
else
  if (size(options.cc,1)~=m || size(options.cc,2)~=ka)
    error('options.cc must be size(x,1) by #factors')
  elseif any(any(isfinite(options.cc)))
    cc1 = true; disp('Equality Constraints on C')
  else
    cc1 = false;
  end
end
if isempty(options.sc)
  sc1 = false;
else
  if (size(options.sc,1)~=ka || size(options.sc,2)~=n)
    error('options.sc must be #factors by size(x,2)')
  elseif any(any(isfinite(options.sc)))
    sc1 = true; disp('Equality Constraints on S')
  else
    sc1 = false;
  end
end
if isempty(options.sclc)
  options.sclc = 1:m;
elseif length(options.sclc)~=m
  options.sclc = 1:m;
end
if isempty(options.scls)
  options.scls = 1:n;
elseif length(options.scls)~=n
  options.scls = 1:n;
end
if isempty(options.ittol)
  options.ittol = 100;
  itmin         = 1e-8;
  itmax         = options.ittol;
elseif options.ittol<1
  itmin         = options.ittol;
  itmax         = 1e6;
  options.ittol = itmax;
elseif options.ittol<0
  error('options.ittol must be positive.')
else
  itmax         = options.ittol;
  itmin         = 1e-8;
end
if isempty(options.tolc) || options.tolc<0
  options.tolc  = 1e-5;  %tolerance on non-negativity for concentrations
end
if isempty(options.tols) || options.tols<0
  options.tols  = 1e-5;  %tolerance on non-negativity for spectra
end

kac     = ones(1,ka); %keep a one for factors with no constraint
if cc1
  %kac(1,find(any(isfinite(cc)))) = 0;
  %change to look only for factors with equality constraints > options.tolc
  for ii=1:ka
    jc    = find(isfinite(options.cc(:,ii)));
    if any(options.cc(jc,ii)>options.tolc)
      kac(1,ii) = 0;
    end
  end
end
if sc1
  %kac(1,find(any(isfinite(sc')))) = 0;
  %change to look only for factors with equality constraints > options.tols
  for ii=1:ka
    jc    = find(isfinite(options.sc(ii,:)));
    if any(options.sc(ii,jc)'>options.tols)
      kac(1,ii) = 0;
    end
  end
end
kac     = find(kac);    %factors without constraints
if cc1, [jc] = find(isfinite(options.cc)); end % locations of C constraints
if sc1, [js] = find(isfinite(options.sc)); end % locations of S constraints

if c0int==1  % initial guess is for C
  if cc1, c0(jc) = options.cc(jc);   end % C equal constraint
  c     = c0;
  switch options.scon % estimate initial guess for S
  case 0
    s0  = c\x;
  case {1,2}
    s0  = c\x;
    s0(find(s0<-options.tols)) = 0; %options.tols;
  end
  if sc1, s0(js) = options.sc(js);   end % S equal constraints
  s     = s0;
  if ~isempty(kac)
    s(kac,:) = normaliz(s(kac,:)); 
  end
else           %initial guess is for S
  if sc1, s0(js) = options.sc(js);   end % S equal constraints
  s       = s0;
  if ~isempty(kac)
    s(kac,:) = normaliz(s(kac,:));
  end
  switch options.ccon %initial guess for C
  case 0
    c0    = x/s;
  case {1,2}
    c0    = x/s;
    c0(find(c0<-options.tolc)) = 0; %options.tolc;
  end
  if cc1, c0(jc) = options.cc(jc);   end % C equal constraint
  c       = c0;
end

it      = 0;
rv=[];
while it<itmax
  switch options.ccon %solve for concentration
  case 0
    c   = x/s;
  case 1
    c   = x/s;
    c(find(c<-options.tolc)) = 0; %options.tolc;
  case 2
    for jj=1:m, c(jj,:) = fastnnls(s',x(jj,:)',options.tolc,c(jj,:)')'; end
  end
  if cc1, c(jc) = options.cc(jc);   end % C equal constraint

  switch options.scon %solve for spectra
  case 0
    s = c\x;
  case 1
    s = c\x;
    s(find(s<-options.tols)) = 0; %options.tols;
  case 2
    for jj=1:n, s(:,jj) = fastnnls(c,x(:,jj),options.tols,s(:,jj)); end
  end
  if sc1, s(js) = options.sc(js);   end % S equal constraints
  if ~isempty(kac)
    s(kac,:) = normaliz(s(kac,:)); 
  end

  it     = it+1; garb = it;
  if (options.ittol<1)&&((it/2-round(it/2))==0)
%change this to only test residuals for locations where
% equality constraints ain't
    resc = 0; ress = 0;
    for ii=1:ka
      resc = resc+norm(c0(:,ii));
      ress = ress+norm(s0(ii,:));
    end
    ress = sqrt(sum(sum((s'-s0').^2),2)/ka/ress);
    resc = sqrt(sum(sum((c-c0).^2),2)/ka/resc);
    if (ress<itmin)&&(resc<itmin)
      it = itmax+1;
    else
      c0 = c;
      s0 = s;
    end
  end
end

s0     = (x-c*s).^2;
switch lower(options.plots)
case 'final'
  figure
	subplot(2,1,1), plot(options.sclc,sum(s0,2))
	xlabel('Concentration Profile'), ylabel('\bfC\rm Sum Squared Residuals')
	title(sprintf('MCR Results after %d iterations',garb))
	subplot(2,1,2), plot(options.scls,sum(s0))
	xlabel('Spectral Profile'),      ylabel('\bfS\rm Sum Squared Residuals')

	figure
	subplot(2,1,1), plot(options.sclc,c)
	xlabel('Concentration Profile'), ylabel('\bfC')
	title(sprintf('MCR Results after %d iterations',garb))
	subplot(2,1,2), plot(options.scls,s)
	xlabel('Spectral Profile'),      ylabel('\bfS')
end
s0     = sum(sum(s0,2));
switch lower(options.display)
case 'on'
  disp(sprintf('residual    %3.6e', s0))
  disp(sprintf('Unmodelled Variance is %g percent', s0*100/sum(sum(x.^2))))
end

%=======================================
function [ndat,norms] = normaliz(dat)
%=======================================
[m] = size(dat,1);
ndat = dat;
norms = zeros(m,1);
for i = 1:m
  if norm(ndat(i,:)) ~= 0
    norms(i) = norm(ndat(i,:));
    ndat(i,:) = ndat(i,:)/norms(i);
  end
end
  



