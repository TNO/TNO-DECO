function [bl] = deco_blcor(aData, aWindowSize)

%-------------------------------------------------------------------------
% function bl = blcor(aData, aWindowSize)
% Baseline estimation function
%--------------------------------------
% Input variables
% aData       : signal vector
% aWindowSize : minimal peak width, 15 is good estimate
%--------------------------------------
% Output variables
% bl          : baseline estimate
%---------------------------------------
% version 0.96
% TNO,Frans van der Kloet, 25 jan 2006
%-------------------------------------------------------------------------

aData = aData(:);
lMin = min(aData);                          % calculate offset (=min)
[bl1]= localBlCor(aData,aWindowSize);       % estimate first baseline points
lidx = find(aData>bl1');                    % find those points that are larger than the baseline
l=length(aData);                            % how many peak points ?
xas = 1:l;                  
r = setdiff(xas,lidx);                      % which points are baseline ?

if(r(1)~=1)
    r = [1,r];                              % add first point (if necessary)
end
if(r(end)~=l)
    r = [r,l];                              % add last point (if necessary)
end

bl1 = interp1(xas(r),aData(r),xas);         % interpolate
bl = localMaFilter(bl1',6*aWindowSize);     % smooth

return;


function [bl,rem] = localBlCor(aData,aWindowSize)

l=length(aData);
ldata = aData(:);
windowSize = aWindowSize;

mstd = localMovingstd(ldata,windowSize);

imstd = find(mstd~=0);
moduslData = localMode(mstd(imstd)); % use the modus for baseline level determination

xas = 1:l;
idx = find(mstd>10*moduslData); % if signal>10*modus the peakdata
rem = setdiff(xas,idx); % else baseline data

if (~isempty(idx) && idx(1)==1)
    rem = [1,rem]; % first point is always baseline
    idx(1) = [];
end

if (~isempty(idx) && idx(end)==l)
    rem = [rem,l]; % last point is always baseline
    idx(end) = [];            
end


blvals = ldata(rem); % copy baseline values
mabl = localMaFilter(blvals,min(floor(length(blvals)/10),6*windowSize)); % smooth baseline values use 6 * windowSize
bl = interp1(xas(rem),mabl,xas); % linear interpolation

return;


%------------------------------------------------------
function [res] = localMaFilter(aData,aWindowSize)
%------------------------------------------------------
% moving average filter
% Frans van der Kloet, 2006

if aWindowSize<=1 
    res=aData;
    return
end

f=ones(aWindowSize,1)/aWindowSize; 

n=size(aData,1);
res=filter2(f,aData);
m2=floor(aWindowSize/2);
n2=ceil(aWindowSize/2)-1;

try  
    res=res([zeros(1,m2)+m2+1,(m2+1):(n-n2),zeros(1,n2)+(n-n2)],:);
catch
end

return % return to main function (blcor)

function s = localMovingstd(x,k,windowmode)
% movingstd: efficient windowed standard deviation of a time series
% usage: s = movingstd(x,k,windowmode)
%
% Movingstd uses filter to compute the standard deviation, using
% the trick of std = sqrt((sum(x.^2) - n*xbar.^2)/(n-1)).
% Beware that this formula can suffer from numerical problems for
% data which is large in magnitude.
%
% arguments: (input)
%  x   - vector containing time series data
%
%  k   - size of the moving window to use (see windowmode)
%        All windowmodes adjust the window width near the ends of
%        the series as necessary.
%
%        k must be an integer, at least 1 for a 'central' window,
%        and at least 2 for 'forward' or 'backward'
%
%  windowmode - (OPTIONAL) flag, denotes the type of moving window used
%        DEFAULT: 'central'
%
%        windowmode = 'central' --> use a sliding window centered on
%            each point in the series. The window will have total width
%            of 2*k+1 points, thus k points on each side.
%        
%        windowmode = 'backward' --> use a sliding window that uses the
%            current point and looks back over a total of k points.
%        
%        windowmode = 'forward' --> use a sliding window that uses the
%            current point and looks forward over a total of k points.
%
%        Any simple contraction of the above options is valid, even
%        as short as a single character 'c', 'b', or 'f'. Case is
%        ignored.
%
% arguments: (output)
%  s   - vector containing the windowed standard deviation.
%        length(s) == length(x)

% check for a windowmode
if (nargin<3) || isempty(windowmode)
  % supply the default: 
  windowmode = 'central';
elseif ~ischar(windowmode)
  error 'If supplied, windowmode must be a character flag.'
end
% check for a valid shortening.
valid = {'central' 'forward' 'backward'};
windowmode = lower(windowmode);
ind = strmatch(windowmode,valid);
if isempty(ind)
  error 'Windowmode must be a character flag: ''c'', ''b'', or ''f''.'
else
  windowmode = valid{ind};
end

% length of the time series
n = length(x);

% check for valid k
if (nargin<2) || isempty(k) || (rem(k,1)~=0)
  error 'k was not provided or not an integer.'
end
switch windowmode
  case 'central'
    if k<1
      error 'k must be at least 1 for windowmode = ''central''.'
    end
    if n<(2*k+1)
      error 'k is too large for this short of a series and this windowmode.'
    end
  otherwise
    if k<2
      error 'k must be at least 2 for windowmode = ''forward'' or ''backward''.'
    end
    if (n<k)
      error 'k is too large for this short of a series.'
    end
end

% we will need the squared elements 
x2 = x.^2;

% split into the three windowmode cases for simplicity
A = 1;
switch windowmode
  case 'central'
    B = ones(1,2*k+1);
    s = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/(2*k+1)))/(2*k));
    s(k:(n-k)) = s((2*k):end);
  case 'forward'
    B = ones(1,k);
    s = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/k))/(k-1));
    s(1:(n-k+1)) = s(k:end);
  case 'backward'
    B = ones(1,k);
    s = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/k))/(k-1));
end

% special case the ends as appropriate
switch windowmode
  case 'central'
    % repairs are needed at both ends
    for i = 1:k
      s(i) = std(x(1:(k+i)));
      s(n-k+i) = std(x((n-2*k+i):n));
    end
  case 'forward'
    % the last k elements must be repaired
    for i = (k-1):-1:1
      s(n-i+1) = std(x((n-i+1):n));
    end
  case 'backward'
    % the first k elements must be repaired
    for i = 1:(k-1)
      s(i) = std(x(1:i));
    end
end

return; % return to caller (localBlCor)


function [M,F,C] = localMode(x,dim)
%MODE   Mode, or most frequent value in a sample.
%   M=MODE(X) for vector X computes M as the sample mode, or most frequently
%   occurring value in X.  For a matrix X, M is a row vector containing
%   the mode of each column.  For N-D arrays, MODE(X) is the mode of the
%   elements along the first non-singleton dimension of X.
%
%   When there are multiple values occurring equally frequently, MODE
%   returns the smallest of those values.  For complex inputs, this is taken
%   to be the first value in a sorted list of values.
%
%   [M,F]=MODE(X) also returns an array F, of the same size as M.
%   Each element of F is the number of occurrences of the corresponding
%   element of M.
%
%   [M,F,C]=MODE(X) also returns a cell array C, of the same size
%   as M.  Each element of C is a sorted vector of all the values having
%   the same frequency as the corresponding element of M.
%
%   [...]=MODE(X,DIM) takes the mode along the dimension DIM of X.
%
%   This function is most useful with discrete or coarsely rounded data.
%   The mode for a continuous probability distribution is defined as
%   the peak of its density function.  Applying the MODE function to a
%   sample from that distribution is unlikely to provide a good estimate
%   of the peak; it would be better to compute a histogram or density
%   estimate and calculate the peak of that estimate.  Also, the MODE
%   function is not suitable for finding peaks in distributions having
%   multiple modes.
%
%   Example: If X = [3 3 1 4
%                    0 0 1 1
%                    0 1 2 4]
%
%   then mode(X) is [0 0 1 4] and mode(X,2) is [3
%                                               0
%                                               0]
%
%   To find the mode of a continuous variable grouped into bins:
%      y = randn(1000,1);
%      edges = -6:.25:6;
%      [n,bin] = histc(y,edges);
%      m = mode(bin);
%      edges([m, m+1])
%      hist(y,edges+.125)
%
%   Class support for input X:
%      float:  double, single
%
%   See also MEAN, MEDIAN, HIST, HISTC.

%   Copyright 2005 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2005/06/21 19:24:03 $

error(nargchk(1,2,nargin, 'struct'))    
if nargin<2 
    % Determine which dimension to use
    dim = find(size(x)~=1, 1);
    if isempty(dim)
      dim = 1;
    end
else
    if ~isscalar(dim) || ~isa(dim,'double') || dim~=floor(dim) ...
                      || dim<1              || ~isreal(dim)
        error('MATLAB:mode:BadDim',...
              'DIM argument must be a scalar specifying a dimension of X.');
    end
end

if ~isfloat(x)
    error('MATLAB:mode:InvalidInput',...
          'X must be an array of double or single numeric values.')
end

dofreq = nargout>=2;
docell = nargout>=3;
wassparse = issparse(x);

sizex = size(x);
if dim>length(sizex)
    sizex = [sizex, ones(1,dim-length(sizex))];
end

sizem = sizex;
sizem(dim) = 1;

% Set up outputs with the proper dimension and type
if wassparse
    M = sparse(sizem(1),sizem(2));  % guaranteed to be 2-D double
else
    M = zeros(sizem,class(x));
end
if dofreq
    F = zeros(sizem);
end
if docell
    C = cell(sizem);
end

% Dispose of empty arrays right away
if isempty(x)
    if docell
        C(:) = {M(1:0)};  % fill C with empties of the proper type
    end
    if prod(sizem)>0
        M(:) = NaN;
        if dofreq
            F(:) = 0;
        end
    end
    return
end

% Convert data to operate along columns of a 2-d array
x = permute(x,[dim, (1:dim-1), (dim+1:length(sizex))]);
x = reshape(x,[sizex(dim),prod(sizem)]);
[nrows,ncols] = size(x);

% Loop over these columns
for j=1:ncols
    v = sort(x(:,j));                        % sorted data
    start = find([1; v(1:end-1)~=v(2:end)]); % start of run of equal values
    freq = [start(2:end);nrows+1] - start;   % frequency of these values
    [maxfreq,firstloc] = max(freq);          % find most frequent
    M(j) = v(start(firstloc));               % smallest most frequent
    if dofreq
        F(j) = maxfreq;                      % highest frequency
    end
    if docell
        if wassparse
            C{j} = sparse(v(start(freq==maxfreq)));
        else
            C{j} = v(start(freq==maxfreq));  % all most frequent
        end
    end
end


return; % return to caller function (localBlCor)