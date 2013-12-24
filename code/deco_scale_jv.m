function [Y,mX,stdX] = scale_jv(X,method,cat)


[m,n ] =size(X); % get a size of the intial matrix

if (nargin<3) % get default category vector
    cat = ones(m,1);
end

if (nargin<2) % set default scaling to meanscale
    method='mean';
end

method = upper(method); % convert method to upper case


idx = find(cat>0); % search all samples whithin model

mX = mean(X(idx,:));   % average of features in model
stdX = std(X(idx,:));  % stdev of features in model
rX = range(X(idx,:));  % range of features in model

stdX(find(stdX==0))=1; % set standard deviation to 1 for constant features
rX(find(rX==0))=1;     % set range to 1 for constant features

switch method % apply the scaling method to the data
    case 'AUTO'  % autoscaleing
        Y = (X-ones(m,1)*mX) ./(ones(m,1)*stdX);
    case 'MEAN' % meancentering
        Y = (X-ones(m,1)*mX);
    case 'RANGE' % range scaling
        Y = (X-ones(m,1)*mX)./ (ones(m,1) * rX);
    otherwise
        Y = X;
        
        
end




        
    