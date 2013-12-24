function z = deco_correlation(block,nr)
%------------------------------------------------------
% calculate regression (R2) between spectrum nr in block
% and all spectra in next block
%------------------------------------------------------
global project;


if block>= project.last || size(project.deco{block+1},1) <=0
    z =[];
    return;
end

% for all peaks in next block
n = size(project.deco{block+1}.sopt,1); % number of peaks in next block
z = zeros(n,1);
h = project.deco{block}.sopt(nr,:);
pw = project.pw;


for i=1:n    
    warning off Matlab:divideByZero
    zr = corrcoef(h,project.deco{block+1}.sopt(i,:)) .^ 2;   
    warning on Matlab:divideByZero
   
    z(i) = zr(1,2);
    if (z(i)>0.95 && project.deco{block}.inc(nr)==1)% highly correlated samples
        % positie in current block
        p1= zeros(project.nfiles,1);
        p2 = zeros(project.nfiles,1);
        %s1 = zeros(project.nfiles,1);
        %s2 = zeros(project.nfiles,1);
        
        % position in current block 
        for file=1:project.nfiles 
            [hgt,pos] = max(project.deco{block}.copt((file-1)*4*pw+1:file*4*pw,nr));
            p1(file)=pos;
        end
        
        % positie in next block
        for file=1:project.nfiles 
            [hgt,pos] = max(project.deco{block+1}.copt((file-1)*4*pw+1:file*4*pw,i));
            p2(file) = pos;
        end    
        
        % determine the median value of the peak positions
        std1 = std(p1);
        std2 = std(p2);
   
        % select which series of peaks to go with
        % A) series with lowest std deviation in peak position 
        % B) else first series
        
        if std1>std2 % select series 2
            project.deco{block+1}.inc(i) = 2;
            project.deco{block}.inc(nr)= 0;
        else
            project.deco{block+1}.inc(i) = 0;
            project.deco{block}.inc(nr) = 3;
        end
    end
end
