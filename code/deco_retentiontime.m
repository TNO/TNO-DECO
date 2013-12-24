function y = deco_retentiontime(block)
%---------------------------------------------------------------
% calculate average (median) retention time of peak in a block
% written by j.t.w.e. vogels tno-zeist (1-06-2008)
%---------------------------------------------------------------
global project;
f=project.interval;
nfiles = project.nfiles;
pw= project.pw;

if size(project.deco{block},1)==0 % zero components
    y=0;
    return;
end
    
N = size (project.deco{block}.copt,2); % number of compounds in block
rt = zeros(N,nfiles);
p  = zeros(N,nfiles);

for file=1:nfiles
    [hgt,pos] = max(project.deco{block}.copt((file-1)*4*pw+1:file*4*pw,:)); % hetzelfde als quality parameters
    rt(:,file) = (f*(pos+(block-1)*2*pw)+project.rt_start(file))/60;
    p(:,file) = pos;    
end

%if block==3
%     p
%end

for i=1:N
    % verwijder pieken op meer dan 
    idx1 = find(p(i,:) >= 5);
    idx2 = find(p(i,:) <= 4*pw-5);
    iid = intersect(idx1,idx2);
    
    if ~isempty(iid)
        pp = median(p(i,iid)); % determine median value      
        rsd=std(p(i,iid));             % get standard deviation on p-position
    else
        pp = [];
        rsd = [];
    end
    
    adx1 = find(p(i,:) >= pw);
    adx2 = find(p(i,:) <= 3*pw);
    aad = intersect(adx1,adx2);
    %if block==3
    %      aad
    %end
    flag =0;
    if isempty(aad)
        flag = 2; % zeker geen piek in hotzone
    end
      
    
    if ~isempty(iid)
        g= median(rt(i,iid));
       
        if (block==3) 
            pp;
        end
        if pp<pw || pp>3*pw && flag==0
            
            flag = 1; % possible exclusion of peak
        end
        project.deco{block}.rt(i) = g; % get mediaan retentiontime
        project.deco{block}.sd(i) = rsd;% store standard deviation
        project.deco{block}.pp(i) = pp;% median position
        project.deco{block}.flag(i) = flag; % 
    else
        project.deco{block}.rt(i) = median(rt(i,:));
        project.deco{block}.sd(i) = std(rt(i,:));
        project.deco{block}.pp(i) = median(p(i,:));
        project.deco{block}.flag(i) = flag;
    end  
end

y=1;


