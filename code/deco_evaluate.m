function y = deco_evaluate()

global project;

% determine retention times
for i=project.first:project.last
    deco_retentiontime(i);
end

%evaluate overlap
for i=project.first:project.last-1
    if size(project.deco{i},1) > 0 % is bock defined ? 
        n= size(project.deco{i}.sopt,1); % number of peaks per block
        s = ones(n,1);
        project.block{i}.inc=s;
        for j=1:n % examine each peak
            deco_correlation(i,j);
        end
    end
end
y=1;
