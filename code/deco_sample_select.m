function  y=deco_sample_select()
y=1;

global project;
% collect all tics
xdata = zeros(project.nfiles,project.ntraces);
for i=1:project.nfiles
    ytic=deco_readfiletic([project.name '.bin'],i,project.ntraces,project.nmasses);
    plot(ytic);
    xdata(i,:) = ytic;
end

[a,b,c] = svd(xdata,0);


ndim = 3;

Y = pdist(a(:,[1:ndim]),'euclidean');
Z = linkage(Y,'average');
% squareform(Y)
figure
dendrogram(Z);
% 
maxc = 10; % define maximum of 10 clusters 
T = clusterdata(a(:,[1:ndim]),'maxclust',maxc)

% find all sample in cluster 2
sel = zeros(maxc,ndim);
for i = 1:maxc
    f=find (T==i);
    f'
    sel(i,:)=a(f(1),[1:ndim]);
end

figure;
hold on;
plot(a(:,1),a(:,2),'o');
plot(sel(:,1),sel(:,2),'xr');
hold off;
