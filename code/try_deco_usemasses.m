function [npeaks,usemasses,c,s]=try_deco_usemasses(data)
maxpeaks=10;
tr_all=zeros(maxpeaks,1);
for npeaks=1:maxpeaks,
tr_all(npeaks)=deco_peakiterator(data,npeaks);
end
ctr_all=sqrt(1-(tr_all.*tr_all));
[maxtr,cidx]=max(tr_all);
[maxtr,usepeaks,c,s]=deco_peakiterator(data,cidx);
c=c(:,usepeaks);
usemasses=deco_usemasses(data,c);
s=s(:,usepeaks,usemasses);
npeaks=length(usepeaks);
end

function [maxtr,usepeaks,c,s]=deco_peakiterator(data,npeaks)
global project;
pw = project.pw;
options = struct('display','off','plots','none');

sp = zeros(npeaks,pw*4,'single');
delta=(3.5/(npeaks+1));
init_pos=floor(((1:npeaks)*delta+0.25)*pw); % space peaks
for sp_init = 1:npeaks,
    sp(sp_init,:) = deco_makegauss(1:pw*4,init_pos(sp_init),pw*1/6); % create initiation gaussian peaks
end
idx = 1:project.nfiles; % index to the files
actfiles=project.nfiles; % count the number of actual files to include
s = zeros(actfiles,npeaks,length(data),'single');
c = zeros(4*pw,npeaks,'single');
for part=1:length(idx),
    [c(:,:),s(idx(part),:,:)]=deco_mcr(data((idx(part)-1)*4*pw+1:idx(part)*4*pw,:),sp',options);
end
[maxtr,usepeaks]=deco_peak(data,c(:,:),s(idx(part),:,:)); % parameter: data has to corrected for multiple samples;
end

function usemasses=deco_usemasses(x,c)
global project;
s=c\x;
for i=1:size(c,2),
    xcs(i,:,:)=c(:,i)*s(i,:);
    xms(i,:)=max(squeeze(xcs(i,:,:)));
    xmsr(i,:)=xms(i,:)/sum(xms(i,:));
end
td=1;
tr=0.1;
k = 1;
mmr = [];
while ~(tr>0.99 || td<1e-10),
    lastmidx=[];
    for i=1:size(c,2),
        midx=find(xmsr(i,:)>=project.thresh*td);
        usemasses=union(midx,lastmidx);
        lastmidx=usemasses;
    end
    [d1,d2,d3]=size(xcs);
    mx=zeros(d2,length(usemasses));
    for i=1:size(c,2),
        mx=mx+squeeze(xcs(i,:,usemasses));
    end
    tr=deco_rvcoeff(x,mx);
    mmr = [mmr; tr, size(usemasses,2)];
    k = k+1;
    td=td*0.1;
end
end

function [maxtr,usepeaks]=deco_peak(x,c,s)
global project;
%s=squeeze(s);
s=reshape(s,size(c,2),size(x,2)); % squeeze, do not work properly for 1x1xm array;
xcs=zeros(size(c,2),size(c,1),size(s,2),'single');
maxxcs=zeros(1,size(c,2));
for i=1:size(c,2),
    xcs(i,:,:)=c(:,i)*s(i,:);
    maxxcs(i)=max(max(squeeze(xcs(i,:,:))));
end

rxcs=maxxcs/max(maxxcs);
prxcs=rxcs/sum(rxcs);

[sprxcs,idx]=sort(prxcs,'descend');
cx=zeros(size(x));
tr=zeros(length(idx),1,'double');
invtr=zeros(length(idx),1,'double');
for i=1:length(idx),
    cx=cx+squeeze(xcs(idx(i),:,:));
    tr(i)=deco_rvcoeff(x,cx);
    invtr(i)=deco_rvcoeff(x,squeeze(xcs(idx(i),:,:)));
end

maxtr=max(tr);
[a,usepeaks] = find(prxcs>=project.thresh*100);
% for i=p:-1:length(usepeaks)+1
%     usepeaks(i)=0;
% end

%npeaks=length(usepeaks);

end

function y=deco_rvcoeff(d1,d2)
s=d1*d1';
t=d2*d2';
y=double(trace(s*t))/sqrt(double(trace(s*s))*double(trace(t*t)));
end


% function [npeaks,usemasses,c,s]=try_deco_usemasses(data)
% global project;
% pw = project.pw;
% options = struct('display','off','plots','none');
% npeaks=10;
% residue=1e-6;
% noisy=true;
% notusemasses=1:length(data);
% while noisy==true,
%     %usemasses=trydeco_criticalmass(data,residue);
%     usemasses=1:length(data);
%     notusemasses=setdiff(notusemasses,usemasses);
%     sp = zeros(npeaks,pw*4,'single');
%     delta=(3.5/(npeaks+1));
%     init_pos=floor(((1:npeaks)*delta+0.25)*pw); % space peaks
%     for sp_init = 1:npeaks,
%         sp(sp_init,:) = deco_makegauss(1:pw*4,init_pos(sp_init),pw*1/6); % create initiation gaussian peaks
%     end
%     idx = 1:project.nfiles; % index to the files
%     actfiles=project.nfiles; % count the number of actual files to include
%     s = zeros(actfiles,npeaks,length(usemasses),'single');
%     rs = zeros(actfiles,npeaks,length(notusemasses),'single');
%     c = zeros(4*pw,npeaks,'single');
%     for part=1:length(idx),
%         [c(:,:),s(idx(part),:,:)]=deco_mcr(data((idx(part)-1)*4*pw+1:idx(part)*4*pw,usemasses),sp',options);
%         rs(idx(part),:,:)=c\data((idx(part)-1)*4*pw+1:idx(part)*4*pw,notusemasses);
%     end
%     [evectorc,evaluec] =eig(c(:,:)'*c(:,:));
%     evaluec=sum(evaluec,1);
%     [noisy,npeaks,c,s]=deco_noisecheck(data((idx(part)-1)*4*pw+1:idx(part)*4*pw,notusemasses),c,s,rs);
% end
% 
% end
% 
% 
% function [y,npeaks,c,s]=deco_noisecheck(x,c,s,rs)
% global project;
% s=squeeze(s);
% sc=size(c,2);
% st=size(c,1);
% ss=size(s,2);
% xcs=zeros(sc,st,ss,'single');
% for i=1:sc,
%     xcs(i,:,:)=c(:,i)*s(i,:);
%     maxxcs(i)=max(max(squeeze(xcs(i,:,:))));
% end
% 
% rxcs=maxxcs/max(maxxcs);
% srxcs=rxcs/sum(rxcs);
% 
% [a,usepeaks] = find(srxcs>=project.thresh*100);
% c=c(:,usepeaks);
% s=s(usepeaks,:);
% npeaks=length(usepeaks);
% 
% tcs=(squeeze(median(s,1)))';
% tncs=(squeeze(median(rs,1)))';
% 
% mtcs=sum(tcs,2);
% mtncs=sum(tncs,2);
% 
% % mx=max(x);
% % figure;
% % plot(mx);
% 
% y=false; % for debug purposes; S. Krishnan
% 
% end
