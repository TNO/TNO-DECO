%================================================================
function np = deco_process_block(block,npeaks,recalc,autoreduce)
%================================================================
% written by R. Jellema, 
% extracted from deco_run for single runs 
% by J.T.W.E. Vogels 7-7-2006
% 23/11/2008 added autoreduce option
%================================wai====================
global project;
pw = project.pw;
numfiles = project.nfiles;
options = struct('display','off','plots','none');
totblock = read_block(block); % read the data

% peak-picking using f_L1_L0 functions
% kappa0 = 0.001;
% maxitr = 20;
% for j = 1:maxitr
%     kappa(j) = kappa0;
%     [peakL0(:,j), npeakL0(j)] = deco_eilersL0(totblock,kappa0);
%     kappa0 = kappa0 + (1/maxitr)*0.01;
% end
% peak-picking using f_L1_L0 functions

% % peak-finding using ALS
% y0=totblock;
% y=y0/max(max(y0));
% [npeaks,usemasses,c,s]=try_deco_usemasses(y);
% % peak-finding using ALS

%block
tic = sum(totblock);
mx = max(tic);
tic = tic / mx;
[a,usemasses] = find(tic>=project.thresh); % which masses do we want to include
project.usemasses(block)=size(usemasses,2);
% remove exluded masses
[v,p] = intersect(usemasses,project.ex_mass-project.minmass+1);
usemasses(p)=[];

% % 2d-peak picking
% y0 = sum(totblock,2);
% y = y0/max(y0); % Normalisation
% [peakL0, npeakL0, C] = deco_eilersL0(y,0.003);
% [tt, tm] = size(totblock);
% mpk = zeros(tt, tm);
% mpkn = zeros(1,tm);
% cb = totblock;
% for j=1:25
%     [omv, omi] = sort(sum(cb),'descend');
%     rm = omi(1);
%     eic = (C * peakL0)*max(y0);
%     err = abs(cb(:,rm) - eic);
%     cb(:,rm) = 0;
%     y = sum(cb,2)/max(y0);
%     [peakL0, npeakL0, C] = deco_eilersL0(y,0.003);
%     [spk,spkx] = sort(peakL0, 'descend');
%     mpk(spkx(1:npeakL0),rm) = peakL0(spkx(1:npeakL0),:);
%     mpkn(rm) = npeakL0;
% end
% % 2d-peak picking

if npeaks ==0    
    project.deco{block}= [] ;
    return;
end

decoresults = project.deco{block};   % get decoresults struct         
%spikes=find(sum(diff(totblock(:,usemasses)))~=0);   % remove spikes 

% changed the condition from greater than 1 to notempty; S. Krishnan
%if length(usemasses)>1 % at least 20 masses should be in the matrix 
if ~isempty(usemasses)
    if  recalc || numpeaks(block)~=size(decoresults.sopt,1) || recalc(block),    
        est=0;
        stop= 0;
        while stop==0
           
            
            %+++++++
            % try code has to be here but have to check!!; S. Krishnan
            %+++++++
            
            
            sp = zeros(npeaks,pw*4,'single');
            delta=(3.5/(npeaks+1));
            init_pos=floor(([1:npeaks]*delta+0.25)*pw); % space peaks
            for sp_init = 1:npeaks,
                sp(sp_init,:) = deco_makegauss(1:pw*4,init_pos(sp_init),pw*1/6); % create initiation gaussian peaks
            end
            idx = 1:project.nfiles; % index to the files
             actfiles=project.nfiles; % count the number of actual files to include
             thesescans=1:pw*4*numfiles;
            s = zeros(actfiles,npeaks,length(usemasses),'single');
            c = zeros(4*pw,npeaks,'single');
                       
            for part=1:length(idx), % for each individual file
               % s = spectrum (height per mass peak) (ncompounds,nmass)
               % c = concentration per spectrum (4* pw, ncompounds)
               % sp= estimated concentration profiles
               [c(:,:),s(idx(part),:,:)]=deco_mcr(totblock((idx(part)-1)*4*pw+1:idx(part)*4*pw,usemasses),sp',options);
            end
            
            % try code begin, S. Krishnan
%             if est==0
%                 [npeaks,usemasses,c,s]=try_deco_usemasses(totblock);
%                 project.numpeaks(block)=npeaks;
%                 % storing results of the number of masses used for de-convolution
%                 % after implementing 'try_deco_usemasses' function. S. Krishnan
%                 project.usemasses(block)=size(usemasses,2);
%                 est=1;
%             end
            % try code end
            
            %compact the spectrum dimension into a single set of spectra
            s=(squeeze(median(s,1))); % s = height per mass peak
            %clear c; % remove the concentrations
            
            sopt=zeros(npeaks,size(totblock,2),'single');
            copt=zeros(4*pw*numfiles,npeaks,'single');
            areaopt=zeros(npeaks,numfiles);
                       
%             [copt(thesescans,:),sopt(:,usemasses),sdopt,suggest_numpeaks,dummy,ni]=...
%                als99_nmf(totblock(thesescans,usemasses),s,actfiles,50,sp,c);       
            n = project.maxiter;
            if n<2 
                n=2;
            end
            
%             %changing paramter: 's' to 'c'; S. Krishnan, use transpose of c
%             %in the function, else ils value will not be assigned.
                           
%             [copt(thesescans,:),sopt(:,usemasses),sdopt,suggest_numpeaks,dummy,ni]=...
%                 deco_als99(block,totblock(thesescans,usemasses),c',actfiles,n,sp,c');

% commented for testing lastest deconvolution procedures

%             [copt(thesescans,:),sopt(:,usemasses),sdopt,ni]=...
%                 deco_als990(block,totblock(thesescans,usemasses),s,actfiles,n); 


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%        
        [ec, es, sdopt,ni] = deco_fixspectra(block, totblock(thesescans,usemasses), npeaks, actfiles, usemasses);
  
        sopt=zeros(npeaks,size(totblock,2),'single');
        copt=zeros(4*pw*numfiles,npeaks,'single');
        areaopt=zeros(npeaks,numfiles);
        
        sopt(:,usemasses) = es;
        copt(thesescans,:) = ec;
               
        s = sum(sopt,2); % recalculate
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

%             [copt(thesescans,:),sopt(:,usemasses),sdopt,suggest_numpeaks,dummy,ni]=...
%                 deco_als99(block,totblock(thesescans,usemasses),s,actfiles,n,sp,c); 
            s = sum(sopt,2); % calculate tic 
              
            stop=1;
            if autoreduce>0
                for i=1:npeaks-1
                    f1 = find(sopt(i,:)>0.05);
                    for j=i+1:npeaks
                      f2 = find(sopt(j,:)>0.05);
                      fi = intersect(f1,f2);
                       if size(fi,2)==size(f2,2) | size(fi,2)==size(f2,2)
                           %sprintf('Block %d (%d %d)[%d %d %d]\n',block,i,j,size(f1,2),size(f2,2),size(fi,2))
                           %disp(['equal compounds' num2str(i) ' ' num2str(j)]);
                           stop=0;
                       end
                    end
                end
            end
            if stop==0
                npeaks = npeaks - 1;
            end
            %autoreduce
        end
        quality = [];
        pn = [];
        pmax = zeros(npeaks,project.nfiles);
        for peak=1:npeaks,
             for file=1:numfiles
                 tmpConc = copt((file-1)*4*pw+1:file*4*pw,peak);   
                 warning off Matlab:divideByZero
                 [eV,sV] = deco_gaussfit(1:4*pw,tmpConc);  
                 warning on Matlab:divideByZero
                 %calculate the area of the peak
                 areaopt(peak,file)=sum(copt((file-1)*4*pw+1:file*4*pw,peak))*s(peak); % calculate the area under the peak
                 [mx,mp] = max(tmpConc);
                 if areaopt(peak,file)>0
                    quality(peak,file).err = (eV/areaopt(peak,file))*100; % percent error value of peak fit
                 else
                    quality(peak,file).err=0;
                 end
                 quality(peak,file).symm = sV; % peak symmetry
                 quality(peak,file).mxpos  = mp; % position of maximum
                 pmax(peak,file) = mp;
             end
            pn{peak} = sprintf('pk(%d)%d',block,peak);
        end     
        incl = ones(npeaks,1); % which peaks are to be included in global struct 
                     
        decoresults=struct('copt',copt,'sopt',sopt,'sdopt',sdopt,'areaopt',areaopt,'ni',ni,'usemass',usemasses,'inc',incl,'quality',quality);
        % store the decoresults
        project.deco{block}=decoresults; 
        project.deco{block}.pmax = pmax; % store the peak maxima
        project.deco{block}.pnames = pn;
        np = npeaks;
   else
        np=0;
    end
% added to properly terminate in case there are no masses; S. Krishnan    
else
    project.deco{block}= [] ;
    return;
end
    
       
%==============================
function y = read_block(i)
%==============================
global project
nf = project.nfiles;
y = deco_readblock([project.name '.bin'],i,project.pw,nf,project.ntraces,project.nmasses);