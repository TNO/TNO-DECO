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

% try code begin, S. Krishnan
[cdata, ncdata]=trydeco_criticalmass(totblock);
totblock=cdata;
% try code end

%block
tic = sum(totblock);
mx = max(tic);
tic = tic / mx;
[a,usemasses] = find(tic>=project.thresh); % which masses do we want to include
% remove exluded masses
[v,p] = intersect(usemasses,project.ex_mass-project.minmass+1);
usemasses(p)=[];

% For comparision purposes between critical and non-critical masses, all
% masses has been used; S. Krishnan
usemasses=1:length(totblock);

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
        stop= 0;
        while stop==0
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
            
            % Re-estimation of 's' with estimated value of 'c' after
            % noise removal; S. Krishnan
            Ntotblock = read_block(block);
            Ntic = sum(Ntotblock);
            Nmx = max(Ntic);
            Ntic = Ntic / Nmx;
            [Na,Nusemasses] = find(Ntic>=project.thresh); % which masses do we want to include
            % remove exluded masses
            [Nv,Np] = intersect(Nusemasses,project.ex_mass-project.minmass+1);
            Nusemasses(Np)=[];
            %rs = zeros(actfiles,npeaks,length(Nusemasses),'single');
            % Re-estimation of 's' with estimated value of 'c' after
            % noise removal; S. Krishnan
            %++++++++++++++++++++++++++++++++++++++
            %Estimating the spectral component 's' in the non-criticalS
            %mass;
            rs = zeros(actfiles,npeaks,length(ncdata),'single');
            %++++++++++++++++++++++++++++++++++++++
            %Estimating the spectral component 's' in the non-critical
            %mass;
            
            for part=1:length(idx), % for each individual file
               % s = spectrum (height per mass peak) (ncompounds,nmass)
               % c = concentration per spectrum (4* pw, ncompounds)
               % sp= estimated concentration profiles
               [c(:,:),s(idx(part),:,:)]=deco_mcr(totblock((idx(part)-1)*4*pw+1:idx(part)*4*pw,usemasses),sp',options);
               % Re-estimation of 's' with estimated value of 'c' after
               % noise removal; S. Krishnan
               %ec=c;
               %[c(:,:),rs(idx(part),:,:)]=deco_mcr(Ntotblock((idx(part)-1)*4*pw+1:idx(part)*4*pw,Nusemasses),ec,options);
               % Re-estimation of 's' with estimated value of 'c' after
               % noise removal; S. Krishnan
               %++++++++++++++++++++++++++++++++++++++
               %Estimating the spectral component 's' in the non-critical
               %mass; The concentration profile 'c' is determined from the
               %critical mass; S. Krishnan
               rs(idx(part),:,:)=c\ncdata;
               %++++++++++++++++++++++++++++++++++++++
            end
            
            % The spectral result is transposed for comparison purposes,
            % very much try; S. Krishnan
            tcs=(squeeze(median(s,1)))';
            tncs=(squeeze(median(rs,1)))';
                       
            % Re-estimation of 's' with estimated value of 'c' after
            % noise removal; S. Krishnan
            %es=(squeeze(median(s,1))); % squeezing the 's' from the first iteration, for debug purpose only; S. Krishnan
            %s=rs;
            %sp=ec';
            usemasses=Nusemasses;
            totblock=Ntotblock;
            % Re-estimation of 's' with estimated value of 'c' after
            % noise removal; S. Krishnan

            %compact the spectrum dimension into a single set of spectra
            s=(squeeze(median(s,1))); % s = height per mass peak                       
            %clear c; % remove the concentrations
            sopt=zeros(npeaks,size(totblock,2),'single');
            copt=zeros(4*pw*numfiles,npeaks,'single');
            areaopt=zeros(npeaks,numfiles);
                       
            %[copt(thesescans,:),sopt(:,usemasses),sdopt,suggest_numpeaks,dummy,ni]=...
            %    als99_nmf(totblock(thesescans,usemasses),s,actfiles,50,sp,c);       
            n = project.maxiter;
            if n<2 
                n=2;
            end
            
            % initialising with the concentration profile 'c' instead of
            % the spectra profile 's', hence next line of the code commented; S. Krishnan
            [copt(thesescans,:),sopt(:,usemasses),sdopt,suggest_numpeaks,dummy,ni]=...
                deco_als99(block,totblock(thesescans,usemasses),c,actfiles,n,sp,c);              
            % initialising with the concentration profile 'c' instead of
            % the spectra profile 's', hence next line of the code commented; S. Krishnan
           
            %[copt(thesescans,:),sopt(:,usemasses),sdopt,suggest_numpeaks,dummy,ni]=...
                %deco_als99(block,totblock(thesescans,usemasses),s,actfiles,n,sp,c);              
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