function [copt,sopt,sdopt,ropt,areaopt,numiter,sstd]=deco_als99(block,d,x0,nexp,nit,sp,initC,tolsigma,isp)

%-------------------------------------------------------------------
%-------------------------------------------------------------------
%** Multivariate Curve Resolution (MCR) - Alternating Least Squares (ALS) *********************
%
%function
%[copt,sopt,sdopt,ropt,areaopt,rtopt]=als(d,x0,nexp,nit,tolsigma,isp,csel,ssel,vclos1,vclos2);
%      where
% INPUT VALUES:
%       d : experimental data matrix
%       x0: initial estimates of the concentration profiles 
%           or of the species spectra 
%     nexp: number of data matrices analyzed simultaneously
%      nit: maximum number of iterations (50 is the default)
% tolsigma: convergence criterion in the difference of sd of residuals 
%           between iterations (0.1% is the default)
%      isp: correspondence among the species in the experiments
% OUTPUT VALUES:
%       copt: optimized species concentrations
%       sopt: optimized species spectra
%       ropt: residuals d - copt*sopt at the optimum
%      sdopt: standard deviation of fitting residuals at the optimum
%    areaopt: areas of concentration profiles (only for quantitation) 	
%      rtopt: ratio of areas (only for quantitation) 
%
%function [copt,sopt,sdopt,ropt,areaopt,rtopt]=als(d,x0,nexp,nit,tolsigma,isp,csel,ssel,vclos1,vclos2);
%
%*****************************************************************************

%*****************************************************************************
% other important variables
% nrow number of rows (spectra) in d
% ncol number of columns (channels, wavelengths) in d
% ils kind of initial estimate provided from efa;
%     ils = 1 initial estimates of concentrations;
%     ils = 2 initial estimates of spectra
% nsign is the total number of significant species
% nexp number of experiments simultaneously analyzed
% nspec number of species in each experiment
% ishape = 0,1,2 data structure (see below)

% Part 1: Check dimensions of initial extimates
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

[nsign, nrow, ncol, abss] = diminitestimates(d, x0);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!% 

% Part 2: Initialisations
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
if nargin<4,nexp=1;end
if isempty(nexp)|| nexp==0,nexp=1;end
if nargin<5,nit=50;end 
if isempty(nit)|| nit==0,nit=50;end
if nargin<8,tolsigma=0.1;end
if isempty(tolsigma) || tolsigma==0,tolsigma=0.1;end
if nargin<9 && nexp==1,isp=ones(1,nsign);end
if nargin<9 && nexp>1,isp=ones(nexp,nsign);end
if isempty(isp) | isp==0
    isp=ones(1,nsign);
end
if nexp==1,
    ncinic(nexp)=1;
    ncfin(nexp)=ncol;
    nrinic(nexp)=1;
    nrfin(nexp)=nrow;
end

scons='y'; % all the spectra matrices the same constraints
ccons='y'; % all the concentration matrices the same constraints
niter=0;%#ok<NASGU> % iterations counter
idev=0;% divergence counter
answer='n'; % default answer
ineg=0;%#ok<NASGU> % used for non-negativity constraints
imod=0;% used for unimodality constraints
matr=1;%#ok<NASGU> % one matrix
matc=1;%#ok<NASGU> % one matrix
numiter=0;
%***************************
% DEFINITION OF THE DATA SET
%***************************
totalconc(1:nsign,1:nexp)=ones(nsign,nexp);
% IN SIMULTANEOUS ANALYSIS OF SEVERAL SAMPLES
% ENTER NUMBER OF SPECTRA
% RHJ - only one matrix is used so this part can be bypassed

if nexp > 1,
   matr = 1;
   ncinic(1)=1;
   ncfin(1)=ncol;
   matc = nexp;		
   nrinic(1)=1;
   nrsol = ones(1,matc)*size(d,1)/matc;    
   for i=1:matc,            
      nrfin(1)=nrsol(1);
        if i>1,
          nrinic(i)=nrfin(i-1)+1;
          nrfin(i)=nrinic(i)+nrsol(i)-1;
        end
        ncinic(i)=1;
        ncfin(i)=ncol;
   end 
else
    %----------------------------------------------------------
    % WHEN ONLY ONE EXPERIMENT IS PRESENT EVERYTHING IS CLEAR
    %----------------------------------------------------------
    nrsol(1)=nrow;
    nrinic(1)=1;
    nrfin(1)=nrsol(1);
    matr = 1;
    matc = 1;
    isp(1,1:nsign)=ones(1,nsign);
    ishape=0;
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%


% Part 3: Setting constrains
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
wcons = [1 2]; % RHJ - non-negativity and unimodality;
ineg=2; % RHJ - Both concentrations and spectra should be positive;
while answer == 'n' ||answer =='N'
    % **************************
    % NON-NEGATIVITY CONSTRAINTS
    % **************************   
    c1 = find(wcons == 1);
    if ~isempty(c1) 
         if ineg==3||ineg ==2        
             if scons=='y' || nexp == 1
                    spneg = ones(nsign,matr);
            end
        end    
        if ineg==1||ineg ==2           
            if ccons=='y' || nexp == 1             
                    cneg = ones(matc,nsign);
            end           
        end
    else       
        cneg=zeros(matc,nsign);
        spneg = zeros(nsign,matr);
    end 
    % **********************
    % UNIMODALITY CONSTRAINT
    % **********************  
    imod=1; % RHJ - unimodality only on concentration profiles
    spmod = ones(matc,1)*ones(1,nsign); % RHJ - result of former line
    rmod=1.0001; % RHJ - Unimodal constraint tolerance for the conc
    cmod=1; % RHJ - Unimodal constraint is horizontal 
    % *******************
    % THREE-WAY STRUCTURE
    % *******************
    c7 = find(wcons == 7);
    if ~isempty(c7) || nexp > 1 % wordt echt gebruikt
        ishape=2;
    end
    % **************************
    % CHECKING THE INITIAL INPUT
    % **************************
     answer='y'; % RHJ - all data correct;
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%


% Part 3: Alternating Least Squares Optimization
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
[copt,sopt,sdopt,ropt,areaopt,numiter,sstd]=als99(block, d, x0, nrow, ncol, nit, tolsigma, cneg, spneg, c1, matc, nrinic, nrfin, nsign, abss, totalconc, ishape, imod, spmod, rmod, cmod, ineg, matr, ncinic, ncfin, isp);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%


end



function [nsign, nrow, ncol, abss] = diminitestimates(d, x0)
    
[nrow,ncol]=size(d);
% check dimensions of initial estimates
[nrow2,ncol2]=size(x0);
if nrow2==nrow && nrow ~=ncol,
    nsign=ncol2;
    ils=1;
end
if ncol2==nrow && ncol~=ncol2
    nsign=nrow2;
    x0=x0';
    ils=1;
end
if ncol2==ncol,
    nsign=nrow2;
    ils=2;
end
if nrow2==ncol && nrow ~= ncol,
    nsign=ncol2;
    x0=x0';
    ils=2;
end

if ils==1,
    conc=x0;
    [nrow,nsign]=size(conc);
    abss=conc\d;
end

[nrow,ncol]=size(d);

if ils==2,
    abss=x0;
    [nsign,ncol]=size(abss);
    if ncol==1,
        conc=d/abss';
    else
        conc=d/abss;
    end
end

end


function [copt,sopt,sdopt,ropt,areaopt,numiter,sstd]=als99(block, d, x0, nrow, ncol, nit, tolsigma, cneg, spneg, c1, matc, nrinic, nrfin, nsign, abss, totalconc, ishape, imod, spmod, rmod, cmod, ineg, matr, ncinic, ncfin, isp)

global project;
nb = project.last-project.first+1;
nphhh = min(size(x0));

waiter=waitbar(0,'Deconvolution','Name',['Block [' num2str(block) ' of ' num2str(nb) '] np= ' num2str(nphhh) ]);

% Part 3.1: Reproductiion of original data by PCA
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% dn is the experimental matrix and d is the PCA reproduced matrix
dn=d;
[u,s,v,d]=deco_pcarep(dn,nsign);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

sstn=sum(sum(dn.*dn));
sst=sum(sum(d.*d));
sigma2=sqrt(sstn);
sigma = sqrt(sst);

niter = 0;
while niter < nit
    niter=niter+1;
    numiter=niter;
    % ***************************************
    % E) ESTIMATE CONCENTRATIONS (ALS solutions)
    % ***************************************
    conc=d/abss;
    abss = abs(abss);
    conc = abs(conc);
    % *****************************************************
    % CONSTRAIN APPROPRIATELY THE CONCENTRATIONS non neg
    % *****************************************************
    if ~isempty(c1)
        for i =1:matc
            kinic=nrinic(i);
            kfin=nrfin(i);
            conc2=conc(kinic:kfin,:);
            for j=kinic:kfin
                if cneg(i,:) == ones(1,size(isp,2))
                    x=deco_fnnls(abss*abss',abss*d(j,:)');
                    conc2(j-kinic+1,:)=x';
                end
            end
            conc(kinic:kfin,:) = conc2;
        end
    end
    if ishape==0 || niter==1,
        for j=1:nsign,
            for inexp=1:matc,
                totalconc(j,inexp)=sum(conc(nrinic(inexp):nrfin(inexp),j));
            end
            if totalconc(j,1)>0,
                rt(j,1:matc)=totalconc(j,1:matc)./totalconc(j,1);
            else
                rt(j,1:matc)=totalconc(j,1:matc);
            end
        end
    end
    % ***********
    % unimodality
    % ***********
    for i = 1:matc
        kinic=nrinic(i);
        kfin=nrfin(i);
        conc2=conc(kinic:kfin,:);
        if imod==1||imod==3,
            for ii=1:nsign,
                if spmod(i,ii)==1,
                    conc2(:,ii)=deco_unimod(conc2(:,ii),rmod,cmod);
                end
            end
        end
        conc(kinic:kfin,:)=conc2;
    end
    % ************************************************
    % QUANTITATIVE INFORMATION FOR THREE-WAY DATA SETS
    % ************************************************
    % recalculation of total and ratio concentrations if ishape=0 or niter=1
    if ishape==0 || niter==1,
        for j=1:nsign,
            for inexp=1:matc,
                totalconc(j,inexp)=sum(conc(nrinic(inexp):nrfin(inexp),j));
            end
            if totalconc(j,1)>0,
                rt(j,1:matc)=totalconc(j,1:matc)./totalconc(j,1);
            else
                rt(j,1:matc)=totalconc(j,1:matc);
            end
        end
    end
    % areas under concentration profiles
    area=totalconc;
    % ********************************
    % ESTIMATE SPECTRA (ALS solution)
    % ********************************
    abss=conc\d;
    % ********************
    % non-negative spectra
    % ********************
    if ineg ==2 ||ineg==3,
        for i = 1:matr
            kinic = ncinic(i);
            kfin = ncfin(i);
            abss2 = abss(:,kinic:kfin);
            for j=kinic:kfin,
                if spneg(:,i)== ones(size(isp,2),1)
                    abss2(:,j-kinic+1)=deco_fnnls(conc'*conc,conc'*d(:,j));
                end
            end
            abss(:,kinic:kfin)=abss2;
        end
    end
    % ************************************
    % constrain the unimodality of spectra
    % ************************************
    for i = 1:matr
        kinic = ncinic(i);
        kfin = ncfin(i);
        abss2 = abss(:,kinic:kfin);
        if imod==2||imod==3,
            for j=1:nsign,
                if spsmod(i,j)==1
                    dummy=deco_unimod(abss2(j,:)',smod,cmod);
                    abss2(j,:)=dummy';
                end
            end
        end
        abss(:,kinic:kfin)=abss2;
    end
    % *******************************
    % CALCULATE RESIDUALS
    % *******************************
    res=d-conc*abss;
    u=sum(sum(res.*res));
    res=[]; % memory saving
    resn=dn-conc*abss;
    un=sum(sum(resn.*resn));
    clear resn; % memory saving
    % ********************************
    % OPTIMIZATION RESULTS
    % *********************************
    sigma=sqrt(u/(nrow*ncol));
    change=((sigma2-sigma)/sigma);
    if change < 0.0,
        idev=idev+1;% disp('FITING IS NOT IMPROVING !!!')
    else
        idev=0;% disp('FITING IS IMPROVING !!!')
    end
    change=change*100; % change in sigma
    waitbar(niter/nit,waiter,['change in sigma:' num2str(change)]);
    sstd(1)=sqrt(u/sst)*100; % lack of fit in PCA
    sstd(2)=sqrt(un/sstn)*100; % lack of fit in exp
    r2=(sstn-un)/sstn; % percent variance explained
    % *************************************************************
    % If change is positive, the optimization is working correctly
    % *************************************************************
    if change>0 || niter==1,
        sigma2=sigma;
        copt=conc;
        sopt=abss;
        sdopt=sstd;
        ropt=res;
        rtopt=[]; %rt';
        itopt=niter;
        areaopt=area;
        r2opt=r2;
    end
    % ******************************************************************
    % test for convergence within maximum number of iterations allowed
    % ******************************************************************
    if abs(change) < tolsigma,
        close(waiter);
        return   % 1st return (end of the optimization, convergence)
    end
    %  finish the iterative optimization if divergence occurs 20 times consecutively
    if idev > 20,
        close(waiter);
        return  % 2nd return (end of optimization, divergence)
    end
    % this end refers to number of iterations initially proposed exceeded
end %%  while loop
close(waiter);
% finish the iterative optimization if maximum number of allowed iterations is exceeded
return          % 3rd return (end of optimization, number of iterations exceeded)

end







    