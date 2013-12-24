function [copt, sopt, sdopt, numiter] = deco_als990(block, d, s0, nexp, nit, tolsigma)

% Part 1: Initialisations
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
if nargin<4,nexp=1;end
if nargin<5,nit=50;end
if nargin<6,tolsigma=0.1;end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% Part 1: Unimodality
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
rmod=1.0001; % RHJ - Unimodal constraint tolerance for the conc
cmod=1; % RHJ - Unimodal constraint is horizontal
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% Part 3: Alternating Least Squares Optimization
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
[copt,sopt,sdopt,numiter,]=als99(block, d, s0, nexp, nit, tolsigma, rmod, cmod);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

end

% SUB-ROUITINES
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

function [copt,sopt,sdopt,numiter]=als99(block, d, s0, nexp, nit, tolsigma, rmod, cmod)

global project;
nb = project.last-project.first+1;
nphhh = min(size(s0));
waiter=waitbar(0,'Deconvolution','Name',['Block [' num2str(block) ' of ' num2str(nb) '] np= ' num2str(nphhh) ]);

[nrow, ncol] = size(d);
nsign = size(s0, 1);

% Part: Non-negativity and unimodality
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
[matr, matc, ncinic, ncfin, nrinic, nrfin] = init(d, nexp);
spneg = ones(nsign,matr);
cneg = ones(matc,nsign);
isp = ones(nexp,nsign);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% Part: Reproductiion of original data by PCA
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% dn is the experimental matrix and d is the PCA reproduced matrix
dn=d;
[u,s,v,d]=deco_pcarep(dn,nsign);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% Part: Error control
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
sstn=sum(sum(dn.*dn));
sst=sum(sum(d.*d));
sigma2=sqrt(sstn);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

idev = 0;  % divergence counter
niter = 0;
while niter < nit
    niter=niter+1;
    numiter=niter;
    % ***************************************
    % E) ESTIMATE CONCENTRATIONS (ALS solutions)
    % ***************************************
    conc=d/s0;
    abss = abs(s0);
    conc = abs(conc);

    conc = nonnegconstrainconc(d, conc, abss, matc, nrinic, nrfin, cneg, isp);
    conc = unimodalityconstrainconc(matc, conc, nrinic, nrfin, nsign, rmod, cmod);
    
    % ********************************
    % ESTIMATE SPECTRA (ALS solution)
    % ********************************
    abss=conc\d;
    abss = nonnegconstrainspec(d, conc, abss, matr, ncinic, ncfin, spneg, isp);
    
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
    % *************************************************************
    % If change is positive, the optimization is working correctly
    % *************************************************************
    if change>0 || niter==1,
        sigma2=sigma;
        copt=conc;
        sopt=abss;
        sdopt=sstd;
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
% finish the iterative optimization if maximum number of allowed iterations is exceeded
close(waiter);
return          % 3rd return (end of optimization, number of iterations exceeded)
end

   
function conc = unimodalityconstrainconc(matc, conc, nrinic, nrfin, nsign, rmod, cmod)

% ***********
% unimodality
% ***********

spmod = ones(matc,1)*ones(1,nsign); % RHJ - result of former line

for i = 1:matc
    kinic=nrinic(i);
    kfin=nrfin(i);
    conc2=conc(kinic:kfin,:);
    for ii=1:nsign,
        if spmod(i,ii)==1,
            conc2(:,ii)=deco_unimod(conc2(:,ii),rmod,cmod);
        end
    end
    conc(kinic:kfin,:)=conc2;
end

end

function conc = nonnegconstrainconc(d, conc, abss, matc, nrinic, nrfin, cneg, isp)

% *****************************************************
% CONSTRAIN APPROPRIATELY THE CONCENTRATIONS non neg
% *****************************************************

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

function abss = nonnegconstrainspec(d, conc, abss, matr, ncinic, ncfin, spneg, isp)

% ********************
% non-negative spectra
% ********************
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

function [matr, matc, ncinic, ncfin, nrinic, nrfin] = init(d, nexp)

[nrow, ncol] = size(d);

if nexp==1,
    ncinic(nexp)=1;
    ncfin(nexp)=ncol;
    nrinic(nexp)=1;
    nrfin(nexp)=nrow;
    matr = 1;
    matc = 1;
elseif nexp > 1,
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
end

end
                           