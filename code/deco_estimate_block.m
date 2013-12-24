function [numsig,npeak,values,s_sim,ie,ind,snsig] = deco_estimate_block(blocknum,level,maxfac)
global project
%------------------------------------------------------------
% deco function to estimate and analyse the number of components per block
% input 
%    blocknum : the block number
% returns 
%  numsig number of significant peaks estimated from rsd
%  npeak  number of peaks estimated from peakpick
%  values eigen values from real data
%  sim    eigen values from simulation data
% written bij J.T.W.E. Vogels 25/03/2008 (last edit)
%----------------------------------------------------------
   pw = project.pw;
   
   nout = nargout;
   block=deco_readblock([project.name '.bin'],blocknum,pw,project.nfiles,project.ntraces,project.nmasses);
  
   %if (project.use_iq) % do we want to filter using iq values
   %      block(:,iq_remove) = []; % remove iq selected masses
   %end   
     
   if nout>6 % use only masses awhich have decent signal2 noise ratio
       x=(1:4*project.pw)';
       sn = zeros(size(block,2),1);
       blsn = block;
       for j=1:size(block,2)
           tot = [];
           for i=1:project.nfiles
               sp = block((i-1)*pw*4+1:i*pw*4,j);
               K = polyfit(x,sp,2); % find (simple) second order baseline
               yout = sp - polyval(K,x); % baseline corrected data
               a = deco_signal2noise(yout);
               if (a<10) % signal2 noise is too low
                   blsn((i-1)*pw*4+1:i*pw*4,j) = 0; % set data for these masses to zero
               end
           end
       end
       blsn = deco_mncn(blsn); % datablock after removing local baselines  
       ssq = sum(sum((deco_mncn(blsn)).^2));
       simblsn=randn(size(block))*(ssq/numel(blsn)).^(1/2);

   end
       
   %disp('deco estimate block');
      
   ssq = sum(sum((deco_mncn(block)).^2));
   simblock=randn(size(block))*(ssq/numel(block)).^(1/2);
       
   [a,b]=size(block);
   if a > b
       block2 = block'*block;
       sblock = simblock'*simblock;
       if nout > 6
           b2 = blsn'*blsn;
           b3 = simblsn'*simblsn;
       end
   else
       block2 = block*block';
       sblock = simblock*simblock';
       if nout>6
           b2 = blsn*blsn';
           b3 = simblsn*simblsn';
       end
   end
       
   % METHOD using eig is much faster (factor 2.5) than using svd(x,0)
   values =(sort(real(abs(eig(block2))),'descend'));
   % do eig on the estimated (simulated) noise data
   s_sim = (sort(real(abs(eig(sblock))),'descend'));  
   snsig = -1;
   if nout > 6
       g_sim  = (sort(real(abs(eig(b2))),'descend'));   % eigenvalues of block baseline corrected data 
       g_sim1 = (sort(real(abs(eig(b3))),'descend'));   % eigenvalues of block baseline corrected data 
       snsig=1;
       while (g_sim(snsig) > g_sim1(snsig)*level)
           snsig=snsig+1;
       end
   end
     
   %estimate the number of peaks for eigenvalues above noise eigenvalue
   numsig=1;
   while (values(numsig) > s_sim(numsig)*level)
        numsig=numsig+1;
   end
   
   numsig;
   snsig;

   if nout > 6
       numsig1=1;
       while (g_sim1(numsig1) > s_sim(numsig1)*level)
            numsig1=numsig1+1;
       end
   end
   
   g = sum(values) - cumsum(values);
   ms = min(a,b);
   
   if maxfac<=0 || maxfac > ms
       maxfac = ms;
   end
       
   %--------------------------------------------
   % Eigenvalue analysis for information 
   % see Brereton chapter 8 page 
   %--------------------------------------------
   % calc some preliminaries to eigenvalue analysis
   sa = zeros(ms,1);
   for i=1:ms
       sa(i)=(a-i+1)*(b-i+1);
   end 
   fn = sum(sa)-cumsum(sa);
   % done
   
   % calculate rsd
   ms= min(min(a,b),maxfac);
   rsd =zeros(ms,1);
   for i=1:ms
       if (b-i>0)
           rsd(i) = sqrt(g(i)/(a*(b-i)));
       end
   end
    %imbedded error
   ie=zeros(ms,1);
   for i=1:ms
       if (b-i)>0 && a>0
           ie(i)=sqrt((g(i)*i)/(a*b*(b-i)));
       end
   end

   % calculate indicator function
   ind = zeros(ms,1);
   for i=1:ms
       if (b-i)>0
           ind(i) = rsd(i)/((b-i)*(b-i));
       end
   end
    
   values = values(1:ms,1);
   s_sim =s_sim(1:ms,1);
      
   %---------------------------------------
   %estimate number of peaks in block
   %---------------------------------------
   k=[]; % empty peakpositions
   for j=1:project.nfiles
        b = block((j-1)*4*pw+1:j*4*pw,:);
        b = sum(b,2);
        np = deco_peakpick(b,-1,0,8,true,8);
        k = [k np];
    end
    k =unique(k); 
    npeak = sum(diff(k)>1)+1;
    