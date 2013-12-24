function y = deco_readblocktic(name,block,file,pw,ntraces,nmasses)
%-------------------------------------------
% y= deco_readblocktic(name,block,pw,nfiles,ntraces,nmasses)
% read blocks (all files within block limits) from file
% construct fro mthos block the tic spectrum for requested file
%-------------------------------------------
% traces is calculated from pw
fp = fopen(name,'rb');

pw_start = (block-1)*2*pw;
pw_end   = (block)*2*pw+2*pw;
esize = 4; % element size (real * 4) 

fileblock = nmasses*ntraces*esize; % size of all data in a single file
blocksize = (block-1)*2*pw*nmasses*esize; % size of data in a single bock

y = single(zeros((pw_end-pw_start),1));

fseek(fp,double((file-1)*fileblock+blocksize),'bof'); % goto begin of filenumber
n=1;
for i=pw_start:pw_end-1 % scan from start trace to end trace
     idx = (file-1)*4*pw+i-pw_start+1; % get index number
     y(n) = single(sum(fread(fp,double(nmasses),'float'))); % read all masses in trace 
     n=n+1;
end   
fclose(fp); % close file

