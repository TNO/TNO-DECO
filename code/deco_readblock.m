function y = deco_readblock(name,block,pw,nfiles,ntraces,nmasses)
%-------------------------------------------
% y= deco_readblock(name,block,pw,nfiles,ntraces,nmasses)
% read block (all files within block limits) from file
%-------------------------------------------
% traces is calculated from pw
fp = fopen(name,'rb');
if fp<0
    msg = ['could not open file: ' name];
    Errordlg(msg);
    y=0;
    return;
end
pw_start = (block-1)*2*pw;
pw_end   = (block)*2*pw+2*pw;
esize = 4; % element size (real * 4) 
%nf=nfiles;
%nm=nmasses;
%nt = ntraces;

fileblock = nmasses*ntraces*esize; % size of all data in a single file
blocksize = (block-1)*2*pw*nmasses*esize; % size of data in a single bock
y = single(zeros((pw_end-pw_start)*nfiles,nmasses));
size(y);

for file = 1: nfiles % read all files
    fseek(fp,double((file-1)*fileblock+blocksize),'bof'); % goto begin of filenumber
    for i=pw_start:pw_end-1 % scan from start trace to end trace
        idx = (file-1)*4*pw+i-pw_start+1; % get index number
        z =  fread(fp,double(nmasses),'float'); % read all masses in trace
        y(idx,:) = z;
    end   
end

fclose(fp); % close file