function [] = deco_resample_set(ff,dir,nfiles)
global project

str='Resample files';
ht = waitbar(0,str);
for i=1:nfiles % process all files in the filesbox    
   waitbar(i/nfiles,ht);   
   name = ff{i};
   filename = [dir name]; % construct full name = directory + filename
   dbcname  = [strtok(filename,'.') '.dbc'];
   deco_resample_file(dbcname,project.interval,project.resample);  
   %load('-mat',dbcname,'theMat','scantimes'); % load from previous run ?
   
end

close(ht);