function []= deco_baseline_set(ff,dir,nfiles)
%--------------------------------------------------
% convert the baseline of the files in the set
% use the project settings as the parameters
%--------------------------------------------------
global project;

for i=1:nfiles % process all files in the filesbox    
   name = ff{i};
   filename = [dir name]; % construct full name = directory + filename
   dbcname  = [strtok(filename,'.') '.dbc'];
  
   load('-mat',dbcname,'theMat','options'); % load from previous run
   v = project.basemethod; % the selected baseline correction method      
   if v==2
      %disp('Groningen');
      options.gsize=project.gsize;
      options.gmove=project.gmove;
      options.gthresh=project.gthresh;
      theMat = deco_baseline(theMat,options,name,2,i,nfiles); %
      save(dbcname,'-append','theMat','options'); % store the baseline corrected file
   elseif v==3
      %disp('Frans');
      options.basewin = project.basewin;
      theMat = deco_baseline(theMat,options,name,3,i,nfiles);
      save(dbcname,'-append','theMat','options'); % store the baseline corrected file
   elseif v==1
      %disp('default option: Eilers');         
      options.asym=project.asym;
      options.lambda=project.lambda;
      options.order=project.order;
      theMat = deco_baseline(theMat,options,name,1,i,nfiles); %
      save(dbcname,'-append','theMat','options'); % store the baseline corrected file
   elseif v==4
       %disp('default option: Convex hull');         
      theMat = deco_baseline(theMat_cdf,options,name,4,i,nfiles); %
      save(dbcname,'-append','theMat','options'); % store the baseline corrected file
   end
end
