function y= deco_prepare()
%----------------------------------
% prepare the structure for analysing
% if building then
% 1) remove baseline
% 2) resample if necessary
% 3) store data in binfile
% 4) estimate number of peaks
% 5) store project data
%--------------------------------------
global project;

% build blocks from the matlab files extracted from the cdf files
% build the matlab files first with
%disp('BUILD PROJECT')

if project.build,    
    % do alignment if necessary 
    % do baseline correction  
    % copy theMat_cdf files to theMat
    %----------------------------------------------------
    % Align the spectra starting from theMat data
    % store it in theMat
    %----------------------------------------------------
    ht = waitbar(0,'loading data');
    for i=1:project.nfiles
      waitbar(i/(project.nfiles+project.ntestfiles),ht);
      dbcname = getname(i);
      load('-mat',dbcname,'align','theMat_cdf','scantimes_cdf'); % save the alignment factor
      if (align ~=0)       
        theMat = circshift(theMat_cdf,align); 
      else
        theMat = theMat_cdf;
      end
      scantimes = scantimes_cdf;
      save(dbcname,'-append','theMat','scantimes');
    end
    close(ht);        
    %---------------------------------------------------------------
    % resample the set if necessary from theMat store in the Mat
    %---------------------------------------------------------------
    deco_resample_set(project.files,project.sdir,project.nfiles); % resample if necessay
    % if this was necessary then the first and last block must be updated
    %--------------------------------------------------------------
    % Correct baseline from the theMat structure store in theMat
    %--------------------------------------------------------------
    deco_baseline_set(project.files,project.sdir,project.nfiles); % do
    %baseline correction if necessary
    %------------------------------------------------
    %deco_iq_filter; (geeft soms nog een foutmelding controleren JV 28-02-2008)
    %--------------------------------------------------
    startmass=0;
    % if necessary remove spikes
    % to add deco_removespike_set
    %-----------------------------------------------------------
    % determine the masses and traces present in the whole set
    %----------------------------------------------------------
    chrom_length = zeros(project.nfiles,1);
    mymass_start = zeros(project.nfiles,1);
    mymass_end   = zeros(project.nfiles,1);
    rt_start     = zeros(project.nfiles,1);
    rt_end       = zeros(project.nfiles ,1);
   
    for i=1:project.nfiles
         name            = getname(i);
         a               = whos('-file',name,'theMat');
         load('-mat',name,'startmass','scantimes');
         chrom_length(i) = a.size(1);
         mymass_start(i) = startmass;
         mymass_end(i)   = startmass+a.size(2)-1;
         rt_start(i)     = scantimes(1);
         rt_end(i)       = max(scantimes);
    end
    
     project.rt_start   = rt_start;
     project.start_mass = mymass_start;
     project.end_mass   = mymass_end;
     project.chrom_len  = chrom_length;
     overall_minmass    = max(mymass_start);
     project.minmass    = overall_minmass;  
     overall_maxmass    = min(mymass_end);
     project.maxmass    = overall_maxmass;
     project.rt_end     = rt_end;
     
    
     start = mean(project.rt_start);
     project.numblocks = floor((min(rt_end) - max(start))/(project.interval*project.pw*2)-1);
                           
     if std(mymass_start)~=0,
         disp('--------------------------------------------------------')
         disp('first mass is not the same for all chroms')
         disp(['minimum mass in all files: ',num2str(overall_minmass)]);
         disp('--------------------------------------------------------')
     end
     
     if std(mymass_end)~=0,
         disp('--------------------------------------------------------')
         disp('last mass is not the same for all chroms')
         disp(['maximum mass in all files: ',num2str(overall_maxmass)]);
         disp('--------------------------------------------------------')
     end
     name=getname(1);
     load('-mat',name,'startmass','theMat');
     % include only those masses that exist in all files
     startmass=max(mymass_start);
     project.theMasses = [1:min(mymass_end)-startmass+1]; %#ok<NBRAK> % store include masses
     %------------------------------------
     % store the number of data blocks
     %------------------------------------
     numfiles =length(project.files); % the number of files
     n1 = overall_minmass-mymass_start(1)+1;
     n2 = n1 + (overall_maxmass-overall_minmass);        
     % calculate and store number of masses
     project.nmasses = n2-n1+1;
     % calculate and store number of traces
     project.ntraces = min(chrom_length);
     
     while ((project.numblocks+1)*project.pw*2>project.ntraces)
        project.numblocks = project.numblocks -1; 
     end
        
     project.iq = 0;
    %----------------------------------------------
    % write to a large binairy file
    %----------------------------------------------
    fp = fopen([project.name '.bin'],'w+b'); 
    h = waitbar(0,'creating binfile');
    for files = 1:numfiles
        waitbar(files/numfiles,h);
          load('-mat',getname(files),'theMat','resample');
          n1 = overall_minmass-mymass_start(files)+1; % minmass
          n2 = n1 + (overall_maxmass-overall_minmass); % maxmass
          t= theMat(1:project.ntraces,n1:n2);  
          fwrite(fp,t','float'); % write transpose to maintain matlab orientation
    end
    fclose(fp);
    close(h);   
    %-------------------------------------------------------------------
    % new estimator for numpeaks based upon sum of squares
    % do this for all blocks(independend of wether this is necessary
    %-------------------------------------------------------------------
    project.numpeaks=[];
    project.level
    if size(project.numpeaks,2)==0, % if project.numpeaks is empty
        nb= project.numblocks;
        project.deco=cell(nb,1); % remove all previous deco results
        project.recalc=ones(nb,1); % set all peaks to be recalculated   
        numpeaks = zeros(nb,1);
        h=waitbar(0,'estimate peaks');
        for i= 1:nb % do for all blocks
           waitbar((i)/(nb),h);
           
           numpeaks(i)= deco_estimate_block(i,project.level,0); % the 2 is an option 
           if (project.maxn>0 && numpeaks(i)>project.maxn)
               disp(['reduce number of compounds for block' num2str(i)]);
               numpeaks(i)= project.maxn;
           end
        end
        close(h);
        project.numpeaks = numpeaks;
        project.build=0; % do not build again   
    end
end

%=====================================================
% PROCESS THE TESTFILES (IF PRESENT)
%=====================================================
if project.ntestfiles> 0
    nact_test = 0; % number of acttual test files
    h = waitbar(0,'aligning spectra');
    
    for i=1:project.ntestfiles
      waitbar(i/project.ntestfiles,h);
      dbcname = lower([project.sdir '\\' project.testfiles{i}]);
      [a,b] = fileparts(dbcname);
      dbcname = [a b '.dbc'];
      load('-mat',dbcname,'align','theMat_cdf','scantimes_cdf'); % save the alignment factor
      if (align ~=0)       
        theMat = circshift(theMat_cdf,align); 
      else
        theMat = theMat_cdf;
      end
      scantimes = scantimes_cdf; %#ok<NASGU>
      save(dbcname,'-append','theMat','scantimes');
      project.rt_start(project.nfiles + i) = scantimes(1);
      %rt_start(i)     = scantimes(1);
    end
    close(h)
    % calc resample resample if necessay
    deco_resample_set(project.testfiles,project.testdir,project.ntestfiles);
    % set the baselines
    deco_baseline_set(project.testfiles,project.testdir,project.ntestfiles);
    %---------------------------------------------------------
    % this is all that is necessary to process the testfiles
    % check for presence of correct min and max masses and 
    % correct ntraces (if wrong then eliminate file)    
    %---------------------------------------------------------
    fp = fopen([project.name '.tst'],'w+b'); 
    n1 = overall_minmass-mymass_start(files)+1; % minmass
    n2 = n1 + (overall_maxmass-overall_minmass); % maxmass
    h = waitbar(0,'creating test binfile');
    incfiles = [];
    for i=1:project.ntestfiles
        filename = [project.testdir project.testfiles{i}];
        dbcname = [strtok(filename,'.') '.dbc'];     
        load('-mat',dbcname,'startmass','theMat');   % get startmass from the dbc file
        minmass = startmass;
        endmass = startmass + size(theMat,2)-1;
        skip=0;
        if (minmass> project.minmass || endmass < project.maxmass)
            disp(['exclude file: ' project.testfiles{i} ' mass out or range'])
            skip=1;
        end
        if skip==0
            [r,c]= size(theMat);
            nact_test = nact_test+1;
            waitbar(i/project.ntestfiles,h);
            if (n2>c)
                theMat = [theMat  zeros(r,n2-c)]; %#ok<AGROW>
            end
            if (project.ntraces>r)
                theMat = [theMat ; zeros(project.ntraces-r,c)];
            end
            t= theMat(1:project.ntraces,n1:n2);  
            incfiles{nact_test} = project.testfiles{i};
            fwrite(fp,t','float'); % write transpose to maintain matlab orientation     
        end
    end
    project.ntestfiles = nact_test;
    project.testfiles = incfiles;
    close(h); % close 
    fclose (fp);
end
y=1;
project.build=0;

%============================================
function name = getname(i)
%============================================
global project
name = lower([project.sdir '\\' project.files{i}]);
[a,b] = fileparts(name);
name = [a b '.dbc'];

