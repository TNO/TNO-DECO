function [] = deco_readset(ff,sdir,nfiles)
global project;

%----------------------------------------------
% batch version of reading files 
% and converting it to a dbc file
%----------------------------------------------
version = project.version;
resample.do       = 0;  % do resample this file
resample.start    = 0.0; % the resample start time
resample.end      = 0.0; % the resample end time 
resample.interval = 0.0; % the resample interval
masscorr          = project.masscorr; % mass correction factor
project.chrom_len = [];

% read a set of spectra
str='Read files';
ht = waitbar(0,str);

for i=1:nfiles % process all files in the filesbox    
    name = ff{i};
    set(ht,'Name',[name ' file:' num2str(i) ' of ' num2str(project.nfiles)]);
    waitbar(i/project.nfiles,ht);
    filename = [sdir name]; % construct full name = directory + filename
    dbcname  = [strtok(filename,'.') '.dbc'];
    if exist(dbcname,'file') % was this file with an old deco version if so then reload
        load('-mat',[strtok(filename,'.') '.dbc'],'version','masscorr');
    end
    % only reload if necessary
    if (exist(dbcname,'file')==0 || project.reload==1 || version ~= project.version || masscorr ~= project.masscorr) % if force load or dbc does not exist then reload   
        [f1,f2,f3] = fileparts(filename);
        align =0;
        spikes = 0;    
        if strcmpi(f3,'.ddf')
            disp('ddfimport')
            [theMat_cdf,a,scantimes] = ddfread(filename);
            %startmass = min(a);      
        elseif strcmpi(f3,'.cdf')     
            [theMat_cdf,startmass,scantimes_cdf] = deco_cdf2mat(strtok(filename,'.'),1.0); % read cdf file
        else
            disp('Could not read data');
            exit(0);
        end
        ntraces = size(theMat_cdf,1);
        scantimes = scantimes_cdf;
        version = project.version;        % reset the version name
     	resample.start = min(scantimes_cdf);
        resample.end   = max(scantimes_cdf);
            
        % calculate the minimum scantime per interval
        interval = min(diff(scantimes_cdf));
       
        resample.interval = interval;  
        %resample
        % here the data is restored in the dbc file                
        theMat = theMat_cdf; %#ok<NASGU>
        options.asym   = 0;
        options.lambda = 0;
        options.order  = 0;
        options.gsize  = 0;
        options.gmove  = 0;
        options.gthresh=0;
        options.basewin=0; 
        save([strtok(filename,'.') '.dbc'],'theMat','theMat_cdf','startmass','scantimes','scantimes_cdf','resample','options','masscorr','align','spikes','masscorr','ntraces','version');                
        clear 'theMat' 'startmass' 'scantimes_cdf' 'scantimes' 'theMat_cdf';           
        project.build=0; % build again     
    end   
end
close(ht)

project.reload=0; % do not reload again