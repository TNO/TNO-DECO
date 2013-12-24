%---------------------------------------------
function y = deco_align(first, maxrange)
%---------------------------------------------
global project;
%---------------------------------------
% automatic alignment of spectra
% first : reference spectrum
% maxrange: max lag of the autocorrelation (in points default 50)
%--------------------------------------
%written by J.T.W.E. Vogels
% 29-01-2006
%-----------------------------------------
if project.nfiles==1
    y=0;
    return;
end

if first<1 || first > project.nfiles
    error('first file not correctly set');
end

ref=[]; % create a reference spectrum
name=project.files{first};
filename = [project.sdir name];
dbcname = [strtok(filename,'.') '.dbc'];
if (exist(dbcname,'file'))
   load('-mat',dbcname,'theMat_cdf','scantimes');  
   ref = sum(theMat_cdf,2); % get the reference tc spectrum
   align =0; %#ok<NASGU>
   save(dbcname,'-append','align');
end

if isempty(ref)
    message('no reference file available');
    return;
end

%-----------------------------------------------------------
%process the model files
%-----------------------------------------------------------
project.change_align = 1;
str = 'Align files';
h=waitbar(0,str);
for i=1:project.nfiles
    waitbar(i/(project.nfiles+project.ntestfiles),h);
    if (i ~= first)
        name=project.files{i};
        filename = [project.sdir name];
        dbcname = [strtok(filename,'.') '.dbc'];
        if exist(dbcname,'file') % only if data file existst
            % start from the original cdf files (clean of everything)
            load('-mat',dbcname,'theMat_cdf'); 
            spec=sum(theMat_cdf,2); % get the tic for the other spectra
            [p,q] = deco_corr(ref,spec,maxrange); % get cross correlation values
            [a,b] = max(p); % determine maximum correlation (a=value, b=pos
            corr = q(b); % the shift needed       
            % the time axis remains the same but the spectrum rotates corr units
            align=corr; %#ok<NASGU> % copy the alignment factor
            % now save the theMat structure for further processing
            save(dbcname,'-append','align'); % save the alignment factor
        end
    end
end

%-----------------------------------------------------------------
% process the test files
%-----------------------------------------------------------------
for i=1:project.ntestfiles
   waitbar((i+project.nfiles)/(project.nfiles+project.ntestfiles),h);
   name=project.testfiles{i};
   filename = [project.sdir name];
   dbcname = [strtok(filename,'.') '.dbc'];
   if exist(dbcname,'file') % only if data file existst
   % start from the original cdf files (clean of everything)
     load('-mat',dbcname,'theMat_cdf'); 
     spec=sum(theMat_cdf,2); % get the tic for the other spectra
     [p,q] = deco_corr(ref,spec,maxrange); % get cross correlation values
     [a,b] = max(p); % determine maximum correlation (a=value, b=pos
     corr = q(b); % the shift needed       
     % the time axis remains the same but the spectrum rotates corr units
     align=corr; %#ok<NASGU> % copy the alignment factor
     % now save the theMat structure for further processing
     save(dbcname,'-append','align'); % save the alignment factor
  end
end

close(h)
            
        
    

