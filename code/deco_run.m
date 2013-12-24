%==================================
% deco_run
%==================================
function []=deco_run()
global project;

%more than 10 peaks may be unrealistic, inspect this later
% do not limit number of peaks per block (added 19-01-2008 JV)
%[a]=find(project.numpeaks>10); 
%numpeaks(a)=10; %force to maximum 10 peaks - should be a user input?

%==============================================================
%==== ALL FILES DECOVOLUTION 
%==============================================================
h=waitbar(0,'Deconvolute blocks');
recalc              = project.recalc;
numpeaks            = project.numpeaks;
excludeblocks       = project.ex_block;

%redetermine first and last block
start = mean(project.rt_start);
project.first =1; % prevents going below the lowest retention time
while ((start + (project.first+1)*(2*project.pw)*project.interval)/60<project.block_start)
      project.first = project.first + 1;
end

project.last = project.first; % prevents goinf over the last possible time
while ((start + (project.last+1)*(2*project.pw)*project.interval)/60<project.block_end && project.last<project.numblocks)
    project.last = project.last + 1;
end

for block = 1:project.numblocks
   project.deco{block} = [];
end

first = project.first;
last  = project.last;
changedpeaks = zeros(project.numblocks,1);

%tic; % start the analysis clock
% for block=project.first:project.last    % for each block       
    tic; % start the analysis clock
% try code begin, S. Krishnan
% Previous for-loop commented and inserted to check specific blocks
 for block=10:10    % for each block
    % try code end
   
    set(h,'Name',[' block:' num2str(block-first) ' of ' num2str(last-first+1)]);
    waitbar((block)/(last-first+1),h);
    clear decoresults;
    if size(project.deco,1)~=(last-first+1) || recalc(block)
        dodeco=1; % redo the deco calculations
    else
        if size(project.deco{block}.sopt,1)~=numpeaks(block) || recalc(block)
            dodeco=1;
        end
    end
    % recalc is set to 1 if a chrom is excluded from the calculations or
    % if the number of peaks is changed in inspect_decoresults
    if numpeaks(block)>0 && ~any(excludeblocks==block) && dodeco,
        changedpeaks(block)=1;
        % do actual deconvolution
        % replaced original code by deco_process_block routine 7-7-2006 JV
        % added fast mode in process 13-7-2006 JV
        deco_process_block(block,numpeaks(block),recalc(block),1);
        deco_predict(block); % do predictions
    else
        project.deco{block} = []; % empty block
    end
    time = toc;  % stop the clock and store the time 
    project.timepermass(block) = time;
end
%time = toc;  % stop the clock and store the time 
close(h); % close the progress waitbar

%project.time = time;
project.numpeaks = numpeaks;
project.changed  = changedpeaks;
