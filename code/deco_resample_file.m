function [] = deco_resample_file(filename,interval,force_resample)
global project;

%--------------------------------------
% deco_resample(filename,interval)
% function to resample a file
% written by J. Vogels
% parameters
%     filename of the dbc file
%     set interval for the whole set
%--------------------------------------
load('-mat',filename,'scantimes','resample'); % load from previous run

minscan = min(diff(scantimes)); %#ok<NODEF>
maxscan = max(diff(scantimes));
range = interval/100; % 1% deviation allowed
%is it necessary to rescan the files ?
% if maximum scanrange or minimum scanrange
% deviates more than 1% from set interval then yes
old = resample.do;

if ( abs(minscan-interval)>range || abs(maxscan-interval)>range)
       resample.do=1;
else
    resample.do=0;
end

%===========================================================
% resample the files
%===========================================================
if (resample.do)
    load('-mat',filename,'theMat'); % load from previous run   
    newtimes = resample.start:interval:resample.end;
    resample.interval = project.interval;
    remat = zeros(length(newtimes),size(theMat,2),'single'); %#ok<NODEF>
    [a,b] = fileparts(filename);
    str = ['resample ' b];
    vt1 = waitbar(0,str);
    for mass=1:size(theMat,2)
        waitbar(mass/size(theMat,2),vt1);
       remat(:,mass) =spline(scantimes,theMat(:,mass),newtimes);
    end
    close(vt1);
    theMat =remat; %#ok<NASGU>
    scantimes =newtimes; %#ok<NASGU>
    save(filename,'-append','theMat','scantimes','resample');
    clear remat;
else
    %disp('no resampling necessary')
    if old ~= resample.do
        save(filename,'-append', 'resample');
    end
end

