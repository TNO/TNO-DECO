function varargout = deco_main(varargin)
% DECO_MAIN M-file for deco_main.fig
%      DECO_MAIN, by itself, creates a new DECO_MAIN or raises the existing
%      singleton*.
% Edit the above text to modify the response to help deco_main
% Last Modified by GUIDE v2.5 23-Jan-2009 11:49:04
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_main_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before deco_main is made visible.
function deco_main_OpeningFcn(hObject, eventdata, handles, varargin)
global project

% Choose default command line output for deco_main
handles.output = hObject;
%clear_project(hObject,eventdata,handles); % reset the content of the project
%reset project to default values  
project.basename  = 'tno1';
project.sdir      = ''; % files directory
project.files     = []; % filenames
project.nfiles    = 0;  % number of files
project.name      = 'tno1';
project.pw        = 20.0; % get linewidth
project.ex_mass   = [];   % excluded  masses
project.ex_block  = [];   % excluded blocks
project.build     = 1;     
project.resample  = 0;    % default do not resample
project.type      = 0;    % GCMS =0 LCMS = 1
project.reload    = 0;    % force a reload of the cdf files
project.lambda    = 1e6;    % baseline lambda
project.asym      = 0.001;  % baseline asymmetry
project.order     = 2;      % baseline order
project.basewin   = 15;     % baseline window
project.gsize     = 1200;   % slice size (groningen)
project.gmove     = 10;     % slice move
project.gthresh   = 200;    % threshhold
project.numpeaks  = [];    % the number of peaks estimated
project.first     = 0;  % first block in set
project.last      = 0;  % last block in set
project.start_mass= []; % start masses per file
project.end_mass  = []; % end masses per file
project.chrom_len = []; % length of traces in each file
project.minchromlen=0;  % minimal chrom length (determined in readdata)
project.minmass   = 0;   % set minmass
project.maxmass   = 0;   % set maxmass
project.theMasses = [];  % masses in the dataset
project.nmasses   = 0;   % number of masses in set
project.ntraces   = 0;   % number of traces in set
project.changed   = [];  % which peaks have changed
project.deco      = cell(0,0); % the deco resuls
project.recalc    = [];  % which blocks need reprocessing
project.version   = 1.8; % deco version
project.interval      = 0;   % global scanning interval
project.masscorr      = 1.0;  % mass correction
project.setbaseline   = 3; % default method = TNO
project.basemethod    = 1; % default baseline method 
project.lastblock     = 0; % zero means all blocks
project.warnings      = 0; % show no warnings
project.masscorr      = 1.0; % the mass correction factor
project.maxrange      = 30;  % maximum alignment range
project.diddeco       = 0;   % did we do the deco
project.smooth        = 0;   % do we smooth the spectra if yes then use this window
project.peakwindow    = 3;   % peaks are the same if position (in points) is x + or - this factor
project.maxiter       = 20;  % maximum number of iterations
project.rebuild       = 0;   % flag to rebuild the prepare structure
project.thresh        = 0.001; % threshold in selection of masses
project.rt_start      = [];  % start times     
project.iq            = [];  % iq values 
project.use_iq        = 1; %  default do not use iq
%project_iq_thresh     = 0.5; % default iq threshold
project.time          = 0.0; % time taken to do deconvolution
project.ntestfiles    =0;    % number of testfiles
project.testdir       =[];  % location of test file directory
project.testfiles     =[]; % list of test files
project.test_ex       = [];  % exclude/include test samples
project.level         = 2;     % noise level corrector in estimate peaks
project.sampleinc     = []; % include samples in calculations
project.change_align  = 0;
project.maxiter       = 50; % maxiteration for als
project.nist_dir      = ''; % default directory
project.full_box      =  1; % default output all dirty mass spectra 
project.clean_box     =  1; % default output clened set of mass  spectra
project.individual_box = 1; % default output areas for individual spectra
project.retention_box  = 1; % default output set of specific retention times 
project.align_ref      = 2;
project.classic        = 0; % classic output
project.sort_output    = 1; % sort the lines in the output on retention time
project.nist_dir       = 'c:\nistdemo\mssearch';
project.autoreduce     = 1;
project.maxn = 10;
project.usemasses  = []; % store the number of usemasses in each block % S. Krishnan
project.timepermass  = []; % store the time for computing each block % S. Krishnan


% reset common values
set(handles.align_button,'Enable','off');
set(handles.run_button,'Enable','off');
set(handles.estimate_button,'Enable','off');
set(handles.inspect_button,'Enable','off');
set(handles.baseline_button,'Enable','off');
set(handles.psavebutton,'Enable','off');
set(handles.prepare_button,'Enable','off');
set(handles.export_button,'Enable','off');   
set(handles.sampleselect,'Enable','off');
set(handles.globalview,'Enable','off');
set(handles.TIC_button,'Enable','off');
set(handles.exclude_tag,'Enable','off');
set(handles.basemethod,'Value',project.setbaseline);
set(handles.pw,'String',num2str(project.pw));
set(handles.lastblock,'Value',project.last);
set(handles.warnings,'Value',project.warnings==0);
set(handles.masscorr,'String',num2str(project.masscorr));
set(handles.project_name,'String',project.name);
set(handles.version,'String',['ver:' num2str(project.version)]);
set_choice(hObject,eventdata,handles);
set(handles.rebuild,'Value',0);
set(handles.arange,'String',num2str(project.maxrange));
set(handles.thresh_tag,'String',num2str(project.thresh));
set(handles.iq_flag,'Value',project.use_iq);
set(handles.timetag,'String','zero');
set(handles.ntest,'String',num2str(project.ntestfiles));
set(handles.show_sigma,'Value',0);
set(handles.max_iter,'String',num2str(project.maxiter));
set(handles.last_save,'String','');
set(handles.align_ref,'String',num2str(project.align_ref));
set(handles.AutoReduce,'Value',project.autoreduce);
set(handles.maxblockn,'String',num2str(project.maxn));
% set highlight on buttons to off
% Update handles structure
guidata(hObject, handles)

function clear_project(hObject,eventdata,handles) %#ok<INUSL,INUSL>
%reset project to default values   
project.sdir      = ''; % files directory
project.files     = []; % filenames
project.nfiles    = 0;  % number of files
project.name      = 'tno1';
project.pw        = 20.0; % get linewidth
project.ex_mass   = [];   % excluded  masses
project.ex_block  = [];   % excluded blocks
project.build     = 1;     
project.resample  = 0;    % default do not resample
project.type      = 0;    % GCMS =0 LCMS = 1
project.reload    = 0;    % force a reload of the cdf files
project.lambda    = 1e6;    % baseline lambda
project.asym      = 0.001;  % baseline asymmetry
project.order     = 2;      % baseline order
project.basewin   = 15;     % baseline window
project.gsize     = 1200;   % slice size (groningen)
project.gmove     = 10;     % slice move
project.gthresh   = 200;    % threshhold
project.numpeaks  = [];    % the number of peaks estimated
project.numblocks = 0;  % the number of blocks
project.first     = 0;  % first block in set
project.last      = 0;  % last block in set
project.start_mass= []; % start masses per file
project.end_mass  = []; % end masses per file
project.chrom_len = []; % length of traces in each file
project.minchromlen=0;  % minimal chrom length (determined in readdata)
project.minmass   = 0;   % set minmass
project.maxmass   = 0;   % set maxmass
project.theMasses = [];  % masses in the dataset
project.nmasses   = 0;   % number of masses in set
project.ntraces   = 0;   % number of traces in set
project.changed   = [];  % which peaks have changed
project.deco      = cell(0,0); % the deco resuls
project.recalc    = [];  % which blocks need reprocessing
project.version   = 1.8; % deco version
project.interval  = 0;   % global scanning interval
project.masscorr     = 1.0;  % mass correction
project.setbaseline  = 3; % default method = TNO
project.basemethod   = 1; % default baseline method 
project.lastblock    = 0; % zero means all blocks
project.warnings     = 0; % show no warnings
project.masscorr     = 1.0; % the mass correction factor
project.maxrange     = 50;  % maximum range
project.diddeco      = 0;   % did we do the deco
project.smooth       = 0;   % do we smooth the spectra if yes then use this window
project.peakwindow   = 3;   % peaks are the same if position (in points) is x + or - this factor
project.maxiter      = 20;  % maximum number of iterations
project.rebuild      = 0;   % flag to rebuild the prepare structure
project.thresh       = 1.00; % threshold in selection of masses
project.fast         = 0;   % use fast als (added 13-7-2006)
project.rt_start     = [];  % start times   
project.iq           = [];
project.use_iq       = 0; %  default do not use iq
%project_iq_thresh    = 0.5; % default iq threshold
project.time         = 0;
project.ntestfiles   = 0;
project.testfiles    = [];
project.testdir      = [];
project.test_ex      = [];
project.level        = 2;  % default level
project.maxiter      = 50; % default number of iterations
project.nist_dir      = '';
project.full_box      =  1; % default output all dirty mass spectra 
project.clean_box     =  1; % default output clened set of mass  spectra
project.individual_box = 1; % default output areas for individual spectra
project.retention_box  = 1; % default output set of specific retention times 
project.align_ref      = 1; % alignment reference spectrum
project.classic        = 0;
project.sort_output    = 1;
project.autoreduce     = 1;
project.maxn           = 10;


% reset common values
set(handles.align_button,'Enable','off');
set(handles.run_button,'Enable','off');
set(handles.estimate_button,'Enable','off');
set(handles.inspect_button,'Enable','off');
set(handles.baseline_button,'Enable','off');
set(handles.psavebutton,'Enable','off');
set(handles.prepare_button,'Enable','off');
set(handles.export_button,'Enable','off');   
set(handles.sampleselect,'Enable','off');
set(handles.globalview,'Enable','off');
set(handles.TIC_button,'Enable','off');
set(handles.exclude_tag,'Enable','off');

% --- Outputs from this function are returned to the command line.---
function varargout = deco_main_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL,INUSL>
%--------------------------------------------------------------------
varargout{1} = handles.output;

% --- Executes on button press in select button1.-------------------
function select_Callback(hObject, eventdata, handles) %#ok<INUSL>
%------------------------------------------------------------------
global project
[files,dir]=uigetfiles('*.cdf','Read NetCDF files');
if size(files,2)==0 % cancel pressed
    return;
end
files = sort(files); % sort files obtained from uigetfiles

project.reload = get(handles.reload,'Value');
project.nfiles = 0;
% whole project should be reset when a new set is read
set(handles.run_button,'Enable','on');
set(handles.baseline_button,'Enable','on');
handles.files = files;
handles.dir   = dir;
project.sdir  = dir;
project.files = files;
project.build = 1;   % reset the building to do a (re)building
project.masscorr = str2double(get(handles.masscorr,'String')); % get entry from main window
project.nfiles = size(files,2); % determine the number of files in set
set(handles.numfiles,'String',num2str(project.nfiles)); % visualize the number of files

%------------------------------------------------------------------
deco_readset(project.files,project.sdir,project.nfiles); % read and convert the dataset
project.ntraces= 0;
project.interval=0;
project.maxtime = 0;
project.mintime = 0;
for i=1:project.nfiles
    name = project.files{i};
    filename = [project.sdir name]; % construct full name = directory + filename
    load('-mat',[strtok(filename,'.') '.dbc'],'resample','ntraces','scantimes');
    
    % determine the global interval between scans
    if (i==1 || resample.interval<project.interval) 
            project.interval = resample.interval;
    end 
    %ntraces
   
    %determine the number of traces
    if (i==1 || ntraces<project.ntraces)
            project.ntraces = ntraces;
    end   
    m = min(scantimes);
    if (i==1 || m > project.mintime)
        project.mintime = m;
    end
    
    m = max(scantimes);
    if (i==1 || m < project.maxtime)
        project.maxtime = m;
    end
    
          
end
%-----------------
pw  = str2double(get(handles.pw,'String')); % get linewidth
% deze waarde is niet goed meer na resampling !
nb = floor(project.ntraces / (2*pw)) -2
endtime = (project.mintime + nb*project.pw*2*project.interval)/60

% beter om de start en eind tijd op te geven ipv van blokken
set(handles.start_time,'String',num2str(project.mintime/60));
set(handles.end_time,'String',num2str(endtime));

set(handles.interval,'String',num2str(project.interval)); % visualize interval time
guidata(hObject,handles);
project.name = get(handles.project_name,'String');
deco_save_project; % save the project 
project.last = nb;
project.first = 1;

set(handles.project_name,'String',project.name);
set(handles.reload,'Value',project.reload);
set(handles.align_button,'Enable','on');
set(handles.run_button,'Enable','on');
set(handles.estimate_button,'Enable','off');
set(handles.inspect_button,'Enable','off');
set(handles.baseline_button,'Enable','on');
set(handles.psavebutton,'Enable','on');
set(handles.prepare_button,'Enable','on');
set(handles.export_button,'Enable','off');
set(handles.sampleselect,'Enable','off');
set(handles.globalview,'Enable','on');
set(handles.TIC_button,'Enable','on');
set(handles.exclude_tag,'Enable','on');

set(handles.firstblock,'String',num2str(project.first));
set(handles.lastblock,'String',num2str(project.last));
set(handles.ntest,'String',num2str(project.ntestfiles));
set(handles.sampleselect,'Enable','on');
set(handles.max_iter,'String',num2str(project.maxiter));


function select_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on button press in pushbutton2.---------------
function baseline_button_Callback(hObject, eventdata, handles)
%------------------------------------------------------------
global project;
% communiceer the base method
 project.basemethod =get(handles.basemethod','Value');
 project.gsize      = str2double(get(handles.slice_size,'String'));
 project.gmove      = str2double(get(handles.slice_move,'String'));
 project.gthresh    = str2double(get(handles.slice_thresh,'String'));
 project.basewin    = str2double( get(handles.basewin,'String'));
 project.lambda     = str2double(get(handles.lambda,'String'));
 project.asym       = str2double(get(handles.asym,'String'));
 project.order      = str2double(get(handles.order_edit,'String'));
 project.thresh     = str2double(get(handles.thresh_tag,'String'));
[a]=deco_base_optimize;
if strcmpi(a,'OK')
    disp('copy');
    set(handles.basemethod,'Value',project.basemethod);
    set_choice(hObject,eventdata,handles) 
    if project.basemethod==1
        set(handles.lambda,'String',num2str(project.lambda));
        set(handles.order_edit,'String',num2str(project.order));
        set(handles.asym,'String',num2str(project.asym));
    elseif project.basemethod==2
        set(handles.slice_size,'String',num2str(project.gsize));
        set(handles.slice_move,'String',num2str(project.gmove));
        set(handles.slice_thresh,'String',num2str(project.gthresh));
    elseif project.basemethod==3
        set(handles.basewin,'String',num2str(project.basewin)); 
    end    
end

% --- Executes on button press in est_peaks.
function est_peaks_Callback(hObject, eventdata, handles)

function project_name_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function project_name_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in rundeco.------------------------
function rundeco_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
%----------------------------------------------------------------

%-----------------------------------------------------------------
function pw_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
%-----------------------------------------------------------------
%pwb  = str2double(get(handles.pw,'String')); % get linewidth
% now we have to recalculate the first and last blocks


% --- Executes during object creation, after setting all properties.
function pw_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------
[files,sdir]=uigetfile({'*.prj','Project files (*.prj)'},'Open project Files');
if (files==0) % no files selected
    return;
end
load('-mat',[sdir files], 'project');

project.sort_output = 1;
project.autoreduce  = 1;


set(handles.numfiles,'String',num2str(project.nfiles));
set(handles.run_button,'Enable','on');
set(handles.pw,'String',num2str(project.pw)); % get the linewidth
set(handles.project_name,'String',project.name); 
set(handles.reload,'Value',project.reload);       % reload samples
set(handles.baseline_button,'Enable','on');
set(handles.interval,'String',num2str(project.interval'));
set(handles.lastblock,'String',num2str(project.last));
set(handles.firstblock,'String',num2str(project.first));
set(handles.thresh_tag,'String',num2str(project.thresh));
set(handles.numtraces,'String',num2str(project.ntraces));
set(handles.nummasses,'String',num2str(project.nmasses));
set(handles.ntest,'String',num2str(project.ntestfiles));
set(handles.max_iter,'String',num2str(project.maxiter));

t = project.time;
hours  = floor(t/3600);  t = t - hours*3600;
minute = floor(t/60);    t = t - minute * 60;
secs   = floor(t);       
set(handles.timetag,'String',[num2str(hours) ':' num2str(minute) ':' num2str(secs)]);

% copy project baseline_button parameters
set(handles.basemethod,'Value',project.basemethod);
set_choice(hObject,eventdata,handles) 
if project.basemethod==1
     set(handles.lambda,'String',num2str(project.lambda));
     set(handles.order_edit,'String',num2str(project.order));
     set(handles.asym,'String',num2str(project.asym));
elseif project.basemethod==2
    set(handles.slice_size,'String',num2str(project.gsize));
    set(handles.slice_move,'String',num2str(project.gmove));
    set(handles.slice_thresh,'String',num2str(project.gthresh));
elseif project.basemethod==3
    set(handles.basewin,'String',num2str(project.basewin));   
end
set(handles.arange,'String',num2str(project.maxrange));
set(handles.rebuild,'Value',0);

set(handles.baseline_button,'Enable','off');
set(handles.align_button,'Enable','off');
set(handles.psavebutton,'Enable','off');
set(handles.prepare_button,'Enable','off');
set(handles.run_button,'Enable','off');

set(handles.estimate_button,'Enable','off');
set(handles.inspect_button,'Enable','off');
set(handles.export_button,'Enable','off');
set(handles.last_save,'String',['last:' datestr(now,'HH:MM:SS')]);

set(handles.start_time,'String',num2str(project.block_start));
set(handles.end_time,'String',num2str(project.block_end));
set(handles.align_ref,'String',num2str(project.align_ref));

if project.nfiles>=1 % different files have been read
    set(handles.sampleselect,'Enable','on');
    set(handles.globalview,'Enable','on');
    set(handles.TIC_button,'Enable','on');
    set(handles.exclude_tag,'Enable','on');    
    set(handles.align_button,'Enable','on');
    set(handles.baseline_button,'Enable','on');
    set(handles.psavebutton,'Enable','on');
    set(handles.run_button,'Enable','on');
    set(handles.prepare_button,'Enable','on');
    if project.build==0 || project.diddeco==1
       set(handles.estimate_button,'Enable','on');
    end
    if project.diddeco==1
          set(handles.inspect_button,'Enable','on');
          set(handles.export_button,'Enable','on');
    end
end

guidata(hObject,handles);

% --- Executes on button press in psavebutton.-------------------
function psavebutton_Callback(hObject, eventdata, handles) %#ok<INUSL,INUSL,DEFNU>
%---------------------------------------------------------------
global project
project.name = get(handles.project_name,'String'); % get project name 
deco_save_project; % save the project
set(handles.project_name,'String',project.name);
set(handles.last_save,'String',['last:' datestr(now,'HH:MM:SS')]);

% --- Executes on button press in resample.
function resample_Callback(hObject, eventdata, handles) %#ok<INUSD,INUSD,DEFNU>
%--------------------------------------------------------------

function edit3_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
%----------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject,'String','Hello');
end

%--------------------------------------------------------------
function numfiles_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------

% --- Executes on button press in method.------------------------
function method_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------

%---------------------------------------------
function save_project(handles)
%---------------------------------------------
deco_save_project; % save the project
set(handles.project_name,'String',project.name);

% --- Executes on button press in reload.
function reload_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------

%---------------------------------------------------------
function lambda_Callback(hObject, eventdata, handles)
%----------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)
%----------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------
function asym_Callback(hObject, eventdata, handles)
%--------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function asym_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%---------------------------------------------------------------
function order_edit_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function order_edit_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in prepare_button.
function prepare_button_Callback(hObject, eventdata, handles)
 %--------------------------------------------------------------- 
global project;
% force a rebuild of the preparation
project.build=1;
do_prepare(hObject,eventdata,handles);

%-----------------------------------------------------
function do_prepare(hObject, eventdata, handles) %#ok<INUSL,INUSL>
%-----------------------------------------------------
global project
% get the set values
%detect if the settings have been changed
%if so reproccess the data    

set(handles.inspect_button,'Enable','off');
set(handles.export_button,'Enable','off');
project.use_iq = get(handles.iq_flag,'Value');
 
str2double(get(handles.lastblock,'String'));
str2double(get(handles.firstblock,'String'));


pw  = str2double(get(handles.pw,'String')); % get linewidth
if (pw ~= project.pw); project.pw=pw;   project.build=1; end

lastblock = str2double(get(handles.lastblock,'String')); % get last block used in testing
if project.last ~= lastblock || lastblock > size(project.numpeaks,1)
    project.last=lastblock;
    if project.last > size(project.numpeaks,1);  project.build=1; end
end

firstblock = str2double(get(handles.firstblock,'String')); % get first block to process
if project.first ~= firstblock
    project.first =firstblock;    
    project.build=1;
end
    
% determine the current set scalingmethod
basemethod = get(handles.basemethod,'Value');
if basemethod ~= project.basemethod; 
    project.basemethod=basemethod; 
    project.build=1; 
end

% get the current set resample interval
interval = str2double(get(handles.interval,'String'));
if project.interval ~= interval; 
    project.interval=interval; 
    project.build=1; 
end;

if project.change_align
    project.build=1;
    project.change_align=0;
end

% was anything changed in the settings for the baseline correction
if project.basemethod==1 % Eilers
    asym = str2double(get(handles.asym,'String'));
    if project.asym ~= asym; project.asym=asym; project.build=1; end;
    lambda=str2double(get(handles.lambda,'String'));
    if project.lambda ~= lambda; project.lambda = lambda; project.build=1; end;
    order = str2double(get(handles.order_edit,'String'));
    if project.order ~= order; project.order=1; project.build=1; end;
elseif project.basemethod==2 % groningen
    gsize = str2double(get(handles.slice_size,'String'));
    if gsize ~=project.gsize; project.gsize = gsize; project.build=1; end;
    gmove = str2double(get(handles.slice_move,'String'));
    if gmove ~= project.gmove; project.gmove =gmove; project.build=1; end;
    gthresh=str2double(get(handles.slice_thresh,'String'));
    if gthresh ~=project.gthresh; project.gthresh=gthresh; project.build=1; end;
else % default
    basewin = str2double(get(handles.basewin,'String'));
    if basewin ~= project.basewin; project.basewin=basewin; project.build=1; end;
end

% reprocess the data
if project.build==1 || get(handles.rebuild,'Value')==1  
    project.build=1;
    project.rebuild   = get(handles.rebuild,'Value');    
    project.pw        = str2double(get(handles.pw,'String')); % get linewidth
    project.lastblock = str2double(get(handles.lastblock,'String')); % get last block used in testing
    project.block_end  = str2double(get(handles.end_time,'String'));
    project.block_start   = str2double(get(handles.start_time,'String'));
    deco_prepare;      % do data preparation step (baseline, resample, peaks estimate)
    deco_save_project; % save the project
    set(handles.rebuild,'Value',0);
    project.rebuild=0;
end

set(handles.estimate_button,'Enable','on'); % enable estimate tool
set(handles.nummasses,'String',num2str(project.nmasses));
set(handles.numtraces,'String',num2str(project.ntraces));
set(handles.firstblock,'String',num2str(project.first));
set(handles.lastblock,'String',num2str(project.last));

project.diddeco=0;
deco_save_project; % save the project

% --- Executes on button press in run_button.--------------------
function run_button_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------
global project
project.pw        = str2double(get(handles.pw,'String')); % get linewidth
project.lastblock = str2double(get(handles.lastblock,'String')); % get last block used in testing
project.last      = project.lastblock;
project.thresh    = str2double(get(handles.thresh_tag,'String'));
project.maxiter = str2double(get(handles.max_iter,'String'));
project.block_start = str2double(get(handles.start_time,'String'));
project.block_end = str2double(get(handles.end_time,'String'));

do_prepare(hObject, eventdata, handles); % run the preparation again (only necessary elements)

if project.warnings==0 % show no warnings
    warning off MATLAB:nearlySingularMatrix;
    warning off MATLAB:rankDeficientMatrix
end
%---------------------------------------------
% do the full deconvolution of all blocks
%--------------------------------------------
deco_run;

t = project.time;
hours  = floor(t/3600);  t = t - hours*3600;
minute = floor(t/60);    t = t - minute * 60;
secs   = floor(t);       
set(handles.timetag,'String',[num2str(hours) ':' num2str(minute) ':' num2str(secs)]);

project.diddeco=1; % we have done a deco deconvolution
deco_save_project; % save the project

% reset the warning settings
if project.warnings==0
    warning on MATLAB:nearlySingularMatrix;
    warning on MATLAB:rankDeficientMatrix
end

%-----------------------------------------------
% enable all the tools in the display
%-----------------------------------------------
set(handles.estimate_button,'Enable','on');
set(handles.inspect_button,'Enable','on');
set(handles.export_button,'Enable','on');

for i=project.first:project.last
    deco_block_evaluate(i); % do prepare for export
end

deco_evaluate;
deco_save_project;

deco_inspect; % switch to visualisation mode

%--------------------------------------------------------------
function frequency_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function frequency_CreateFcn(hObject, eventdata, handles)
%---------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in resample.-------------------------
function radiobutton6_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------

%-------------------------------------------------------------
function interval_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function interval_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in spikes. -----------------------
function spikes_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------

%---------------------------------------------------------------
function masscorr_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function masscorr_CreateFcn(hObject, eventdata, handles)
%---------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in baselineset.---------------------
function baseset_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------------

%---------------------------------------------------------------
function lastblock_Callback(hObject, eventdata, handles)
global project;
lastblock = str2double(get(hObject,'String'));
firstblock = str2double(get(handles.firstblock,'String'));
nb = floor(project.ntraces/(2*project.pw))-2;
if (lastblock>nb || lastblock < firstblock)
    lastblock=nb; 
    set(handles.lastblock,'String',num2str(lastblock));
end;
endtime = (project.mintime + lastblock*project.pw*2 *project.interval)/60;
set(handles.end_time,'String',num2str(endtime));

%--------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function lastblock_CreateFcn(hObject, eventdata, handles)
%----------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in inspect_button.-------------------
function inspect_button_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------
deco_inspect

% --- Executes on selection change in basemethod.-----------------
function basemethod_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------
set_choice(hObject,eventdata,handles);

%----------------------------------------------   
function set_choice(hObject,eventdata,handles)
%--------------------------------------------
n = get(handles.basemethod,'Value');
set(handles.lambda,'Visible','Off');
set(handles.order_edit,'Visible','Off');
set(handles.asym,'Visible','Off');
set(handles.slice_size,'Visible','Off');
set(handles.slice_move,'Visible','Off');
set(handles.slice_thresh,'Visible','Off');
set(handles.basewin,'Visible','Off');

if n==1
    set(handles.text_1,'String','lambda');
    set(handles.text_a,'String','asym');
    set(handles.text_3,'String','order');
    set(handles.lambda,'Visible','On');
    set(handles.order_edit,'Visible','On');
    set(handles.asym,'Visible','On');
    
elseif n==2
    set(handles.text_1,'String','size');
    set(handles.text_a,'String','move');
    set(handles.text_3,'String','thresh');
    set(handles.slice_size,'Visible','On');
    set(handles.slice_move,'Visible','On');
    set(handles.slice_thresh,'Visible','On');
elseif n==3
    set(handles.text_1,'String','window');
    set(handles.text_a,'String','');
    set(handles.text_3,'String','');
    set(handles.basewin,'Visible','On'); 
end

% --- Executes during object creation, after setting all properties.
function basemethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------------
function basewin_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function basewin_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-------------------------------------------------------------
function slice_size_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function slice_size_CreateFcn(hObject, eventdata, handles)
%---------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------------
function slice_move_Callback(hObject, eventdata, handles)
%------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function slice_move_CreateFcn(hObject, eventdata, handles)
%---------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------
function slice_thresh_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function slice_thresh_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in warnings.------------------
function warnings_Callback(hObject, eventdata, handles)
%------------------------------------------------------------

% --- Executes on button press in rebuild.------------------
function rebuild_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------

% --- Executes on button press in align_button.----------------------
function align_button_Callback(hObject, eventdata, handles) 
%-------------------------------------------------------------
global project;
first = str2double(get(handles.align_ref,'String'));
maxrange = str2double(get(handles.arange,'String'));
if (first < 0) 
    first=1;
end

if (first > project.nfiles)
    first = project.nfiles;
end

project.align_ref = first;
set(handles.align_ref,'String',num2str(first));

if  project.nfiles>1 
     deco_align(first,maxrange); % project.align_range
end

%---------------------------------------------------------------
function arange_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function arange_CreateFcn(hObject, eventdata, handles)
%----------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in estimate_button.---------------
function estimate_button_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------
deco_estimate_tool;

%------------------------------------------------------------------
function nummasses_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function nummasses_CreateFcn(hObject, eventdata, handles)
%-----------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%---------------------------------------------------------
function numtraces_Callback(hObject, eventdata, handles)
%--------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function numtraces_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-----------------------------------------------------------
function firstblock_Callback(hObject, eventdata, handles)
global project;
firstblock = str2double(get(hObject,'String'));
lastblock = str2double(get(handles.lastblock,'String'));
if (firstblock<1 || firstblock>=lastblock)
    firstblock=1;
    set(hObject,'String',num2str(firstblock));
end
set(handles.start_time,'String',num2str((project.mintime + project.pw*2*project.interval*(firstblock-1))/60));

%------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function firstblock_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-----------------------------------------------------------------
function smooth_window_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function smooth_window_CreateFcn(hObject, eventdata, handles)
%---------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in smooth.
function smooth_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------

% --- Executes on button press in export_button.------------------
function export_button_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------------
global project;
ht = waitbar(0,'creating Exportdata');
deco_export;  

if (project.classic)
    temp_deco_clean;  % first clean the data
    waitbar(0.3,ht,'cleaning data');
end
waitbar(0.6,ht,'export v2.0');
%deco_exp
waitbar(1.0,ht,'ready export');

pause(1);
close(ht);

%disp('export done')

function start_time_Callback(hObject, eventdata, handles)
global project;
s_time = str2double(get(hObject,'String'))*60; % start time
if (s_time<project.mintime) 
    s_time = project.mintime;
    set(hObject,'String',num2str(s_time));
end
first_block = floor((s_time - project.mintime)/(project.pw*2*project.interval))+1;
set(handles.firstblock,'String',num2str(first_block));

    
% --- Executes during object creation, after setting all properties.
function start_time_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function end_time_Callback(hObject, eventdata, handles)
global project;
s_time = str2double(get(hObject,'String')); % start time
nb = floor(project.ntraces/(2*project.pw))-2;
last_time = (project.mintime + nb*project.pw*2*project.interval)/60;
if (s_time>last_time) 
    s_time = last_time;
    set(hObject,'String',num2str(s_time));
end
blks= floor((s_time*60 - project.mintime) / (project.pw*2*project.interval));
set(handles.lastblock,'String',num2str(blks));
    


% --- Executes during object creation, after setting all properties.
function end_time_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thresh_tag_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function thresh_tag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in TIC_button.-----------------
function TIC_button_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------
global project;
if project.nfiles>0 
    deco_tic;
end;
    
% --- Executes on button press in fast_method.----------------
function fast_method_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------

% --- Executes on button press in exclude_tag.------------------
function exclude_tag_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------
global project;
if (project.nfiles>0) 
    deco_exclude; 
end;

%---------------------------------------------------------------
function new_button_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------
global project;
[a,b]=deco_projectname('filename',project.name); % collect answer and new projectname
if a(1)=='O'
    project.name = b;
    set(handles.project_name,'String',project.name); 
    clear_project(hObject,eventdata,handles);
end

  
function iq_flag_Callback(hObject, eventdata, handles)

function timetag_Callback(hObject, eventdata, handles)

function timetag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function select_test_Callback(hObject, eventdata, handles)
global project
[files,dir]=uigetfiles('*.cdf','Read NetCDF files');
if size(files,2)==0 % cancel pressed
    return;
end

files = sort(files); % sort files obtained from uigetfiles
project.reload = get(handles.reload,'Value');
project.testdir  = dir;
project.ntestfiles = size(files,2);
project.testfiles = files;
deco_readset(project.testfiles,dir,project.ntestfiles); % read and convert the dataset
set(handles.ntest,'String', num2str(project.ntestfiles));
set(handles.numfiles,'String', num2str(project.nfiles));

function sampleselect_Callback(hObject, eventdata, handles)
global project;
if project.nfiles + project.ntestfiles > 0
    deco_select_dialog();
    set(handles.ntest,'String', num2str(project.ntestfiles));
    set(handles.numfiles,'String', num2str(project.nfiles));
else
    disp('no data present')
end

function edit26_Callback(hObject, eventdata, handles)

function edit26_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function show_sigma_Callback(hObject, eventdata, handles)

function ntest_Callback(hObject, eventdata, handles)

function ntest_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function globalview_Callback(hObject, eventdata, handles)
deco_overview;

function max_iter_Callback(hObject, eventdata, handles)

function max_iter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in global_button.
function global_button_Callback(hObject, eventdata, handles)
deco_global_params;

function align_ref_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function align_ref_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in AutoReduce.
function AutoReduce_Callback(hObject, eventdata, handles)






function maxblockn_Callback(hObject, eventdata, handles)
% hObject    handle to maxblockn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxblockn as text
%        str2double(get(hObject,'String')) returns contents of maxblockn as a double


% --- Executes during object creation, after setting all properties.
function maxblockn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxblockn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


