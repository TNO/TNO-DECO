function varargout = deco_estimate_tool(varargin)
global start_block;
global disable_save;
% DECO_ESTIMATE_TOOL M-file for deco_estimate_tool.fig
%      DECO_ESTIMATE_TOOL, by itself, creates a new DECO_ESTIMATE_TOOL or raises the existing
%      singleton*.
%
%      H = DECO_ESTIMATE_TOOL returns the handle to a new DECO_ESTIMATE_TOOL or the handle to
%      the existing singleton*.
%
%      DECO_ESTIMATE_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECO_ESTIMATE_TOOL.M with the given input arguments.
%
%      DECO_ESTIMATE_TOOL('Property','Value',...) creates a new DECO_ESTIMATE_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deco_estimate_tool_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deco_estimate_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help deco_estimate_tool
% Last Modified by GUIDE v2.5 13-Aug-2008 10:20:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_estimate_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_estimate_tool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
start_block=1;
disable_save =0;

if nargin > 1
    disable_save= 1;   
end

if nargin >0 
    start_block=varargin{1};
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before deco_estimate_tool is made visible.
function deco_estimate_tool_OpeningFcn(hObject, eventdata, handles, varargin)
global start_block;
global values;
global disable_save;
% Choose default command line output for deco_estimate_tool
handles.output = hObject;
if start_block > 0
    handles.block=start_block;
else
    handles.block = 1;
end
values =[];
% Update handles structure
guidata(hObject, handles);
prepare_block(handles.block); % read the first data block
set(handles.screeplot,'Value',0);
set(handles.imbedded,'Value',1);
set(handles.indicator,'Value',1);
if (disable_save) 
    set(handles.save_button,'Enable','off');
end
       
plot_data(hObject,eventdata,handles);

% --- Outputs from this function are returned to the command line.
function varargout = deco_estimate_tool_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL,INUSL>
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in block_plus.
function block_plus_Callback(hObject, eventdata, handles) %#ok<DEFNU>
global project;
if (handles.block<project.last)
    handles.block = handles.block+1;
    guidata(hObject,handles);
    prepare_block(handles.block); % read the first data block
    plot_data(hObject,eventdata,handles);
end

% --- Executes on button press in block_min.
function block_min_Callback(hObject, eventdata, handles) %#ok<DEFNU>
if (handles.block>1)
    handles.block = handles.block-1;
    guidata(hObject,handles);
    prepare_block(handles.block); % read the first data block
    plot_data(hObject,eventdata,handles);
end

% --- Executes on button press in zoom.
function zoom_Callback(hObject, eventdata, handles) %#ok<INUSD,INUSD,DEFNU>
zoom;

% --- Executes on button press in full.
function full_Callback(hObject, eventdata, handles) %#ok<DEFNU,INUSD,INUSD>
axis auto;

% --- Executes on button press in show_log.
function show_log_Callback(hObject, eventdata, handles) %#ok<DEFNU>
plot_data(hObject,eventdata,handles);

function plot_data(hObject,eventdata,handles) %#ok<INUSL,INUSL>
global project;
global values;
global ssim;
global nfac;
global datablock;
global ie;
global ind;
global s2n;

dd = values;
ds = ssim;
legendstr = [];
n=0;
cla;
hold on;
plot(dd,'-b'); 
n=n+1; legendstr{n} = 'data'; 
plot(ds,'ob');
n=n+1; legendstr{n} = 'noise'; 
ie = ie / max(ie) * max(values);
ind = ind / max(ind) * max(values);
if get(handles.imbedded,'Value') ==1 
    plot(ie,'g');
    n=n+1; legendstr{n} = 'ie'; 
end
if get(handles.indicator,'Value') ==1
    plot(ind,'m');
    n=n+1; legendstr{n} = 'ind'; 
end
idx = find(dd>=1.0);
if get(handles.screeplot,'Value') == 1
    plog = log(dd(idx));
    mxlog = max(plog);
    plot(idx,plog/mxlog*max(dd),'r')
    n=n+1; legendstr{n} = 'scree';
end
yrange=ylim;
line([nfac nfac],[yrange(1) yrange(2)],'Color','b');
line([s2n s2n],[yrange(1) yrange(2)],'Color','k');
hold off;
k = [];
for j=1:project.nfiles
     b = datablock((j-1)*4*project.pw+1:j*4*project.pw,:);
     b = sum(b,2); % get tic
     np = deco_peakpick(b,-1,0,8,true,8);
     k = [k np];
end
k =unique(k); 
np = sum(diff(k)>1)+1;
rt_mean = mean(project.rt_start)/60;
rt_st = ((handles.block-1)*2*project.pw*project.interval/60 + rt_mean); % time in seconds
rt_ed = rt_st + (4*project.pw*project.interval)/60;
title(['estimate peaks block ' num2str(handles.block) ' nsig:' num2str(nfac) ' np:' num2str(np) ' [' num2str(rt_st) '-' num2str(rt_ed) ']' ]);
legend(legendstr);

% estimate the number of compounds in the block
function prepare_block(blocknum)
global project;
global datablock;
global blocktic;
global values;
global ssim;
global nfac;
global ie;
global ind;
global s2n;

datablock=deco_readblock([project.name '.bin'],blocknum,project.pw,project.nfiles,project.ntraces,project.nmasses);
blocktic=deco_readblocktic([project.name '.bin'],blocknum,project.pw,project.nfiles,project.ntraces,project.nmasses);

[nfac,np,values,ssim,ie,ind,s2n]=deco_estimate_block(blocknum,2,20);
%s2n=0;
nfac = project.numpeaks(blocknum); % this was the value calculated

%numsig,npeak,values,s_sim,ie,g_sim
% --- Executes on button press in plus_button.
function plus_button_Callback(hObject, eventdata, handles)
global project;
global nfac;
nfac = nfac +1;
plot_data(hObject,eventdata,handles);
    
% --- Executes on button press in min_button.
function min_button_Callback(hObject, eventdata, handles)
global project;
global nfac;
nfac = nfac-1;
plot_data(hObject,eventdata,handles);

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
global datablock;
global project;
global nfac;
project.numpeaks(handles.block)=nfac;

fname = ['block' num2str(handles.block)];
tics=zeros(project.nfiles,4*project.pw);
for i=1:project.nfiles
    f1 = sum(datablock((i-1)*4*project.pw+1:i*4*project.pw,:),2);
    tics(i,:) = f1;
end

save(fname,'datablock','tics');

% --- Executes on button press in screeplot.
function screeplot_Callback(hObject, eventdata, handles)
plot_data(hObject,eventdata,handles);

% --- Executes on button press in imbedded.
function imbedded_Callback(hObject, eventdata, handles)
plot_data(hObject,eventdata,handles);

% --- Executes on button press in indicator.
function indicator_Callback(hObject, eventdata, handles)
plot_data(hObject,eventdata,handles);

% --- Executes on button press in individual.
function individual_Callback(hObject, eventdata, handles)
plot_data(hObject,eventdata,handles);


