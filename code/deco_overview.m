function varargout = deco_overview(varargin)
% DECO_OVERVIEW M-file for deco_overview.fig
%      DECO_OVERVIEW, by itself, creates a new DECO_OVERVIEW or raises the existing
%      singleton*.
%
%      H = DECO_OVERVIEW returns the handle to a new DECO_OVERVIEW or the handle to
%      the existing singleton*.
%
%      DECO_OVERVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECO_OVERVIEW.M with the given input arguments.
%
%      DECO_OVERVIEW('Property','Value',...) creates a new DECO_OVERVIEW or raises the
%      existing singleton*.  Starting from the left_button, property value pairs are
%      applied to the GUI before deco_overview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deco_overview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help deco_overview

% Last Modified by GUIDE v2.5 04-Sep-2008 11:06:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_overview_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_overview_OutputFcn, ...
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

% --- Executes just before deco_overview is made visible.
function deco_overview_OpeningFcn(hObject, eventdata, handles, varargin)
global project;
global first_use; % first use

first_use = 1;
set(handles.corr_button,'Value',1);
set(handles.orig_button,'Value',1);
set(handles.ref_tag,'String',num2str(project.align_ref));

% Choose default command line output for deco_overview
handles.output = hObject;
set(handles.lin_tag,'Value',0);
set(handles.spec_tag,'Value',1);
set(handles.cur_spec,'String','1');


% Update handles structure
guidata(hObject, handles);
load_spec(hObject,eventdata,handles);
plot_func(hObject,eventdata,handles);


% --- Outputs from this function are returned to the command line.
function varargout = deco_overview_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% plot the data 
function plot_func(hObject,eventdata,handles)
global project;
global first_use;
global tics;
global gmx;

axes(handles.linview)
curr = str2double(get(handles.cur_spec,'String'));

show_corr = get(handles.corr_button,'Value');
show_orig = get(handles.orig_button,'Value');
%ref = str2double(get(handles.ref_tag,'String'))
colors = zeros(60,3);

f = 1.0/16;
for i=1:8; colors(i,:) = [ 0 0 i*f+0.5]; end;
for i=1:16; colors(i+8,:) = [0 i*f 1]; end;
for i=1:16; colors(i+24,:) = [i*f 1 1-i*f]; end;
for i=1:16; colors(i+40,:) = [1 1-i*f 0]; end;
for i=1:4; colors(i+56,:) = [1-i*f 0 0]; end;
colors(60,:) = [1 1 1];
cla;
hold on;

if (get(handles.lin_tag,'Value')==1)
    set(handles.linview,'Color',[0 0 0]); % black background
    for i=1:size(tics,1) % plot in linview mode
        for j=1:size(tics{i,1})
            c = floor(log(abs(tics{i,1}(j)) +1.0)/gmx * 59)+1;
            if (c>1)
                line([tics{i,3}(j) tics{i,3}(j)],[i-0.5 i+0.5],'Color',colors(c,:));
            end
        end
    end
else % plot in tic mode
    set(handles.linview,'Color',[1 1 1]); % white background
    n=0;
    f = [];
    if show_orig==1 
        n=n+1; f{n} = 'orig';
    end
    if (show_corr==1)
        n=n+1; f{n} = 'corr';
    end
    for i=1:project.nfiles+project.ntestfiles
        mx = max(tics{i,1});
        text(0,i+0.5,tics{i,5},'HorizontalAlignment','Right');
        if (show_orig==1)
            plot(i+tics{i,1}/mx,'Color',[0.5 0.5 0.5]); % plot the original data in gray
        end
        if (show_corr==1)
            if curr==i
                plot(i+tics{i,2}/mx,'r');
            else
                plot(i+tics{i,2}/mx,'b');                   % plot the shifted data
                
            end
            
        end
    end
    legend(f);
end

if (first_use==1)
    axis tight;
    first_use =0;
end


function load_spec(hObject,eventdata,handles)
global project;
global tics;
global gmx; % global max
%curr = str2double(get(handles.cur_spec,'String'));
n = project.nfiles + project.ntestfiles;
f = project.nfiles;
tics = cell(n,5);
ref = str2double(get(handles.ref_tag,'String'));

gmx=0;
for i=1:project.nfiles
   
    filename = [project.sdir project.files{i}];
    dbcname = [strtok(filename,'.') '.dbc'];
    load('-mat',dbcname,'theMat_cdf','scantimes','align');
    tics{i,1} = sum(theMat_cdf,2); % original data
    tics{i,2} = sum(circshift(theMat_cdf,align),2); % the corrected data
    tics{i,4} = align; % the alignment factor
    [a,b]= deco_corr(tics{i,1},tics{i,2},10);
    tics{i,5} = strtok(project.files{i},'.');
    mx = max(log(abs(tics{i,1})+1.0));
    if (mx>gmx)
        gmx=mx; 
    end; % determine highest max value
    tics{i,3} = scantimes/60; % the x - axis (times)
end
for i=1:project.ntestfiles
    filename = [project.sdir project.testfiles{i}];
    dbcname = [strtok(filename,'.') '.dbc'];
    load('-mat',dbcname,'theMat_cdf','scantimes','align');
    tics{f+i,1} = sum(theMat_cdf,2); % original data
    tics{f+i,2} = sum(circshift(theMat_cdf,align),2); % the corrected data
    tics{f+i,4} = align; % the alignment factor
    tics{f+i,5} = strtok(project.testfiles{i},'.');
    mx = max(log(abs(tics{i,1})+1.0));
    if (mx>gmx)
        gmx=mx; 
    end; % determine highest max value
    tics{f+i,3} = scantimes/60; % the x - axis (times)
end
if project.nfiles+project.ntestfiles>0
    set(handles.cur_spec,'String',num2str(1));
    set(handles.cur_shift,'String',num2str(tics{1,4}));
end

% --- Executes on button press in zoom_button.
function zoom_button_Callback(hObject, eventdata, handles)
zoom on;

% --- Executes on button press in full_button.
function full_button_Callback(hObject, eventdata, handles)
axis tight;

% --- Executes on button press in spec_tag.
function spec_tag_Callback(hObject, eventdata, handles)
set(handles.lin_tag,'Value',0);
disp('spec_tag');
plot_func(hObject,eventdata,handles);

% --- Executes on button press in lin_tag.
function lin_tag_Callback(hObject, eventdata, handles)
set(handles.spec_tag,'Value',0);
%disp('lin_tag')
plot_func(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function linview_CreateFcn(hObject, eventdata, handles)

function ref_tag_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ref_tag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in use_button.
function use_button_Callback(hObject, eventdata, handles)

% --- Executes on button press in plus_button.
function plus_button_Callback(hObject, eventdata, handles)
global project;
global tics;
curr = str2double(get(handles.cur_spec,'String'));
if (curr<project.nfiles+project.ntestfiles) 
    curr = curr+1;
    %load('-mat',dbcname,'theMat_cdf');
    set(handles.cur_spec,'String',num2str(curr));
    set(handles.cur_shift,'String',num2str(tics{curr,4}));
    plot_func(hObject,eventdata,handles);

end

% --- Executes on button press in min_button.
function min_button_Callback(hObject, eventdata, handles)
global tics;
curr = str2double(get(handles.cur_spec,'String'));
if (curr>1) 
    curr = curr-1;
    set(handles.cur_spec,'String',num2str(curr));
    set(handles.cur_shift,'String',num2str(round(tics{curr,4})));
    plot_func(hObject,eventdata,handles);
end
    

function cur_spec_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function cur_spec_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in left_button.
function left_button_Callback(hObject, eventdata, handles)
global project;
global tics;
current = str2double(get(handles.cur_spec,'String'));
if current>0 && current <= project.nfiles + project.ntestfiles
    cshift = double(tics{current,4});
    tics{current,2} = circshift(tics{current,1},cshift+1);
    tics{current,4} = cshift+1;
    set(handles.cur_shift,'String',num2str(round(tics{current,4})));
    plot_func(hObject,eventdata,handles);
end


% --- Executes on button press in right_button.
function right_button_Callback(hObject, eventdata, handles)
global tics;
global project;
current = str2double(get(handles.cur_spec,'String'));
if current>0 && current <= project.nfiles + project.ntestfiles
    cshift  = double(tics{current,4});
    tics{current,2} = circshift(tics{current,1},cshift-1);
    tics{current,4} = cshift-1;
    set(handles.cur_shift,'String',num2str(round(tics{current,4})));
    plot_func(hObject,eventdata,handles);
end

function cur_shift_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function cur_shift_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in orig_button.
function orig_button_Callback(hObject, eventdata, handles)
plot_func(hObject,eventdata,handles);

% --- Executes on button press in corr_button.
function corr_button_Callback(hObject, eventdata, handles)
plot_func(hObject,eventdata,handles);


