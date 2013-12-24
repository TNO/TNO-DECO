function varargout = deco_global_params(varargin)
% DECO_GLOBAL_PARAMS M-file for deco_global_params.fig
%      DECO_GLOBAL_PARAMS, by itself, creates a new DECO_GLOBAL_PARAMS or raises the existing
%      singleton*.
%
%      H = DECO_GLOBAL_PARAMS returns the handle to a new DECO_GLOBAL_PARAMS or the handle to
%      the existing singleton*.
%
%      DECO_GLOBAL_PARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECO_GLOBAL_PARAMS.M with the given input arguments.
%
%      DECO_GLOBAL_PARAMS('Property','Value',...) creates a new DECO_GLOBAL_PARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deco_global_params_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deco_global_params_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help deco_global_params

% Last Modified by GUIDE v2.5 23-Nov-2008 18:14:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_global_params_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_global_params_OutputFcn, ...
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

% --- Executes just before deco_global_params is made visible.
function deco_global_params_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
global project;
set(handles.nist_dir,'String',project.nist_dir); %project.nist_dir);
set(handles.full_box,'Value',project.full_box);
set(handles.clean_box,'Value',project.clean_box);
set(handles.individual_box,'Value',project.individual_box);
set(handles.retention_box,'Value',project.retention_box);
set(handles.classic_box,'Value',project.classic);
%set(handles.sort_box,'Value',project.sort_output);
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = deco_global_params_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function nist_dir_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in nistbutton.
function nistbutton_Callback(hObject, eventdata, handles)
global project;
project.nist_dir=uigetdir('c:\');
set(handles.nist_dir,'String',project.nist_dir);

function nist_dir_Callback(hObject, eventdata, handles)

% --- Executes on button press in full_box.
function full_box_Callback(hObject, eventdata, handles)
global project;
project.full_box = num2str(get(hObject,'String'));

% --- Executes on button press in clean_box.
function clean_box_Callback(hObject, eventdata, handles)
global project;
project.clean_box = get(hObject,'Value');

% --- Executes on button press in individual_box.
function individual_box_Callback(hObject, eventdata, handles)
global project;
project.individual_box = get(hObject,'Value');

% --- Executes on button press in retention_box.
function retention_box_Callback(hObject, eventdata, handles)
global project;
project.retention_box = get(hObject,'Value');

% --- Executes on button press in classic_box.
function classic_box_Callback(hObject, eventdata, handles)
global project;
project.classic = get(hObject,'Value');

% --- Executes on button press in sort_box.
function sort_box_Callback(hObject, eventdata, handles)
global project;
%project.sort_output = get(hObject,'Value');

