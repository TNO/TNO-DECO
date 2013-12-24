function varargout = deco_exclude(varargin)
% DECO_EXCLUDE M-file for deco_exclude.fig
%      DECO_EXCLUDE, by itself, creates a new DECO_EXCLUDE or raises the existing
%      singleton*.
%
%      H = DECO_EXCLUDE returns the handle to a new DECO_EXCLUDE or the handle to
%      the existing singleton*.
%
%      DECO_EXCLUDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECO_EXCLUDE.M with the given input arguments.
%
%      DECO_EXCLUDE('Property','Value',...) creates a new DECO_EXCLUDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deco_exclude_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deco_exclude_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help deco_exclude
% Last Modified by GUIDE v2.5 18-Jan-2008 16:01:13
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_exclude_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_exclude_OutputFcn, ...
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

% --- Executes just before deco_exclude is made visible.
function deco_exclude_OpeningFcn(hObject, eventdata, handles, varargin)
global project;

handles.output = hObject;
w=  num2str(project.ex_block);
set(handles.exblock,'String',w);
w = num2str(project.ex_mass);
set(handles.ex_mass,'String',w);

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = deco_exclude_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in okbutton.
function okbutton_Callback(hObject, eventdata, handles)
global project;
w = get(handles.exblock,'String');
if ~isempty(w)
    project.ex_block=[];
    w=strrep(w,',',' ');
    w=strrep(w,';',' ');
    h= textscan(w,'%s');
    n = size(h{1},1);
    for i=1:n
        project.ex_block(i)=str2double(h{1}(i));
    end
end

w = get(handles.ex_mass,'String');
project.ex_mass=[];
if ~isempty(w)
    w=strrep(w,',',' ');
    w=strrep(w,';',' ');
    h= textscan(w,'%s');
    n = size(h{1},1);
    for i=1:n
        project.ex_mass(i)=str2double(h{1}(i));
    end
end
delete (handles.figure1);


function exblock_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of exblock as text
%        str2double(get(hObject,'String')) returns contents of exblock as a double

% --- Executes during object creation, after setting all properties.
function exblock_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ex_mass_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ex_mass_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


