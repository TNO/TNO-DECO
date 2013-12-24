function varargout = deco_projectname(varargin)
% DECO_PROJECTNAME M-file for deco_projectname.fig
%      DECO_PROJECTNAME by itself, creates a new DECO_PROJECTNAME or raises the
%      existing singleton*.
%
%      H = DECO_PROJECTNAME returns the handle to a new DECO_PROJECTNAME or the handle to
%      the existing singleton*.
%
%      DECO_PROJECTNAME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECO_PROJECTNAME.M with the given input arguments.
%
%      DECO_PROJECTNAME('Property','Value',...) creates a new DECO_PROJECTNAME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deco_projectname_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deco_projectname_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help deco_projectname
% Last Modified by GUIDE v2.5 28-Jan-2006 15:13:04
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_projectname_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_projectname_OutputFcn, ...
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

function deco_projectname_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = 'Yes';
guidata(hObject, handles);
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
         case 'title'
          set(hObject, 'Name', varargin{index+1});
         case 'string'
          set(handles.text1, 'String', varargin{index+1});
         case 'filename'
           set(handles.projectname,'String',varargin{index+1});
        end
    end
end

FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);
    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);
load dialogicons.mat
IconData=questIconData;
questIconMap(256,:) = get(handles.figure1, 'Color');
IconCMap=questIconMap;
Img=image(IconData, 'Parent', handles.axes1);
set(handles.figure1, 'Colormap', IconCMap);
set(handles.axes1, ...
    'Visible', 'off', ...
    'YDir'   , 'reverse'       , ...
    'XLim'   , get(Img,'XData'), ...
    'YLim'   , get(Img,'YData')  ...
    );
set(handles.figure1,'WindowStyle','modal')
uiwait(handles.figure1);

function varargout = deco_projectname_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
varargout{2} = get(handles.projectname,'String');
delete(handles.figure1);

function pushbutton1_Callback(hObject, eventdata, handles)
handles.output = get(hObject,'String');
guidata(hObject, handles);
h = [get(handles.projectname,'String') '.prj'];
if exist(h,'file')
    a=questdlg('Save anyway','Error! Project Exists','Yes','No','Yes'); % Yes or No answer
    if a(1)=='N' 
        return;
    end
end
uiresume(handles.figure1);

function pushbutton2_Callback(hObject, eventdata, handles)
handles.output = get(hObject,'String');
guidata(hObject, handles);
uiresume(handles.figure1);

function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);
else
     delete(handles.figure1);
end

function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    handles.output = 'No';
    guidata(hObject, handles);
    uiresume(handles.figure1);
end      
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    

function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function projectname_Callback(hObject, eventdata, handles)

function projectname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


