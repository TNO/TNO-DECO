function varargout = deco_align_tool(varargin)
% DECO_ALIGN_TOOL M-file for deco_align_tool.fig
% Edit the above text to modify the response to help deco_align_tool
% Last Modified by GUIDE v2.5 23-Jan-2006 16:31:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_align_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_align_tool_OutputFcn, ...
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

% --- Executes just before deco_align_tool is made visible.
function deco_align_tool_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = deco_align_tool_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in zoom_in.
function zoom_in_Callback(hObject, eventdata, handles)
zoom on;

% --- Executes on button press in full.
function full_Callback(hObject, eventdata, handles)
axis auto;

%---------------------------------------------
function read(hObject,eventdata,handles)
%----------------------------------------------
global project;
global spec;
global times;

spec= zeros(project.nfiles);
times=zeros(project.nfiles);
for i=1:project.nfiles
    name = project.files{i};
    dbcname = [ strtok(name,'.') '.dbc'];
    load('-mat',[project.sdir dbcname],'theMat_cdf','scantimes');
    spec{i} = sum(theMat_cdf,2);
    times{i}=scantimes; 
end

plot_figure(hObject,eventdata,handles);

% --- Executes on button press in ReadSet.
function ReadSet_Callback(hObject, eventdata, handles)
read(hObject,eventdata,handles);

%------------------------------------------------------
function plot_figure(hObject,eventdata,handles)
%-----------------------------------------------------
global project;
global spec;
global times;

if project.nfiles<=1
    return
end

n = str2double(get(handles.target,'String'));
if n>project.nfiles
    n=project.nfiles;
end
if (n<1)
    n=1;
end

ref = str2double(get(handles.refspec,'String'));
if ref>project.nfiles
    ref=project.nfiles;
end
if (ref<0)
    ref=1;
end

[p,q] = deco_corr(spec{ref},spec{n},50); % get cross correlation of spec1 and 2
[a,b] = max(p); % determine maximum correlation (a=value, b=pos)
corr = q(b); % the shift needed (left or right)

news = circshift(spec{n}',corr)'; % new spectrum points
times{n}(1:5)
newi = circshift(times{n},corr); % new retention times
newi(1:5)

cla;
hold on;
disp('align tool')
%plot(times{ref},spec{ref},'k'); % plot reference
plot(times{n},spec{n},'g'); % plot raw
%plot(times{n},news,'r');    % plot corrected

legend(project.nfiles,'ref','raw','corrected');
hold off;

% --- Executes on button press in plusbutton.
function plusbutton_Callback(hObject, eventdata, handles)
global project;
n = str2double(get(handles.target,'String'));
if (n<project.nfiles) 
    n=n+1;
    set(handles.target,'String',num2str(n));
end
plot_figure(hObject,eventdata,handles);

% --- Executes on button press in minbutton.
function minbutton_Callback(hObject, eventdata, handles)
n = str2double(get(handles.target,'String'));
if (n>1) 
    n=n-1;
    set(handles.target,'String',num2str(n));
end
plot_figure(hObject,eventdata,handles);

function refspec_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function refspec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function target_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function target_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


