function varargout = deco_select_dialog(varargin)
% DECO_SELECT_DIALOG M-file for deco_select_dialog.fig
%      DECO_SELECT_DIALOG by itself, creates a new DECO_SELECT_DIALOG or raises the
%      existing singleton*.
%
%      H = DECO_SELECT_DIALOG returns the handle to a new DECO_SELECT_DIALOG or the handle to
%      the existing singleton*.
%
%      DECO_SELECT_DIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECO_SELECT_DIALOG.M with the given input arguments.
%
%      DECO_SELECT_DIALOG('Property','Value',...) creates a new DECO_SELECT_DIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deco_select_dialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deco_select_dialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help deco_select_dialog

% Last Modified by GUIDE v2.5 08-Apr-2008 11:06:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_select_dialog_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_select_dialog_OutputFcn, ...
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

function deco_select_dialog_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = 'Yes';
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
         case 'title'
          set(hObject, 'Name', varargin{index+1});
         case 'string'
          set(handles.text1, 'String', varargin{index+1});
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
set(handles.baseline,'Value',0);
% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')
% UIWAIT makes deco_select_dialog wait for user response (see UIRESUME)

basecalc(hObject,eventdata,handles);
guidata(hObject, handles);
plot_data(hObject,eventdata,handles);

uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = deco_select_dialog_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in AcceptButton.
function AcceptButton_Callback(hObject, eventdata, handles)
global TS;
global project;
handles.output = get(hObject,'String');
%disp('Accept button')
lvls = str2double(get(handles.level,'string'));
modelset = [];
testset = [];
nt=0;
for i=1:lvls
    f=find (TS==i);
    modelset{i} = project.files{f(1)};
    for j=2:size(f)
        nt=nt+1;
        testset{nt} = project.files{f(j)};
    end
end
project.nfiles = size(modelset,2);
project.files = modelset;
project.ntestfiles = size(testset,2);
project.testfiles = testset;
project.build=1;

guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in Cancelbutton.
function Cancelbutton_Callback(hObject, eventdata, handles)
handles.output = get(hObject,'String');
guidata(hObject, handles);
%disp('Cancel button')
uiresume(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);
else
    delete(handles.figure1);
end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    guidata(hObject, handles);
    uiresume(handles.figure1);
end    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    

function npc_Callback(hObject, eventdata, handles)

function npc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function level_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.

function level_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function baseline_Callback(hObject, eventdata, handles)
basecalc(hObject,eventdata,handles);
plot_data(hObject,eventdata,handles);

function processButton_Callback(hObject, eventdata, handles)
recalc(hObject,eventdata,handles);
plot_data(hObject,eventdata,handles);

%-----------------------------------------------------
%BASECALC calculates a PCA on the whole data
%-----------------------------------------------------
function basecalc(hObject,eventdata,handles)
global project;
global scores;
global eigval;

project.nfiles = project.nfiles + project.ntestfiles;
project.files = [project.files project.testfiles];
project.ntestfiles =0;
project.testfiles =[];

set(handles.level,'String',num2str(floor(project.nfiles/3)));
set(handles.nsamples,'String',num2str(project.nfiles));
basecorr = get(handles.baseline,'Value');
deco_resample_set(project.files,project.sdir,project.nfiles);
nlength = 0;
for i=1:project.nfiles
    name = deco_get_name(i);
    a = whos('-file',name,'theMat');
    if i==1
        nlength = a.size(1);
    else
        nlength = min(nlength,a.size(1));
    end
end

xdata = zeros(project.nfiles,project.ntraces);
for i=1:project.nfiles
    name = deco_get_name(i);
    load('-mat',name,'theMat');
    ytic = sum(theMat(1:nlength,1:end),2);
    if (basecorr)
        ytic = blcor(ytic,15); % default choice
    end
    %baseline correction
    xdata(i,:) = ytic;
end
% 
xdata = deco_scale_jv(xdata,'MEAN');
[scores,b] = svd(xdata,0);
eigval = diag( b );
tvar = sum(eigval);
idx= find( (eigval*100/tvar) >= 3);
size(idx,1);
set(handles.npc,'String',num2str(size(idx,1)));
recalc(hObject,eventdata,handles);


%-------------------------------------------------------------
%RECALC: recalculation to determine the clustering of data
%-------------------------------------------------------------
function recalc(hObject,eventdata,handles)
global project;
global scores;
global hierscore;
global eigval;
global sel;
global TS;
global drawlvl;

ndim = str2double(get(handles.npc,'string'));
lvls = str2double(get(handles.level,'string'));

if ndim<2
    ndim=2;
end

if ndim> min(size(scores))
    ndim = min(size(scores))-1;
end

perc = sum(eigval(1:ndim)*100/sum(eigval));
set(handles.totvar,'String',num2str(perc));
set(handles.npc,'String',num2str(ndim));
set(handles.nsamples,'String',num2str(project.nfiles));

if lvls<2
   lvls=2; % must be at least 2 levels
end

if lvls>project.nfiles
    lvls = project.nfiles;
end

set(handles.level,'String',num2str(lvls));

Y = pdist(scores(:,1:ndim),'euclidean');
hierscore = linkage(Y,'average');

x = project.nfiles - lvls;
drawlvl = hierscore(x,3 );
TS = cluster(hierscore,'maxclust',lvls);

% find all samples in cluster 2
sel = zeros(lvls,ndim);
for i=1:lvls
    f=find(TS==i);
    sel(i,:)=scores(f(1),1:ndim);
end


%------------------------------------------------
% plot the result of the data analysis
%------------------------------------------------
function plot_data(hObject,eventdata,handles)
global scores;
global project;
global drawlvl;
global hierscore;
global sel;
global eigval;
axes(handles.pcaplot);
cla;
hold on;
plot(scores(:,1),scores(:,2),'o','MarkerSize',8);
% plot labels
for i=1:project.nfiles
    text(scores(i,1),scores(i,2),[' ' num2str(i)]);
end
% plot selected files
plot(sel(:,1),sel(:,2),'o','MarkerEdgeColor','k',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',8);                    
title('pca - meancenter')
% plot axis labels
p=eigval(1)*100/sum(eigval);
xlabel(['pc1:' num2str(p) '%']);
p=eigval(2)*100/sum(eigval);
ylabel(['pc2:' num2str(p) '%']);
hold off;

axes(handles.hierplot);
dendrogram(hierscore);
x=1:project.nfiles;
y = ones(1,project.nfiles) * drawlvl;
line(x,y,'Color','r','LineStyle','--');
