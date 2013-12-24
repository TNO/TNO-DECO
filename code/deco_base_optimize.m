function varargout = deco_base_optimize(varargin)
% DECO_BASE_OPTIMIZE M-file for deco_base_optimize.fig
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 01-Jun-2006 16:23:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_base_optimize_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_base_optimize_OutputFcn, ...
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

% --- Executes just before deco_base_optimize is made visible.-------------
function deco_base_optimize_OpeningFcn(hObject, eventdata, handles, varargin)
%--------------------------------------------------------------------------
global project;
global base;
% This function has no output args, see OutputFcn.
% Choose default command line output for deco_base_optimize
handles.output = 'Cancel';
handles.cf=1; % currentfile
handles.cm=1; % currentmass

% set value to values define in deco main
set(handles.current_lambda,'String',num2str(project.lambda));
set(handles.current_asym,'String',num2str(project.asym));
set(handles.current_order,'String',num2str(project.order));
set(handles.basewin,'String',num2str(project.basewin));
set(handles.slize_size,'String',num2str(project.gsize));
set(handles.slize_move,'String',num2str(project.gmove));
set(handles.slice_thresh,'String',num2str(project.gthresh));
set(handles.data_thresh,'String',num2str(project.thresh));

set(handles.original,'Value',1);
set(handles.baseline,'Value',1);
set(handles.corrected,'Value',0);

set(handles.BaseMethod,'Value',project.basemethod);
set_base_method(hObject,eventdata,handles);
set(handles.figure1,'WindowStyle','modal'); % create a modal window
%numblocks=floor(min(chrom_length)/(2*project.pw))-2
   
% Update handles structure
guidata(hObject, handles);
base =[]; %#ok<NASGU>
readfile(handles,1);
base =correct_base(hObject,eventdata,handles);
plot_func(hObject,eventdata,handles);
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.--------
function varargout = deco_base_optimize_OutputFcn(hObject, eventdata, handles) 
%-------------------------------------------------------------------------
varargout{1} = handles.output;
delete(handles.figure1);

%----------------------------------------------------------------
function current_lambda_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function current_lambda_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%---------------------------------------------------------------
function current_asym_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function current_asym_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%----------------------------------------------------------------
function current_order_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function current_order_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in massbox.-----------
function massbox_Callback(hObject, eventdata, handles)
%--------------------------------------------------------
global base;
global startmass;
m=get(handles.massbox,'Value');
handles.cm =m;
guidata(hObject,handles);
set(handles.slider_mass,'Value',m + startmass );
base =correct_base(hObject,eventdata,handles);
plot_func(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function massbox_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in zoomin_button.-------------
function zoomin_button_Callback(hObject, eventdata, handles)
%------------------------------------------------------------
zoom on;

% --- Executes on button press in full_button.-----------------
function full_button_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------
axis auto;

% --- Executes when figure1 is resized.-------------------
function figure1_ResizeFcn(hObject, eventdata, handles)
%---------------------------------------------------------

% --- Executes on button press in file_plus.--------------
function file_plus_Callback(hObject, eventdata, handles)
global project;
global base;
%---------------------------------------------------------
if handles.cf<project.nfiles
    handles.cf = handles.cf +1;    % increase file counter
    readfile(handles,handles.cf);  % read the data
    base =correct_base(hObject,eventdata,handles); % calculate a new baseline for this mass
    plot_func(hObject,eventdata,handles); % plot the result
    guidata(hObject,handles);
end

% --- Executes on button press in file_min.-----------------
function file_min_Callback(hObject, eventdata, handles)
global base;
%-----------------------------------------------------------
if handles.cf>1
    handles.cf = handles.cf -1;   % decrease the file counter
    readfile(handles,handles.cf); % read the data
    base =correct_base(hObject,eventdata,handles); % calcuate new baseline
    plot_func(hObject,eventdata,handles); % plot the result
    guidata(hObject,handles);
end

% --- Executes on button press in original.-------------
function original_Callback(hObject, eventdata, handles)
%-------------------------------------------------------
plot_func(hObject,eventdata,handles);

% --- Executes on button press in baseline.-------------
function baseline_Callback(hObject, eventdata, handles)
%-------------------------------------------------------
plot_func(hObject,eventdata,handles);

% --- Executes on button press in corrected.---------------
function corrected_Callback(hObject, eventdata, handles)
%----------------------------------------------------------
plot_func(hObject,eventdata,handles);

%------------------------------------------------
function plot_func(hObject,eventdata,handles)
%------------------------------------------------
global theMat;
global scantimes;
global startmass;
global name;
global base;
global tic;
global iq;
global project

ov      = get(handles.original,'Value'); % show original
bv      = get(handles.baseline,'Value');  % show fitted baseline
cv      = get(handles.corrected,'Value'); % show corrected orig-base
showb   = get(handles.show_blocks,'Value'); % show blocks
noise   = str2double(get(handles.noise,'String'));
pp      = get(handles.peakpick,'Value'); % show peakpicking
tic_but = get(handles.TIC_button,'Value'); % show tic spectrum
smooth_button = get(handles.smooth_button,'Value'); % show smoothed spectrum

cla;

hold on;
maxy = max(theMat(:,handles.cm)); % determine the highest peak
timeax = scantimes / 60 ;

pb = project.pw; % plot blocks
if showb
    for i=1:project.numblocks
        line([i*2*pb i*2*pb],[0 maxy],'Color','g')
    end
end

if (ov) % show original spectrum
    plot(timeax,theMat(:,handles.cm),'b'); % draw original spectrum
end

thr = str2double(get(handles.data_thresh,'String'));
if (thr<maxy)
    line(xlim,[thr thr]); % draw threshold
end

if bv % show calculated baseline
     % gives some times some nan 's 
    plot(timeax,base,'g'); % draw baseline
end

if cv  %show baseline corrected spectrum
    corr = theMat(:,handles.cm) - base;
    plot(timeax,corr,'k'); 
end

if (pp)% show peak picking on baseline corrected data
    [px,yh,wp]= deco_peakpick(theMat(:,handles.cm)-base,-1,-1,noise,1);
    plot(timeax(px),theMat(px,handles.cm),'bo');
end

if tic_but% show the tic spectrum
    tm = max(tic);
    t = max(theMat(:,handles.cm)) / tm; % same vertical display size as original specturm
    plot(timeax,tic .* t,'m');
end

if smooth_button % show smoothed spectrum
    sm=deco_smooth(theMat(:,handles.cm),8);
    plot(timeax,sm,'r'); % smoothed spectrum
end

mx = abs(maxy)/100;
plot([min(timeax) max(timeax)],-mx);

title([name ' mass:' num2str(handles.cm+startmass-1) ' iq:' ,num2str(iq(handles.cm))]); 
axis tight;
hold off


%--------------------------------------------------
function readfile(handles,i) % read original cdf file
%--------------------------------------------------
global project
global theMat
global filename;
global name;
global startmass;
global scantimes; %#ok<NUSED>
global base;
global tic;
global iq;

if i<1
    return; 
end

name = project.files{i};
filename = strcat(project.sdir, name);
dbcname  = strcat(strtok(filename,'.'), '.dbc'); % generate dbc name
load('-mat',dbcname,'theMat','scantimes','startmass'); % load from previous run

%size(theMat)
%v=min(theMat);
%size(v)
%theMat = theMat - repmat(v,size(theMat,1),1) ;

start = startmass;
m1 = startmass+1;
m2 = size(theMat,2)+startmass-1;
r = size(theMat,2)-2;
tic = sum(theMat,2); % generate tic spectrum

[m,n] = size(theMat);

iq=zeros(n,1);
for i=1:n
    sm=deco_smooth(theMat(:,i),20);
    ent = deco_entropy(sm);
    if (ent)
        iq(i) = 1 / ent; 
    end;
end

[miq,mpos] = max(iq); % find max iq value
iq = iq / miq; % normalize on largest iq value
[miq,mpos] = max(iq); % find max iq value

set(handles.slider_mass,'Value',m1,'Min',m1,'Max',m2,'SliderStep',[1.0/r 1.0/r]);
set(handles.massbox,'String',num2str([start:m2]'));
base =[];

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over massbox.
%function massbox_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on slider movement.---------------------------
function slider_mass_Callback(hObject, eventdata, handles)
%global theMat;
global base;
%global startmass;
%------------------------------------------------------------
m=round(get(hObject,'Value'))-get(hObject,'Min');
if (m>0) 
    handles.cm=m;
    set(handles.massbox,'Value',handles.cm)
end
guidata(hObject,handles);
% calculate the new baseline for this mass
base =correct_base(hObject,eventdata,handles);
plot_func(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function slider_mass_CreateFcn(hObject, eventdata, handles)
%------------------------------------------------------------------
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.----------------------------
function slider_order_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------
val=get(hObject,'Value');
set(handles.current_order,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function slider_order_CreateFcn(hObject, eventdata, handles)
%------------------------------------------------------------------
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.---------------------------
function slider_lambda_Callback(hObject, eventdata, handles)
%------------------------------------------------------------
%val=get(hObject,'Value')
set(handles.current_lambda,'String',num2str(10^val));

% --- Executes during object creation, after setting all properties.
function slider_lambda_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.--------------------------
function slider_asym_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------
val=get(hObject,'Value');
set(handles.current_asym,'String',num2str(10 ^ val));

% --- Executes during object creation, after setting all properties.
function slider_asym_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on selection change in BaseMethod.------------
function BaseMethod_Callback(hObject, eventdata, handles)
%------------------------------------------------------------
set_base_method(hObject,eventdata,handles); 
%base = correct_base(hObject,eventdata,handles);
plot_func(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
%----------------------------------------------------------------
function BaseMethod_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
%---------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------------------
function y = correct_base(hObject,eventdata,handles)
%------------------------------------------------------------------
global theMat;
options_bl.lambda  = str2double(get(handles.current_lambda,'String')); %#ok<ST2NM>
options_bl.p       = str2double(get(handles.current_asym,'String'));
options_bl.d       = str2double(get(handles.current_order,'String'));
options_bl.basewin = str2double(get(handles.basewin,'String'));
options_bl.gsize   = str2double(get(handles.slize_size,'String'));
options_bl.gmove   = str2double(get(handles.slize_move,'String'));
options_bl.gthresh = str2double(get(handles.slice_thresh,'String'));

g = get(handles.BaseMethod,'Value');
if g==4 % Convex hull baseline
    y = deco_convex_baseline(theMat(:,handles.cm),10)';
elseif g==3  % TNO method of baseline correction
    %disp('TNO');
    y = deco_blcor(double(theMat(:,handles.cm)),options_bl.basewin);
elseif g==2 % groningen method of baseline
    %disp('groningen');
    y = deco_medianbaseline(theMat(:,handles.cm),options_bl.gsize,options_bl.gmove,options_bl.gthresh);
elseif (g==1) % 
    % Eilers method
    y = single(deco_asysm(double(theMat(:,handles.cm)),options_bl.lambda,options_bl.p,options_bl.d));
else
    set(handles.corrected,'Value',0);
    y= [];
end   

%---------------------------------------------------------------
function basewin_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------

%---------------------------------------------------------------
function basewin_CreateFcn(hObject, eventdata, handles)
%---------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%---------------------------------------------------------------
function slize_size_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------

% ------------------------------------------------------------------
function slize_size_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------
function slize_move_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------

% -------------------------------------------------------------------
function slize_move_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-----------------------------------------------------------------------
function slice_thresh_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function slice_thresh_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------
function set_base_method(hObject,eventdata,handles) %#ok<INUSL,INUSL>
%------------------------------------------------------
set(handles.current_lambda,'Visible','Off');
set(handles.current_asym,'Visible','Off');
set(handles.current_order,'Visible','Off');
set(handles.lambda_text,'Visible','Off');
set(handles.asym_text,'Visible','Off');
set(handles.order_text,'Visible','Off');
set(handles.basewin,'Visible','Off');
set(handles.base_text,'visible','Off');  
set(handles.slize_size,'visible','Off');
set(handles.slize_move,'visible','Off');
set(handles.slice_thresh,'visible','Off');
set(handles.size_text,'visible','Off');
set(handles.move_text,'visible','Off');
set(handles.thresh_text,'visible','Off');

n = get(handles.BaseMethod,'Value');
if n==1
    set(handles.current_lambda,'Visible','On');
    set(handles.current_asym,'Visible','On');
    set(handles.current_order,'Visible','On');
    set(handles.lambda_text,'Visible','On');
    set(handles.asym_text,'Visible','On');
    set(handles.order_text,'Visible','On');
elseif n==2
    set(handles.slize_size,'visible','On');
    set(handles.slize_move,'visible','On');
    set(handles.slice_thresh,'visible','On');
    set(handles.size_text,'visible','On');
    set(handles.move_text,'visible','On');
    set(handles.thresh_text,'visible','On');
elseif n==3
    set(handles.basewin,'Visible','On');
    set(handles.base_text,'visible','On');       
end

% --- Executes on button press in tic.-----------------------
function tic_Callback(hObject, eventdata, handles)
%------------------------------------------------------------
plot_func(hObject,eventdata,handles); % redraw using new settings


% --- Executes on button press in peakpick.------------------
function peakpick_Callback(hObject, eventdata, handles)
%------------------------------------------------------------
plot_func(hObject,eventdata,handles); % redraw using new settings

% --- Executes on button press in usethis.-----------------
function usethis_Callback(hObject, eventdata, handles)
%----------------------------------------------------------
global project;
g = get(handles.BaseMethod,'Value');
project.basemethod=g;
if (g==3) % TNO Method
    project.basewin    = str2double( get(handles.basewin,'String'));
elseif (g==2) % Groningen method
    project.gsize      = str2double(get(handles.slize_size,'String'));
    project.gmove      = str2double(get(handles.slize_move,'String')); 
    project.gthresh    = str2double(get(handles.slice_thresh,'String'));
elseif (g==1) % Eilers method
    project.lambda     = str2double(get(handles.current_lambda,'String'));
    project.asym       = str2double(get(handles.current_asym,'String'));
    project.order      = str2double(get(handles.current_order,'String'));   
end
handles.output='OK';
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in Cancel.-------------------
function Cancel_Callback(hObject, eventdata, handles)
%----------------------------------------------------------
handles.output='Cancel';
guidata(hObject,handles);
uiresume(handles.figure1);

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------------

% --- Executes when user attempts to close figure1.------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------
if isequal(get(handles.figure1,'waitstatus'),'waiting')
    uiresume(handles.figure1);
else
    delete(hObject);
end

%----------------------------------------------------------------------
function noise_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function noise_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in smooth_button.-----------------
function smooth_button_Callback(hObject, eventdata, handles)
%-----------------------------------------------------------------
plot_func(hObject,eventdata,handles); % redraw using new settings

% --- Executes on button press in show_blocks.------------------
function show_blocks_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------
plot_func(hObject,eventdata,handles); % redraw using new settings

% --- Executes on button press in TIC_button.
function TIC_button_Callback(hObject, eventdata, handles)
plot_func(hObject,eventdata,handles); % redraw using new settings

function data_thresh_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function data_thresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


