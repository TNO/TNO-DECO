function varargout = deco_tic(varargin)
% DECO_TIC M-file for deco_tic.fig
%      DECO_TIC, by itself, creates a new DECO_TIC or raises the existing
%      singleton*.
%      H = DECO_TIC returns the handle to a new DECO_TIC or the handle to
%      the existing singleton*.
%      DECO_TIC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECO_TIC.M with the given input
%      arguments.
%      DECO_TIC('Property','Value',...) creates a new DECO_TIC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deco_tic_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deco_tic_OpeningFcn via varargin.
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help deco_tic
% Last Modified by GUIDE v2.5 13-Aug-2008 11:54:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_tic_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_tic_OutputFcn, ...
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

% --- Executes just before deco_tic is made visible.
function deco_tic_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.n = 1; % start with first file
% Update handles structure
guidata(hObject, handles);
readfile(handles);
set(handles.reconstruct_tag,'Value',1);
set(handles.flag_tag,'Value',1);
set(handles.peaks_tag,'Value',0);
set(handles.block_tag,'Value',0);
plot_func(hObject,eventdata,handles);


% --- Outputs from this function are returned to the command line.
function varargout = deco_tic_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%---------------------------------------------
function plot_func(hObject,eventdata,handles)
%---------------------------------------------
global scantimes;
global name;
global tic;
global project;
global max_wp;
global most_wp;

cla
hold on;
timeax = scantimes / 60;

[base1] = deco_convex_baseline(tic,most_wp);
%base2 = deco_blcor(tic,15); % default blcor correction

legendstr{1} = 'tic';
legendstr{2} = 'base';

plot(timeax,tic,'b');
plot(timeax,base1,'r');
pw = project.pw;

%------------------------------------------------------
% plot the calculated peaks under the original tic
%------------------------------------------------------    
show_peaks = get(handles.peaks_tag,'Value');
show_reco  = get(handles.reconstruct_tag,'Value');
show_mark  = get(handles.markers,'Value');
flag_set   = get(handles.flag_tag,'Value');
show_block = get(handles.block_tag,'Value');

glob_tic = zeros(length(timeax),1);
g = [];
gx=[];
gy =[];
np=0;
if project.diddeco
    
    %nn = project.last - project.first + 1; % determine number of active blocks
    cmap = hsv(project.last+1);     
    for i=project.first:project.last       
        %rt_start = (i-1)*2*pw*project.interval + project.rt_start(handles.n);
        st = (i-1)*2*pw;
        ax = timeax(st+1:st+4*pw);
        
        if ~isempty(project.deco{i})
            npeaks = size(project.deco{i}.copt,2);   
                  for j=1:npeaks % met alle pieken
                    
                    if flag_set ==1 && project.deco{i}.flag(j) ~=0
                        show = 0;
                    else
                        show=1;
                    end
                    
                    if project.deco{i}.inc(j) ~=0 && show
                        peak = sum(project.deco{i}.copt((handles.n-1)*4*pw+1:handles.n*4*pw,j)*project.deco{i}.sopt(j,:),2);
                        if (show_peaks)
                            plot(ax,peak,'Color',cmap(i,:));
                        end
                        [mval,mpos] = max(peak);
                        mval =double(mval);
                        mpos = double(mpos);
                        glob_tic(st+1:st+4*pw,:) = glob_tic(st+1:st+4*pw,:) + peak;             
                        if show_mark==1
                            line([ax(mpos) ax(mpos)],[0 mval] ,'Color',cmap(i,:));
                            name = flipud(char(project.deco{i}.pnames(j)));
                            v = num2str(ax(mpos));
                            ga = double(ax(mpos));
                            gb = double(mval);
                            np=np+1;
                            g{np} = [name ':' v];
                            gx(np) = ga;
                            gy(np) = gb;
                        end
                    end
               end           
        end
    end
end

if np > 0 &&  show_mark==1
    text(gx,gy,g,'rotation',90);
end

%------------------------------------------
% plot the reconstructed chromatogram
%------------------------------------------
if show_reco
    legendstr{3} = 'reco';
    plot(timeax,glob_tic,'g--')
end
title(['File:' project.files{handles.n} ' (' num2str(handles.n) ' of ' num2str(project.nfiles) ') Median Peakwidth:' num2str(most_wp) ' Max:' num2str(max_wp) ']'],'Interpreter','none');
ylim('auto')
legend(legendstr);

% plot the odd block indicators
if (show_block)
    st = 1;
    g = ylim;
    x = [timeax(st) timeax(st) inf];
    y = [g(1) g(2) inf];
    x1= [];
    y1= [];
    for i=1:project.last
         if (rem(i,2)==0) % odd blocks
             x = [x timeax(st + 2*i*pw) timeax(st+2*i*pw) inf timeax(st + 2*i*pw+4*pw) timeax(st+2*i*pw+4*pw) inf ] ;
             y = [y g(1) g(2) inf g(1) g(2) inf];            
         else
             text(double(timeax(st + 2*i*pw)),max(ylim),num2str(i),'VerticalAlignment','Top','HorizontalAlignment','Center');
             x1 = [x1 timeax(st + 2*i*pw) timeax(st+2*i*pw) inf timeax(st + 2*i*pw+4*pw) timeax(st+2*i*pw+4*pw) inf];
             y1 = [y1 g(1) g(2) inf g(1) g(2) inf];           
         end
    end  
    plot(x,y,'--','color',[.8 .8 .8]); % plot the division lines
    plot(x1,y1,'-.','color',[.8 .8 .8]); % plot the division lines

end
hold off

%-----------------------------------
function readfile(handles)
%----------------------------------
global project;
global name;
global scantimes;
global base;

global tic;
global xp;
global theMat;
global max_wp;
global most_wp;

name = project.files{handles.n};
filename = [project.sdir name];
dbcname = [strtok(filename,'.') '.dbc'];
load('-mat',dbcname,'theMat','scantimes');
tic = sum(theMat,2);
tic=tic(:);

base = deco_convex_baseline(tic,10);
base=base(:);
[xp,yp,wp] = deco_peakpick(tic-base,0,-1,5);

most_wp = median(wp);
max_wp = max(wp);
idx = find(wp < 15);
xp(idx) = [];


% --- Executes on button press in zoom_button.
function zoom_button_Callback(hObject, eventdata, handles)
zoom on;

% --- Executes on button press in full_button.
function full_button_Callback(hObject, eventdata, handles)
global scantimes;
timeax = scantimes / 60;
mi = min(timeax);
mx = max(timeax);
d = (mx-mi)/20;

axis([ min(timeax)-d max(timeax)+d 0 1e6]);
axis 'auto y'


% --- Executes on button press in prev_file.
function prev_file_Callback(hObject, eventdata, handles)
if handles.n > 1 
    handles.n=handles.n-1;
    guidata(hObject, handles);
    readfile(handles);
    plot_func(hObject,eventdata,handles);
end

% --- Executes on button press in next_file.
function next_file_Callback(hObject, eventdata, handles)
global project;
if handles.n < project.nfiles
     handles.n = handles.n + 1;
     guidata(hObject, handles);
     readfile(handles);
     plot_func(hObject,eventdata,handles);
end


% --- Executes on button press in reconstruct_tag.
function reconstruct_tag_Callback(hObject, eventdata, handles)
plot_func(hObject,eventdata,handles)


% --- Executes on button press in peaks_tag.
function peaks_tag_Callback(hObject, eventdata, handles)
plot_func(hObject,eventdata,handles)

% --- Executes on button press in block_tag.
function block_tag_Callback(hObject, eventdata, handles)
plot_func(hObject,eventdata,handles)

% --- Executes on button press in move_left.
function move_left_Callback(hObject, eventdata, handles)
%disp('move right') 
global scantimes;
z = xlim;
d = (z(2)-z(1))*0.9; % move 90 of previous position
if z(1)+d < scantimes(end)/60
    xlim([z(1)+d z(2)+d]); % only move if axis remains in scantimes range
end


% --- Executes on button press in move_right.
function move_right_Callback(hObject, eventdata, handles)
%disp('move left')
global scantimes;
z = xlim;
d = (z(2)-z(1))*0.9; % move 90 of previous position
if (z(2)-d) > scantimes(1)/60
    xlim([z(1)-d z(2)-d]); % only move if axis remains in scantimes range
end

% --- Executes on button press in markers.
function markers_Callback(hObject, eventdata, handles)
plot_func(hObject,eventdata,handles)

% --- Executes on button press in flag_tag.
function flag_tag_Callback(hObject, eventdata, handles)
plot_func(hObject,eventdata,handles)


