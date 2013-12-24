function varargout = deco_inspect(varargin)
% DECO_INSPECT M-file for deco_inspect.fig
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help deco_inspect
% Last Modified by GUIDE v2.5 23-Nov-2008 18:28:27
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deco_inspect_OpeningFcn, ...
                   'gui_OutputFcn',  @deco_inspect_OutputFcn, ...
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

% --- Executes just before deco_inspect is made visible.----------------
function deco_inspect_OpeningFcn(hObject, eventdata, handles, varargin)
%-----------------------------------------------------------------------
global project;
global origset;
global gHandles;

% Choose default command line output for deco_inspect
handles.output = hObject;
handles.n = project.first;
handles.sel = 1;
handles.ticnum =1;
handles.twee=1;
handles.plotfirst=0;
origset = 0;
set(handles.oldnc,'String',num2str(0));
set(handles.area_button,'Value',1);
set(handles.error_button,'Value',0);
set(handles.symm_button,'Value',0);
set(handles.pos_button,'Value',0);
if size(project.deco{handles.n},1)>0 && size(project.deco{handles.n}.pnames,2)>0
    set(handles.name_edit,'String',flipud(char(project.deco{handles.n}.pnames(handles.sel))));
end
% Update handles structure
guidata(hObject, handles);
gHandles = handles;
deco_plot(hObject,eventdata,handles);

   
% --- Outputs from this function are returned to the command line.
function varargout = deco_inspect_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
%-----------------------------------------------------------------------
varargout{1} = handles.output;

% --- Executes on next file button press in block.
function block_Callback(hObject, eventdata, handles)
%-----------------------------------------------------
global project;
set(handles.oldnc,'String',num2str(0));
skip = get(handles.skip_blocks,'Value');
if (handles.n<project.last) % move to the next block
    handles.n = handles.n + 1;
    if skip &&  (project.numpeaks(handles.n)==0 || isempty(project.deco{handles.n}))
        while (project.numpeaks(handles.n)==0 || isempty(project.deco{handles.n})) && handles.n < project.last
           handles.n = handles.n + 1;
        end
    end
    deco_plot(hObject,eventdata,handles);     % plot the new data
    guidata(hObject,handles)
end

% --- Executes on button press in block_min.------------------------
function block_min_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------
global project;
set(handles.oldnc,'String',num2str(0));
skip = get(handles.skip_blocks,'Value');
if (handles.n>project.first)
    handles.n = handles.n - 1;
    if skip && (project.numpeaks(handles.n)==0 || isempty(project.deco{handles.n}))
        while (project.numpeaks(handles.n)==0 || isempty(project.deco{handles.n})) && handles.n>1
            handles.n = handles.n - 1;
        end
    end
    deco_plot(hObject,eventdata,handles);  
    guidata(hObject,handles)
end

%----------------------------------------------
function deco_plot(hObject,eventdata,handles) %#ok<INUSL>
%----------------------------------------------
global project;

if ~isempty(project.deco{handles.n} )
     npeaks = size(project.deco{handles.n}.copt,2);
else
    npeaks=0;
end

if handles.n>=1 && handles.n<=project.last 
    
    cmap = hsv(npeaks); % create a color map
    pw = project.pw;
    numfiles = project.nfiles;
    masses =project.theMasses;    
     
    set(handles.ncomp,'String',num2str(npeaks));
    th = str2double(get(handles.thresh,'String'));
    block = read_block(handles.n);
    
    % which file do we want to plot
    % the plot with the lowest regression between calc and original
    rsquared = zeros(project.nfiles,1);
    if (npeaks>0)
        for i=1:project.nfiles
          ds = sum(block((i-1)*4*project.pw+1:i*4*project.pw,masses),2);
          kk = project.deco{handles.n}.copt;
          ll = project.deco{handles.n}.sopt;
          ab = kk( (i-1) *4*pw+1:i*4*pw,:) * ll(:,masses);
          r= corrcoef(sum(ab,2),ds) .^ 2;
          rsquared(i) = r(1,2);
        end
    end
    
    %------------------------------------------
    % prepare quality parameters
    %------------------------------------------
    ppp=zeros(npeaks,numfiles);
    perr=zeros(npeaks,numfiles);
    psym=zeros(npeaks,numfiles);
    medpp = zeros(npeaks,1);
    for j=1:npeaks
        for i=1:numfiles
            ppp(j,i) = project.deco{handles.n}.quality(j,i).mxpos;
            perr(j,i) = project.deco{handles.n}.quality(j,i).err;
            psym(j,i) = project.deco{handles.n}.quality(j,i).symm;
        end
        medpp(j) = project.deco{handles.n}.pp(j);
    end
   
   %---------------------------------------------------------------
   % plot the main Title area
   %------------------------------------------------------------
   rt_mean = mean(project.rt_start);
   rt_st = ((handles.n-1)*2*project.pw*project.interval + rt_mean); % time in seconds
   rt_ed = rt_st + (4*project.pw*project.interval);
   if (npeaks>0)
        set(handles.model_text,'String',['block',num2str(handles.n),' of ' num2str(project.last) ' from:',num2str(rt_st/60),' to:',num2str(rt_ed/60) ' ni:' num2str(project.deco{handles.n}.ni) ]);
   else
        set(handles.model_text,'String',['block',num2str(handles.n),' of ' num2str(project.last) ' from:',num2str(rt_st/60),' to:',num2str(rt_ed/60) ' ni:' num2str(0) ]);
   end
   
   %-------------------------------------------------------
   %plot the individual peaks in topleft window
   %------------------------------------------------------
   xas = (rt_st + (1:4*project.pw) * project.interval)/60;
   axes(handles.model); % select topleft
   cla;
  
   
   if npeaks > 0
        hold on;
        [lowest,low] = min(rsquared);
        %select least regressed peak (to illustrate worst case)
        if (handles.ticnum>0) % first entry into routine
            low = handles.ticnum;
        else % when using mover buttons
            handles.ticnum=low;
        end
        s= project.deco{handles.n}.copt((low-1)*4*pw+1:low*4*pw,:);
        mmy = max(double(s));
        ylim([0 max(mmy)]*1.1);
        xlim([xas(1) xas(4*pw)]);
        hh = plot(xas,s);
    
        meeneem = ones(npeaks,1);
        for i=1:npeaks  
            %if project.deco{handles.n}.inc(i)==0
            %   set(hh(i),'color',cmap(i,:),'LineStyle','--'); 
            % draw peaks outside hot zone dotted (to be added)
            %else 
            set(hh(i),'Color',cmap(i,:));
            if sum(s(:,i))>0
                [err,symm,area,yp,sig] = deco_gaussfit(1:project.pw*4,s(:,i));
            else
                area=0;
            end
            %end
            hh2 = plot(xas,yp);
            set(hh2,'Color',cmap(i,:),'LineStyle','--');
            mx = project.deco{handles.n}.quality(i,low).mxpos;
            mxh = double (s(mx,i));
           
            if (area>0) 
                val1 = err / area * 100.0;
                val  = sig*6; % breedte v/d lijn in punten (2*3*sigma)
                val2 = sum(abs(yp'-s(:,i)))/area * 100;
               if (val2>30)
                    meeneem(i)=0;
                end
            else
                meeneem(i)=0;
            end
            
            
            text(xas(mx),mxh*1.05,[num2str(val1,'%.2f') ' ' num2str(val2,'%.2f') '(' num2str(meeneem(i),'%d') ')'],'HorizontalAlignment','Center','Color',cmap(i,:));

        end
        
        [path,name]= fileparts(project.files{low});
        line([xas(pw) xas(pw)],ylim,'color',[.8 .8 .8],'LineStyle','--');
        line([xas(pw*3) xas(pw*3)],ylim,'color',[.8 .8 .8],'LineStyle','--');
        line([xas(pw*2) xas(pw*2)],ylim,'color',[.8 .8 .8],'LineStyle','--');
        
        
        for i=1:npeaks
            line([(rt_st + project.deco{handles.n}.pp(i)*project.interval)/60 
                  (rt_st + project.deco{handles.n}.pp(i)*project.interval)/60],ylim,'color',cmap(i,:));
           
            mx = project.deco{handles.n}.quality(i,low).mxpos;
            mxh = double (s(mx,i));
%            
%              if get(handles.pos_button,'Value')==1
%                 val = project.deco{handles.n}.quality(i,low).mxpos-medpp(i); % positional difference from median position
%                 if val~=0
%                     text(xas(mx),mxh*1.05,num2str(val,'%.0f'),'HorizontalAlignment','Center','Color',cmap(i,:));
%                 end
%             end
%             if get(handles.symm_button,'Value')==1
%                 val = project.deco{handles.n}.quality(i,low).symm; % symmetry of the peak
%                 if val~=1
%                     text(xas(mx),mxh*1.05,num2str(val,'%.2f'),'HorizontalAlignment','Center','Color',cmap(i,:));
%                 end
%             end
%             if get(handles.error_button,'Value')==1 
%                 val = project.deco{handles.n}.quality(i,low).err; % fitting error 
%                 if val > 1
%                     text(xas(mx),mxh*1.05,num2str(val,'%.2f'),'HorizontalAlignment','Center','Color',cmap(i,:));
%                 end
%             end
%             
            
        end
                
        title(['(' num2str(low) ') ' name],'Interpreter','none');
        hold off; 
   end
         
   
    %------------------------------------------------
    % plot original and reconstructed spectra
    %-----------------------------------------------
    axes(handles.original); % middle left
    cla;
    xlim([xas(1) xas(4*pw)]);
    pw = project.pw;
    % original spectrum
    hold on;   
    ds = sum(block((handles.ticnum-1)*4*project.pw+1:handles.ticnum*4*project.pw,masses),2);
    plot(xas,ds,'k');
    meeneem;
    if (npeaks>0)
        ab = zeros(project.pw*4,1);
        for j=1:npeaks
            if meeneem(j)==1
                ab = ab + sum(project.deco{handles.n}.copt((low-1)*4*pw+1:low*4*pw,j) * project.deco{handles.n}.sopt(j,masses),2);
            end
        end
        %size(ab)
        plot(xas,ab,'r-.')
        legend('orig','reconstruct');
        text(xas(2*pw),max(ylim),['R^2:' num2str(rsquared(low),'%.3f')],'VerticalAlignment','Top','HorizontalAlignment','Center');
    end
    hold off;
    
    %-----------------------------------------------
    % DRAW RESIDUALS or individual mass traces
    %-----------------------------------------------
    axes(handles.residuals);
    cla; 
    hold on;
    
    %[highmasses,idxhighmasses]=sort(max(sum(block(:,project.theMasses)),1));
    i=handles.ticnum;
    spc = block((i-1)*pw*4+1:i*pw*4,:);
    
    z = zeros(size(spc,2),1);
    x=(1:4*project.pw)';

    for i=1:size(spc,2)
        K = polyfit(x,spc(:,i),2);
        v = spc(:,i) - polyval(K,x);
        z(i) = deco_signal2noise(v);
    end
    
    minsn = 6;
   
    title(['individual masses (sn >' num2str(minsn) ')'])
    hold off;
    axis tight;

    xlim([xas(1) xas(4*pw)]);
    idx = z>minsn; 
    if (length(find(idx)>0))
         plot(spc(:,idx))
    end
       
      
    %-------------------------------
    %plot the mass spectra
    %------------------------------
    axes(handles.calcspec);
    cla;
    %-----------------------------------------------------------------
    % count the number of peaks spectra have in common between them ?
    %-----------------------------------------------------------------
    hold on;
    if npeaks>0
        hold on;
        lastm =0;  
        for i=1:npeaks
            idx = find(project.deco{handles.n}.sopt(i,:)>th);
            if (max(idx)>lastm) 
                lastm = max(idx);
            end
        end
        if lastm
            maxm = masses(lastm)+project.minmass-1;
        end
        
        for i=1:npeaks;
            x = project.deco{handles.n}.sopt(i,:);
            mx = max(x); 
            if (mx) % removed plot bug mx=0 for empty spectra 9-6-2006 JV
                idx = find(x>th);
                x = x /mx + i;
                z = [];
                y = [];
                for j=1:length(idx)
                    z = [z masses(idx(j)) masses(idx(j)) 0]; %#ok<AGROW>
                    y = [y i x(idx(j)) nan]; %#ok<AGROW>
                end
                plot(z+project.minmass-1,y,'Color',cmap(i,:));
                plot([0 maxm],[i i],'Color',cmap(i,:)); % plot x axis       
                xp = double(masses(idx)+project.minmass-1);
                yp = double(x(idx));
                text(xp,yp,num2str(xp'));
                g = xlim;
                pp = project.deco{handles.n}.rt(i);
                text((g(1)+g(1))/2,i+0.8,[' RT=' num2str(pp) ' ' flipud(char(project.deco{handles.n}.pnames(i)))  ]);
            end
        end
        axis tight;
        title('mass spectra')
        hold off;
    end
    
    %---------------------------------------------------------------------
    % plot the areas and the concentration prediction of the model and
    % testset
    %--------------------------------------------------------------------
    axes(handles.predict); % bottom left
    cla;
    hold on;
    if get(handles.area_button,'Value')==1
        if npeaks>0
            if project.ntestfiles>0
                g = [project.deco{handles.n}.areaopt project.deco{handles.n}.pred];
            else
                g = project.deco{handles.n}.areaopt;
            end
            bar(g');
            colormap(cmap);
            title('sample concentration'); 
            xlim([0  size(g,2)+1]);
            ylim([0  max(max(g))*1.1]);  
        end
    end
    if get(handles.error_button,'Value')==1
        if npeaks>0
            g = perr;
            bar(g');
            colormap(cmap);
            title('area relative fit error'); 
            xlim([0  size(g,2)+1]);
            ylim([0  max(max(g))*1.1]);  
            line(xlim,[1.0 1.0],'Color',[1. 0. 0.],'LineStyle','--');
        end
    end
     if get(handles.symm_button,'Value')==1
        if npeaks>0
            g = psym;
            bar(g');
            colormap(cmap);
            title('peak symmetry error'); 
            xlim([0  size(g,2)+1]);
            ylim([0  max(max(g))*1.1]);  
            line(xlim,[0.8 0.8],'Color',[1. 0. 0.],'LineStyle','--');
        end
     end
     if get(handles.pos_button,'Value')==1
        if npeaks>0
            g = ppp-repmat(medpp,1,size(ppp,2));         
            bar(g');
            colormap(cmap);
            title('relative peak position (median)'); 
            xlim([0  size(g,2)+1]);
            ylim([min(min(g))*1.1  max(max(g))*1.1]);  
            line(xlim,[-3 -3],'Color',[1. 0. 0.],'LineStyle','--');
            line(xlim,[ 3 3],'Color',[1. 0. 0.],'LineStyle','--');
        end
     end
      if get(handles.regression_button,'Value')==1
        if npeaks>0
            bar(rsquared');
            colormap(cmap);
            %plot(rsquared)
            %title('Regression original / fit'); 
            xlim([0.5  project.nfiles+0.5]);
            ylim([min(rsquared)*0.9 1.0]);  
            line(xlim,[-3 -3],'Color',[1. 0. 0.],'LineStyle','--');
            line(xlim,[ 3 3],'Color',[1. 0. 0.],'LineStyle','--');
        end
      end
     line([handles.ticnum handles.ticnum],ylim,'color',[.8 .8 .8],'LineStyle','--');
    hold off;   
end

% --- Executes during object creation, after setting all properties.
function residuals_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------

%==============================
function y = read_block(i)
%==============================
global project
nf = project.nfiles;
y = deco_readblock([project.name '.bin'],i,project.pw,nf,project.ntraces,project.nmasses);
 
%-------------------------------------------------------
function thresh_Callback(hObject, eventdata, handles)
%-------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
%--------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------------
function ncomp_Callback(hObject, eventdata, handles)
%---------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function ncomp_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reprocess.----------------
function reprocess_Callback(hObject, eventdata, handles) %#ok<DEFNU>
%-----------------------------------------------------------
global origdeco
global project
global origset
global npeaks

warning off MATLAB:nearlySingularMatrix;
warning off MATLAB:rankDeficientMatrix;
npeaks = str2double(get(handles.ncomp,'String')); 

if (size(project.deco{handles.n},1)>0)
    nn = size(project.deco{handles.n}.copt,2);
else
    nn=0;
end

if npeaks ~= nn
      
    origdeco = project.deco{handles.n};
    set(handles.oldnc,'String',num2str(nn));
    origset=1;
    autoreduce = get(handles.AutoReduce,'Value');
    
    if (npeaks>0)
        deco_process_block(handles.n,npeaks,1,autoreduce);
        project.numpeaks(handles.n)=npeaks;
    else 
        project.deco{handles.n} = [];
        project.numpeaks(handles.n)=0;
    end
        
    deco_predict(handles.n);
    deco_block_evaluate(handles.n);
    deco_evaluate;
    deco_plot(hObject,eventdata,handles);
    deco_save_project;
 end

% --- Executes on button press in resetbutton.
%------------------------------------------------------------
function resetbutton_Callback(hObject, eventdata, handles)
global origdeco
global origset
global project
global npeaks
%------------------------------------------------------------
if origset
    nn = size(project.deco{handles.n}.copt,2);
    set(handles.oldnc,'String',num2str(nn));
    t=project.deco{handles.n};
    project.deco{handles.n}=origdeco;
    origdeco=t;
    nn = size(origdeco.copt,2);
    if handles.n>=1 && handles.n<=project.last
            project.numpeaks(handles.n)=npeaks;
    end       
    set(handles.oldnc,'String',num2str(nn));
    deco_plot(hObject,eventdata,handles);
end

%-------------------------------------------------------
function oldnc_Callback(hObject, eventdata, handles)
%-------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function oldnc_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in fast_deco.-----------------
function fast_deco_Callback(hObject, eventdata, handles)
%------------------------------------------------------------

function goto_block_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function goto_block_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in GotoButton.
function GotoButton_Callback(hObject, eventdata, handles)
global project;

n =  str2double(get(handles.goto_block,'String'));
if n>=1 && n<=project.last 
    handles.n = n;
    deco_plot(hObject,eventdata,handles);  
    guidata(hObject,handles);
end

% --- Executes on button press in nexttic.
function nexttic_Callback(hObject, eventdata, handles)
global project
if handles.ticnum < project.nfiles
    handles.ticnum = handles.ticnum +1;
    guidata(hObject,handles);
    deco_plot(hObject,eventdata,handles);
end

% --- Executes on button press in prevtic.
function prevtic_Callback(hObject, eventdata, handles)
if (handles.ticnum>1)
    handles.ticnum = handles.ticnum -1;
    guidata(hObject,handles);
    deco_plot(hObject,eventdata,handles);
end


% --- Executes on button press in next_example.
function next_example_Callback(hObject, eventdata, handles)
global project
if handles.twee < project.nfiles
    handles.twee = handles.twee+1;
    deco_plot(hObject,eventdata,handles);
end

% --- Executes on button press in prev_example.
function prev_example_Callback(hObject, eventdata, handles)
if handles.twee>1
    handles.twee = handles.twee-1;
    deco_plot(hObject,eventdata,handles);
end


% --- Executes on button press in area_button.
function area_button_Callback(hObject, eventdata, handles)
set(handles.area_button,'Value',1);
set(handles.error_button,'Value',0);
set(handles.symm_button,'Value',0);
set(handles.pos_button,'Value',0);
set(handles.regression_button,'Value',0);
deco_plot(hObject,eventdata,handles);


% --- Executes on button press in error_button.
function error_button_Callback(hObject, eventdata, handles)
set(handles.area_button,'Value',0);
set(handles.error_button,'Value',1);
set(handles.symm_button,'Value',0);
set(handles.pos_button,'Value',0);
set(handles.regression_button,'Value',0);
deco_plot(hObject,eventdata,handles);


% --- Executes on button press in symm_button.
function symm_button_Callback(hObject, eventdata, handles)
set(handles.area_button,'Value',0);
set(handles.error_button,'Value',0);
set(handles.symm_button,'Value',1);
set(handles.pos_button,'Value',0);
set(handles.regression_button,'Value',0);
deco_plot(hObject,eventdata,handles);


% --- Executes on button press in pos_button.
function pos_button_Callback(hObject, eventdata, handles)
set(handles.area_button,'Value',0);
set(handles.error_button,'Value',0);
set(handles.symm_button,'Value',0);
set(handles.pos_button,'Value',1);
set(handles.regression_button,'Value',0);
deco_plot(hObject,eventdata,handles);


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes on mouse press over axes background.
function calcspec_ButtonDownFcn(hObject, eventdata, handles)
global project;
currPt=get(gca,'CurrentPoint');
y = floor(currPt(1,2));
if y<0 
    y=1;
end;
handles.sel=y;
guidata(hObject, handles);
set(handles.name_edit,'String',flipud(char(project.deco{handles.n}.pnames(y))));

function name_edit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function name_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in assign_button.
function assign_button_Callback(hObject, eventdata, handles)
global project;
name = get(handles.name_edit,'String');
project.deco{handles.n}.pnames{handles.sel} = name;
deco_plot(hObject,eventdata,handles);
deco_save_project;

% --- Executes on button press in estimate_button.
function estimate_button_Callback(hObject, eventdata, handles)
deco_estimate_tool(handles.n,0);


% --- Executes on button press in reco_button.
function reco_button_Callback(hObject, eventdata, handles)

% --- Executes on button press in regression_button.
function regression_button_Callback(hObject, eventdata, handles)
set(handles.area_button,'Value',0);
set(handles.error_button,'Value',0);
set(handles.symm_button,'Value',0);
set(handles.pos_button,'Value',0);
set(handles.regression_button,'Value',1);
deco_plot(hObject,eventdata,handles);

% --- Executes on button press in ident_tag.
function ident_tag_Callback(hObject, eventdata, handles)
if handles.sel>0
    deco_identify(handles.n,handles.sel);
end

% --- Executes on mouse press over axes background.
function model_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on mouse press over axes background.
function original_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on mouse press over axes background.
function residuals_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on mouse press over axes background.
function predict_ButtonDownFcn(hObject, eventdata, handles)
global project;
currPt=get(gca,'CurrentPoint');
x = round(currPt(1,1));
if x>0 && x<=project.nfiles
    handles.ticnum=x;
    deco_plot(hObject,eventdata,handles);
end
    
% --- Executes on key press with focus on name_edit and no controls selected.
function name_edit_KeyPressFcn(hObject, eventdata, handles)

function printfig(src,evnt)
evnt;


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
 deco_save_project;


% --- Executes on button press in AutoReduce.
function AutoReduce_Callback(hObject, eventdata, handles)

% --- Executes on button press in skip_blocks.
function skip_blocks_Callback(hObject, eventdata, handles)

