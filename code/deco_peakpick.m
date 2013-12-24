function [xp,yp,wp] = deco_peakpick(y,noiselvl,th,mult,dpos,block,outfile)
%##########################################################################
%#  function [xp, yp, wp] = peakpick(y,noise,th,mult,pos)
%# 
%# Purpose : produce a line list from a spectrum
%#
%# Input parameters :
%# y      input spectrum a vector
%# noise  estimated signal to noise (if noise <=0) it is estimated from the
%#        data using the noisejv algorithm (default -1)
%# th     threshold ( default 0.0) when th<0 th = 3* noise 
%# mult   noise multiplier (default 3)
%# dpos   if set only positive peaks are picked (default false)
%# block  size of the block to estimate noise (default 16 )
%# 
%# Output Parameters:
%# xp   list of x-coordinates of detected peaks
%# yp   list of heights of detected peaks
%# wp   list of estimated line widths of peaks
%# example pts  = peakpickjv(y,-1,0.0,5,false);
%##########################################################################
%# ------------------------------------------------------------------------
%# Adapted by J.T.W.E. Vogels 6/7/2004
%# ------------------------------------------------------------------------
%
if ~exist('noiselvl','var'); noiselvl = -1;   end;
if ~exist('th','var');    th = 0.0;     end;
if ~exist('dpos','var');  dpos = 0; end;
if ~exist('mult','var');  mult = 3;     end;
if ~exist('outfile','var'); outfile='peaks.txt'; end;
if ~exist('block','var'); block = 16; end;

xp=[];
yp=[];
wp=[];
    
y = real(y);           % only use real part of the spectrum to peak pick
y= y(:);               % make a column vector
if (size(find(y),1)<=1)
    return;
end

%y = (y / max(y)) *100; % scale highest value to 100 units

[m] = size(y,1);
   
if (noiselvl<=0)
    noiselvl = deco_bnoise(y,block);
end   

if (th<0)
    th = noiselvl*5; % if noise < 0 set default threshold level to 3 times signal to noise ratio
end;

nw = 1;
%index =0;
state = 0; % SEARCH = 0 UP =1 DOWN = 2 MIN = 3 MAX = 4
last = y(1);
noiselim = mult * noiselvl;

localmax = last;
localmin = last;
llamp = last;
localmaxp = nw;
localminp = nw;
llpos = nw;
np=0;

 while (nw < m) 
    point = y(nw);
    found = 0; 
   switch state
        case 0  % SEARCH
                if (point > last)   % SEARCH
                    if (point-localmin)>noiselim
                        state = 4; % max
                    end;
                    if point>localmax 
                        localmax = point;
                        localmaxp =nw;
                    end;
                else
                    if (localmax-point)>noiselim
                        state = 3; % min
                    end;
                    if point < localmin
                        localmin = point;
                        localminp=nw;
                    end;
                end;
        case 1 % UP
               if  (point<last)
                   if localmax < last && localmaxp > localminp
                       localmax  = last;
                       localmaxp = nw-1;
                   end;
                   if  (localmax-point) > noiselim
                       state = 3;
                       if (localmax > th) && (localmaxp ~= llpos)
                           found = 1;
                           llamp = localmax;
                           llpos = localmaxp;
                       end;
                       if localmaxp > localminp
                           localmin  = point;
                           localminp = nw;
                       end;
                   else
                    state  = 2;
                   end
               else
                   if (point - localmin) > noiselim
                       state = 4;
                     
                       if localmin < -th && dpos == 0
                           found = 1;
                           llamp = localmin;
                           llpos = localminp;
                       end;
                       if localmaxp < localminp
                            localmax  = point;
                            localmaxp = nw;
                       end;
                  end; 
              end;
        case 2  %DOWN         
            if (point>last) %#ok<ALIGN>
                if (localmin>last) && (localminp>localmaxp)
                    localmin = last;
                    localminp = nw-1;
                end;
                if (point -localmin) > noiselim
                    state = 4;
                    if localmin < -th && dpos == 0
                        found = 1;
                        llamp = localmin;
                        llpos = localminp;
                    end;
                    if (localmaxp<localminp)
                        localmax = point;
                        localmaxp=nw;
                    end;
                else 
                    state = 1;
                end;
            else % point <= last
             if (localmax - point) > noiselim
                state = 3;
                if (localmax>th) && (localmaxp ~= llpos)
                    found = 1;
                    llamp =localmax;
                    llpos = localmaxp;
                end;
                if (localmaxp > localminp)
                    localmin = point;
                    localminp=nw;
                end;
            end;
        end;
   
        case 4  %MAXS
                if point<last %#ok<ALIGN>
                if (localmax<last)
                    localmax = last;
                    localmaxp = nw-1;
                end;
                if (last-point)>noiselim
                    state = 3;
                    if last > th && nw-1 ~= llpos
                        found = 1;
                        llamp = last;
                        llpos = nw-1;
                    end;
                    if (localmaxp > localminp)
                        localmin = point;
                        localminp=nw;
                    end;
                else
                    state = 2;
                end;                                               
            end;
        case 3  %MINS
            if (point > last) %#ok<ALIGN>
                if (localmin>last)
                    localmin = last;
                    localminp = nw-1;
                end;
                if (point-last) > noiselim
                    state = 4;
                    if last<-th && nw-1 ~=llpos && dpos==0
                        found=1;
                        llamp = last;
                        llpos = nw-1;
                    end;
                    if (localmaxp<localminp)
                        localmax = point;
                        localmaxp = nw;
                    end;
                else
                    state = 1;
                end;
          end; % end min              
    end; % end switch
    
    if (found && abs(llamp) > th) 
        np = np +1;
        xp(np) = llpos;
        yp(np) = llamp;
    end;
    last = point;
    nw = nw + 1;
end;

%sm=deco_smooth(y,250);

% --------------------------------------------------
% find the halfheight width of the detected peaks
% --------------------------------------------------
 wp = ones(1,np); % get peak width at half height
 for i=1:np
     xloc = xp(i); % get position of the i-th peak
     hgt = y(xloc)/2;
     nl = 1; 
     while ( (xloc-nl)>=1 && y(xloc-nl)> hgt); nl = nl + 1; end; % search left end
     nr = 1;
     while ( (nr+xloc)<=m && y(xloc+nr)> hgt); nr = nr + 1; end; % search right end
     w = min(nr,nl);
     wp(i)=2*w;
 end
 

