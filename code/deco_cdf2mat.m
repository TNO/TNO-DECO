function [matrix,startmass,starttime]=deco_cdf2mat(cdf_file,massfactor);
%-------------------------------------------------------------------
% transform the linear signal of netcdf to a matrix
% written by R.J. Jellema
% adapted by J.T.W.E. Vogels 
%-------------------------------------------------------------------
% last modification 05-07-2006
%-------------------------------------------------------------------
% input:      nc2matfile the name of the input file
% info:       structured information matrix
%-------------------------------------------------------------------
waitsigns=['*'];

multifactor = 1/massfactor; % correction for LCMS data (massfactor = 1.001)

fname = cdf_file;

[dir,name]=fileparts(cdf_file);
str = ['reading file:' name];


[V_scan_index]= cdf_extract_d([cdf_file '.cdf'],'scan_index');      
[V_intensity_values V_mass_values D_scan_number D_point_number V_time]= cdf_extract_s([cdf_file '.cdf'],'intensity_values','mass_values','scan_number','point_number','scan_acquisition_time');      

max_mass=round(multifactor*max(V_mass_values));
min_mass=round(multifactor*min(V_mass_values(find(V_mass_values>0))));

num_masses = (max_mass - min_mass)+1; % added by JV 12/7/2005 to replace full matrix
matrix=zeros(D_scan_number,num_masses,'single');

for scan=1:double(D_scan_number), % for all scans
    if (scan < D_scan_number) 
        to_scan = V_scan_index(scan+1); % get from current to next scan number
    else
        to_scan = single(D_point_number); % get last scan number
    end
    
    if (to_scan>size(V_mass_values,1)) % bug fix in single conversion fixed JV 5-7-2006
        to_scan=size(V_mass_values);
    end
   
    if V_scan_index(scan) ~= to_scan  % added check on empty scans JV 12/7
        massas_in_scan = V_scan_index(scan)+1:to_scan;   
        idx = find(V_mass_values(massas_in_scan)==0);
      
        massas_in_scan(idx)=[];
        mss = round(multifactor*V_mass_values(massas_in_scan))';          
        % deco_unique  fixes the problem of multiple peaks per mass
        % 3/8/2006 J.V.
        [v,w] = deco_unique_mass(mss,V_intensity_values(massas_in_scan));
        nm(1) = length(mss);
        nm(2) = length(V_intensity_values(massas_in_scan));
        nm(3) = length(v);
        nm(4) = length(w);
        
        %matrix(scan,[round(multifactor*V_mass_values(massas_in_scan))-min_mass+1])=V_intensity_values(massas_in_scan)';
        matrix(scan,[v-min_mass+1])=w';      
        %k= nm';
    end
    %waitbar(scan/double(D_scan_number),ht);
end


%close(ht); % close waitbar

% last signal is 0, remove this;
deco_keep3 matrix min_mass V_time; % remove all data except matrix and min_mass
  % pack the data 

startmass=min_mass;
starttime=V_time;

