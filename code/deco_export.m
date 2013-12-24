function y = deco_export()
%--------------------------------------------------------------------
% Export the areas and peak information as afinal step from deco
% adapted from deco_export + deco_clean by R. Jellema
% corrected output 23/11/2008 
% output is now exactly the same as in deco classic approach
%--------------------------------------------------------------------
% Written 23/08/2008 J. Vogels
%-------------------------------------------------------------------
global project;

clc; % clear screen

pw = project.pw;         % project average peak width
f = project.interval;    %#ok<NASGU> % project average interval
y = 1;

globt = []; % global times 
globa = []; % global areas
globs = []; % global indexes of files to include
names = []; % global peak names

f = project.interval;
n = 0;

% results for all spectra are collected 
% using the peakpositions and area's of all spectra (model and test)
% in a future version this may be replaced by another version where the
% selection of which peaks to print is limited to only the model spectra
% for now we left the complete version in. J.V. 1-08-2008
lastb = project.last;

disp(['limited to ' num2str(lastb) ' blocks'])
for block=project.first:lastb %project.last
   
    if (isempty(project.deco{block} ) == 0) % block is not empty then process
         N = size(project.deco{block}.copt,2); % number of peaks per block
    else 
        N=0;
    end
   
    if (N>0) % ignore empty blocks
        p = zeros(N,project.nfiles+project.ntestfiles); 
        rt = zeros(N,project.nfiles+project.ntestfiles);
        for file=1:project.nfiles % for all files in project
            % new version
            pos = project.deco{block}.pmax(:,file)';
            rt(:,file) = (f*(pos+(block-1)*2*pw)+project.rt_start(file))/60;
            p(:,file) = pos;
            % pos is the position of the line in the current block
            
        end
       
        for file=1:project.ntestfiles % collect th results from the test files
            pos = project.deco{block}.tmax(:,file);
            rt(:,project.nfiles+file) = (f*(pos+(block-1)*2*pw)+project.rt_start(file))/60;
            p(:,project.nfiles+file) = pos; 
        end     
         
        % hot zone is pw - 3*pw peaks outside this region are cold !
        idx1 = find(p<project.pw);    % find elements with pos < pw
        if (block<project.last) %do include peaks from end of last block
            idx2 = find(p>=3*project.pw); % find elements with pos > 3* pw
        else
            idx2 = find(p>=3.9*project.pw);% last block goes to 90% of width of last block
        end
         
        if (project.ntestfiles > 0)
            area = [project.deco{block}.areaopt project.deco{block}.pred];
        else
            area = project.deco{block}.areaopt;
        end
        area(idx2)=0; 
        area(idx1)=0; % set areas for exclude peaks to zero
        rt(idx2)=0;   rt(idx1)=0;     % remove peaks outside central region
       
        p(idx1)=0; % set peaks below area to zero
        p(idx2)=0; % set peaks above area to zero
                      
        % we are going to change this in future to include peaks from same
        % block and remove peaks from different blocks where they overlap
        % i.e. where they have the same spectrum
        pc = sum(rt,2);               % sum retention time per peak
        idx = find(pc ~= 0);          % get the indices for peaks where rt is valid
        globs{block} = idx;           % store the indexes of peaks to include
        globt= [globt; rt(idx,:)];    % store the accepted retention times
        project.deco{block}.incl = zeros(N,1); % set peak inclusion to zero
        project.deco{block}.gpos = zeros(N,1); % store reference to global area matrix
                
        present = sum(rt>0,2);        % number of files with valid retention time
        st = size(globa,1);           % how big is the global block at the moment
        % mark valid peaks with number of valid files per peak        
        present(idx);
        project.deco{block}.incl(idx)=present(idx);
        for j=1:length(idx)
            project.deco{block}.gpos(idx(j)) = st+j ;  % store pointer to global area block
            n=n+1; % counter on number of names (in sync with globa globt etc)
            names{n} = project.deco{block}.pnames(idx(j)) ;
        end 
        globa = [globa ; area(idx,:)]; % collate valid area's        
    end
end

exclude = [];
for j=project.first:lastb-1%project.last-1      
   
    if isempty(project.deco{j}) ==0  
        inc = project.deco{j}.incl;
    else 
        inc = [];
    end
        
    for a=1:length(inc)   
        
       if (inc(a)>0 && inc(a)<project.nfiles+project.ntestfiles && isempty(project.deco{j+1})==0)  % not all peaks found now match with next block
           %disp ('peaks available in next block')
           g = project.deco{j+1}.incl; % which peaks are available 
           g1 = find(g>0); % these are the potential candidates in next block
           r = zeros(length(g),1); % initalize array                   
           for k=1:length(g1)
               r(g1(k))=deco_spec_angle(project.deco{j}.sopt(a,:),project.deco{j+1}.sopt(g1(k),:));    
           end
           
           [mx mp] = max(r); % get maximum overlap 
           
           if (mx>0.95)  % if > 0.95 then combine lines
               if (j+1)<lastb
                   project.deco{j+1}.incl(mp)=0;  
               end; % set include for other peak to zero
               disp(sprintf('combine %d block %d with %d from block %d (%d,%d)',a,j,mp,j+1,project.deco{j}.gpos(a),project.deco{j+1}.gpos(mp))); 
               if (j+1)<lastb
                   exclude = [exclude project.deco{j+1}.gpos(mp)];  %#ok<AGROW>
               end; 
               %project.deco{j+1}.sopt(mp,:) = 1; % this is done in classic
               gp1 = project.deco{j}.gpos(a);
               gp2 = project.deco{j+1}.gpos(mp);
               if (j+1)<lastb 
                   for k=1:project.nfiles
                        if globa(gp1,k) == 0  
                             %disp(sprintf('%d %d %e',gp1,k,globa(gp2,k)))
                             globa(gp1,k) = globa(gp2,k);
                             globt(gp1,k) = globt(gp2,k); 
                        end
                   end
               end
           end
       end
   end     
end

glob_mean = sum(globt,2) ./ sum(globt>0,2); % calculate average RT per compound (corrected for occurance)

disp('temporary disabled raw spectra')
% 
% %--------------------------------------------------------------------------
% % write raw spectra export (when requested) (checked and OK 23/11/2008 JV)
% %--------------------------------------------------------------------------
% n=0;
% if (project.full_box ~= 0)
%     fid = fopen([project.name,'_FULL.msp'],'w+t');
%     disp('Write full msp')
%     for j=1:size(globs,2)
%         idx = globs{j};
%         for i=1:length(idx)
%             spec = project.deco{j}.sopt(idx(i),:);
%             [a,b] = find(spec>0.00001); % masses above threshold
%             n=n+1;
%             if ~isempty(a),        
%                 fprintf(fid,'Name Rt= %.4f Peak %d %s\n',glob_mean(j),n,project.deco{j}.pnames{idx(i)});
%                 fprintf(fid,'NumPeaks %d\n',length(b));
%                 fprintf(fid, '%f %f\n', [project.minmass+project.theMasses(b)-1; spec(b)]);
%             else
%                 numberofpeaks=1;
%                 fprintf(fid,'Name: Rt %.4f Peak: %d not usable peak\n',glob_mean(j),n);
%                 fprintf(fid,'%s %f \r\n','Num Peaks:',numberofpeaks);
%                 fprintf(fid, '1 1\n');
%             end
%         end
%     end
%     fclose(fid);
% end
% % Clean up the specta
% globa(exclude,:) = [];
% globt(exclude,:) = []; % global retention times
% names(exclude) = [];
% glob_mean = sum(globt,2) ./ sum(globt>0,2); % calculate average RT per compound

% disp('diasbles cleaned spectra')
%-----------------------------------------------
%Write cleaned spectra (when requested)
%-----------------------------------------------
if (project.clean_box ~= 0)
    disp('write clean data')
    fid = fopen([project.name,'clean_.msp'],'w+t');
    n=0;
    for j=1:size(globs,2)
        idx = globs{j};
        for i=1:length(idx)
            if project.deco{j}.incl(idx(i)) ~= 0
                spec = project.deco{j}.sopt(idx(i),:);
                n=n+1;
                [a,b] = find(spec>0.00001); % masses above threshold
                if ~isempty(a),        
                    fprintf(fid,'Name Rt= %.4f Peak %d %s\n',glob_mean(j),n,char(names{n}));
                    fprintf(fid,'NumPeaks %d\n',length(b));
                    fprintf(fid, '%f %f\n', [project.minmass+project.theMasses(b)-1; spec(b)]);
                end
            end
        end
    end
    fclose(fid);
end

%-----------------------------------------------
% Export global area data (when requested)
%-----------------------------------------------
fres = fopen([project.name,'_AREA.xls'],'w+t');
fprintf(fres,'nr\t rt(min)\t ident');
for i=1:project.nfiles 
    [path,name] = fileparts(project.files{i});
    fprintf(fres,'\t%s',name);
end
for i=1:project.ntestfiles 
    [path,name] = fileparts(project.testfiles{i});
    fprintf(fres,'\t%s',name);
end

fprintf(fres,'\n');
[val,idx] = sort(glob_mean);

for i=1:size(globa,1)
    if project.sort_output==1 
        id = idx(i);
    else
        id = i; %idx(i); % don't sort entries
    end
    fprintf(fres,'%03d',i); % pieknummer
    fprintf(fres,'\t%6.2f',glob_mean(id));
    fprintf(fres,'\t%10s',char(names{id}));
    for j=1:project.nfiles
        fprintf(fres,'\t%10.3E',globa(id,j));
    end
    for j=1:project.ntestfiles
        fprintf(fres,'\t%10.3E',globa(id,j+project.nfiles));
    end
    fprintf(fres,'\n');
end
fclose(fres);

disp('temporary disable global retention times')
% %-----------------------------------------------
% % Export global retention time data
% %-----------------------------------------------
% if (project.retention_box ~=0)
%     fres = fopen([project.name,'_RT.xls'],'w+t');
%     fprintf(fres,' nr\trt(min)\tident');
%     for i=1:project.nfiles
%         [path,name] = fileparts(project.files{i});
%         fprintf(fres,'\t%s',name);
%     end
%     fprintf(fres,'\n');
%     [val,idx] = sort(glob_mean);
%     for i=1:size(globa,1)
%         id = i; %idx(i);
%         fprintf(fres,'%03d\t',i); % pieknummer
%         fprintf(fres,'%6.3f\t',glob_mean(id));
%         fprintf(fres,'%10s\t',char(names{id}));
%         for j=1:project.nfiles
%             fprintf(fres,'%.4f\t',globt(id,j));
%         end
%         fprintf(fres,'\n');
%     end
%     fclose(fres);
% end

disp('temporary disabled individual files')
% %------------------------------------------
% % export integrals of individual files
% %------------------------------------------
% if (project.individual_box ~= 0)
%     for i=1:project.nfiles
%         [path,name] = fileparts(project.files{i});
%         fp = fopen([name ,'_are.txt'],'w+t');
%         fprintf(fp,'Nr. Rt(min)  Area     %%Area\n');
%         sms = sum(globa(:,i)); % sum off all integrals
%         [val,idx] = sort(globt(:,i)); % sort by retention time
%         for j=1:size(globa,1)
%             if (globa(idx(j),i)>0.0)
%                 % the peak number is for compatability to msp file 
%                 % intentionally left to point to the original peak 
%                 fprintf(fp,'%03d %5.3f ',idx(j),globt(idx(j),i));
%                 fprintf(fp,'%10.3E %6.2f%%\n',globa(idx(j),i),globa(idx(j),i)*100.0/sms);
%             end
%         end
%         fclose (fp);
%     end
% end
% fclose('all');
