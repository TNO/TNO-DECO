%------------------------------------------------------
% deco clean
% modified to include peak assignments 23/11/2008 JV
%------------------------------------------------------

disp('CLEAN')
global project;

Area_thisfile =[];
Rt_thisfile = [];
MS_thisfile=[];
Block_thisfile=[];
SpecInBloc_thisfile=[];

pw = project.pw;
f= project.interval;

SpecInBlock_thisfile=[];
Area_correct = [];
Copt_thisfile = zeros(project.nfiles,min(project.chrom_len));
TIC_thisfile = zeros(project.nfiles,min(project.chrom_len));

testarea = [];

names = [];
lastb = project.last;

disp([' = is limited to ' num2str(lastb) ' blocks'])
for block = project.first:lastb %project.last;
    block 
    if ~any(project.ex_block==block) && project.numpeaks(block)>0,
        mypeaks=zeros(size(project.deco{block}.copt,2),project.nfiles);
        myrt=mypeaks;
        for file=1:project.nfiles;              
            [height,idx]=max(project.deco{block}.copt((file-1)*4*pw+1:file*4*pw,:));
            
            [temp,savepeaks] = find(idx>=pw & idx<3*pw);
            mypeaks(savepeaks,file)=1;
            myrt(savepeaks,file)=(f*(idx(savepeaks)+(block-1)*2*pw)+project.rt_start(file))'/60;
            Block_thisfile_temp(savepeaks,file) = block;
            SpecInBlock_thisfile_temp(savepeaks,file) = savepeaks';
        end
               
        good=find(sum(Block_thisfile_temp,2)~=0);
        Block_thisfile=[Block_thisfile; Block_thisfile_temp(good,:)];
        SpecInBlock_thisfile=[SpecInBlock_thisfile; SpecInBlock_thisfile_temp(good,:)];
        clear Block_thisfile_temp SpecInBlock_thisfile_temp;
        any_info = find(sum(mypeaks')) % which peaks are meaningfull for any of the files
        
        no_info = find(~sum(mypeaks')); % which peaks are not meaningfull for all files  
        if any(any_info),
            %any_info
            %size(project.deco{block}.pred)
            
            % mypeaks 0/1  voor binnen of buiten pw<x<3*pw
            Area_thisfile = [Area_thisfile ; project.deco{block}.areaopt(any_info,:).*mypeaks(any_info,:)];
            names = [names project.deco{block}.pnames(any_info)];
           
            if project.ntestfiles>0
                temp = project.deco{block}.pred;
                idx = find(project.deco{block}.tmax < pw);
                idx2 = find(project.deco{block}.tmax>=3*pw);
                r = project.deco{block}.pred;
                r(idx)=0; r(idx2)=0;
                testarea = [testarea ; r(any_info,:)];
            end
            
            local_area = project.deco{block}.areaopt(any_info,:).*mypeaks(any_info,:);
            % wat doet dit ? Volgens mij is dit gewoon fout
            for u=1:length(any_info) 
                csum = sum(project.deco{block}.sopt(any_info(u),:));
                local_area(u,:) = local_area(u,:) .* csum;
            end
            % groot vraagteken
            
            Area_correct = [Area_correct ; local_area]; % stick area's together  
            Rt_thisfile   = [Rt_thisfile; myrt(any_info,:)];
            MS_thisfile   = [MS_thisfile; project.deco{block}.sopt(any_info,:)];
            %Block_thisfile = [Block_thisfile; ones(length(any_info),1)*block];
            for spec = 1:length(any_info),
                a = find(mypeaks(any_info(spec),:));
                for file = 1:length(a),
                    Copt_thisfile(a(file),(block-1)*2*pw+1:block*2*pw+2*pw) = project.deco{block}.copt((a(file)-1)*4*pw+1:a(file)*4*pw,any_info(spec))'+ ...
                        Copt_thisfile(a(file),(block-1)*2*pw+1:block*2*pw+2*pw);
                    TIC_thisfile(a(file),(block-1)*2*pw+1:block*2*pw+2*pw) = sum((project.deco{block}.copt((a(file)-1)*4*pw+1:a(file)*4*pw,any_info(spec))*...
                        project.deco{block}.sopt(any_info(spec),:)),2)'+ ...
                        TIC_thisfile(a(file),(block-1)*2*pw+1:block*2*pw+2*pw);
                end
            end
        end
    end
end

%========================================================
% RENGER METHOD 
%========================================================
for i=1:size(Area_thisfile,1);
    a=find(Rt_thisfile(i,:)~=0);
    Rt_mean(i)=mean(Rt_thisfile(i,a)');
end
% 
% fid = fopen([project.name,'_NIST_export_classic_full.msp'],'w');
% export_numbers=[1:size(Area_thisfile,1)];
% for j=1:length(export_numbers),
%     [a,b]=find(MS_thisfile(export_numbers(j),:) > 0.00001);
%     numberofpeaks=length(b);
%     if ~isempty(a),        
%         fprintf(fid,'%s %s %s %s %s\r\n','Name: Rt=', num2str(Rt_mean(j)),'Peak: ',num2str(j),names{j});
%         fprintf(fid,'%s %f \r\n','Num Peaks:',numberofpeaks);
%         temp_spec=MS_thisfile(export_numbers(j),b);
%         fprintf(fid, '%f %f\r\n', [project.minmass+project.theMasses(b)-1; temp_spec]);
%     else
%         numberofpeaks=1;
%         fprintf(fid,'%s %s %s %s %s\r\n','Name: Rt=', num2str(Rt_mean(j)),'Peak: ',num2str(j),'not usable peak');
%         fprintf(fid,'%s %f \r\n','Num Peaks:',numberofpeaks);
%         temp_spec=1;
%         fprintf(fid, '%f %f\r\n', [1; temp_spec]);
%     end
% end
% fclose(fid);
% equal number of entries in classic and normal retention times are
% different

%---------------------------------------------------------------------
% should be the same above this line (barring the exact rt position)
%---------------------------------------------------------------------
% copy the matrices to retain the originals

Area_thisfile_2=Area_thisfile;
Rt_thisfile_2=Rt_thisfile;
MS_thisfile_2=MS_thisfile;
Block_thisfile_2=Block_thisfile;
SpecInBlock_thisfile_2=SpecInBlock_thisfile;
for i=1:size(Block_thisfile,1),
    temp = find(Block_thisfile(i,:));
    if ~isempty(temp),
        Block_total(i) = Block_thisfile(i,temp(1));
    end
end


remove_results=[];
[a,b]=find(prod(Area_thisfile')==0); % find empty 'cells'

for i=1:length(a)-1,
    % find mass spectra in next block;
    matchMS = find(Block_total==(Block_total(b(i))+1));
  
    bi = b(i);
    hh = Block_total(b(i)+1);
      
    if ~isempty(matchMS),
        temp=zeros(2,project.nfiles);
        these = [b(i) matchMS]; % make matrix with the relevant spectra
        
        thespecsimm=[temp_calcsimm(double(MS_thisfile_2(these,:)))-eye(length(these))]; % calculate match  
        
        [e,f]=max(thespecsimm(1,:)); % get maximum match
                  
        if f~=1, % if the spectrum to be matched is only zeros, no match is possible (all NaN)
            maxcorr = matchMS(f-1);
        else
            e=0;
        end
       
        temp(1,find(Area_thisfile(b(i),:)~=0))=1;
        temp(2,find(Area_thisfile(maxcorr,:)~=0))=1;
        if max(e) >0.95;% &  sum(prod(temp)~=0),
            %disp(['row: ',num2str(b(i)), ' matches with: ',num2str(maxcorr)]);
            maxarea(1)=max(Area_thisfile(b(i),find(temp(1,:)==1)));
            maxarea(2)=max(Area_thisfile(maxcorr,find(temp(2,:)==1)));
            doublepeaks=find((prod(temp)==1));
            [mymax,maxidx]=max(maxarea); 
            if maxidx==1
                maxidx=b(i);
            else
                maxidx=maxcorr;
            end
            [mymin,minidx]=min(maxarea); 
            if minidx==1
                minidx=b(i);
            else
                minidx=maxcorr;
            end
           
        end
        
        if max(e) >0.95 && isempty(doublepeaks), % if there is a good match
            % find which files do have results for this spectrum
            matchcomp = find(Area_thisfile(b(i),:)==0);
%           find zeros in matching spectrum
            targetcomp = find(Area_thisfile(maxcorr,:)~=0);
            % find 'values' that correspond with 'zeros'
            a=intersect(matchcomp,targetcomp);
            if project.ntestfiles>0
                matchtest = find(testarea(b(i),:)==0); 
                targettest= find(testarea(maxcorr,:)~=0);
                 atest = intersect(matchtest,targetcomp);
            end       
            if ~isempty(a); % if any match fill in the zeros with matching values
                %disp('correct data')
                            
                Area_thisfile_2(b(i),a)=Area_thisfile_2(maxcorr,a);
                % correction of Retention times (added 18/12/2008 J.V.)
                Rt_thisfile_2(b(i),a) = Rt_thisfile_2(maxcorr,a); 
                
                if project.ntestfiles>0
                    testarea(b(i),atest)   = testarea(maxcorr,atest);
                end
                
                Area_correct(b(i),a)   = Area_correct(maxcorr,a);
                Block_thisfile_2(b(i),a)=Block_thisfile_2(maxcorr,a);
                SpecInBlock_thisfile_2(b(i),a)=SpecInBlock_thisfile_2(maxcorr,a);
                % add variable to the list of variables to be removed
                MS_thisfile_2(maxcorr,:)=1;
               
                remove_results = [remove_results maxcorr];
            end        
        end
    end
end

% clean the matrices, removing the problematic variables
Area_thisfile_2(remove_results,:)=[];
Area_correct(remove_results,:)=[];
Rt_thisfile_2(remove_results,:)=[];
MS_thisfile_2(remove_results,:)=[];
Block_thisfile_2(remove_results,:)=[];
SpecInBlock_thisfile_2(remove_results,:)=[];
names(remove_results) = [];

% recalculate the mean retention times
for i=1:size(Area_thisfile_2,1);
    a=find(Rt_thisfile_2(i,:)~=0);
    Rt_mean_2(i)=mean(Rt_thisfile_2(i,a)');
end

if (project.ntestfiles>0)
    testarea(remove_results,:) = [];
    testsizeremove= size(testarea);
    areasizeremove= size(Area_thisfile_2);
end

% % export the spectra to a NIST compatible file
% % checked 23/11/2008 JV
% fid = fopen([project.name,'_NIST_classic_clean.msp'],'w');
% export_numbers=[1:size(Area_thisfile_2,1)];
% for j=1:length(export_numbers),
%     [a,b]=find(MS_thisfile_2(export_numbers(j),:) > 0.00001);
%     numberofpeaks=length(b);
%     %fprintf(fid,'%s %s\r\n','Name: Rt=', num2str(Rt_mean(j)),' Cluster: ',num2str(Block_thisfile(j)));
%     if ~isempty(a),        
%         fprintf(fid,'%s %s %s %s %s\r\n','Name: Rt=', num2str(Rt_mean_2(j)),'Peak: ',num2str(j),names{j});
%         fprintf(fid,'%s %f \r\n','Num Peaks:',numberofpeaks);
%         temp_spec=MS_thisfile_2(export_numbers(j),b);
%         fprintf(fid, '%f %f\r\n', [project.minmass+project.theMasses(b)-1; temp_spec]);
%     else
%         numberofpeaks=1;
%         fprintf(fid,'%s %s %s %s %s\r\n','Name: Rt=', num2str(Rt_mean_2(j)),'Peak: ',num2str(j),'not usable peak');
%         fprintf(fid,'%s %f \r\n','Num Peaks:',numberofpeaks);
%         temp_spec=1;
%         fprintf(fid, '%f %f\r\n', [1; temp_spec]);
%     end
% end
% fclose(fid);

% create sample labels
slbls=[];
num_files=length(project.nfiles);
for i=1:project.nfiles;
    slbls=strvcat(slbls,strtok(project.files{i},'c.'));
end

data=struct('start_runs',project.rt_start,...
    'Block_thisfile_2',Block_thisfile_2,...
    'SpecInBlock_thisfile_2',SpecInBlock_thisfile_2,...
    'TIC_thisfile',TIC_thisfile,...
    'slbls',slbls,...
    'num_files',num_files,...
    'Area_thisfile_2',Area_thisfile_2,...
    'Rt_mean_2',Rt_mean_2);

save(['results_',project.name],'data');

[m,n] = size(Area_thisfile_2);

%------------------------------------------------------
% output changed 23/11/2008 to include peak assignments 
% checked against deco_export
%------------------------------------------------------
fres = fopen([project.name,'_classic_AREA.xls'],'w');
fprintf(fres,'nr\tmean\tident');
for i=1:project.nfiles
    [path,name] = fileparts(project.files{i});
    fprintf(fres,'\t%s',name);
end

for i=1:project.ntestfiles
     [path,name] = fileparts(project.testfiles{i});
     fprintf(fres,'\t%s',name);
end
fprintf(fres,'\n');

for i=1:m
    fprintf(fres,'%d',i); % pieknummer
    fprintf(fres,'\t%8.3f\t%s',Rt_mean_2(1,i),names{i});
    for j=1:n
        fprintf(fres,'\t%10.3E',Area_thisfile_2(i,j));
    end
    for j=1:project.ntestfiles
        fprintf(fres,'\t%10.3E',testarea(i,j));
    end
    fprintf(fres,'\n');
end
fclose(fres);

%-------------------------------------
% Export for all files the area's
% only export peaks with intensity > 0
%-------------------------------------
% for i=1:n % for all files export the area's separately
%     [path,name] = fileparts(project.files{i});
%     ffile = fopen([name,'_areas.txt'],'w+t');
%     sms = sum(Area_correct(:,i));
%     sms1=sum(Area_thisfile_2(:,i));
%         
%     fprintf(ffile,'%s\n',name);
%     fprintf(ffile,'Nr. Rt(min)  Area     %%Area\n');
%     z=0;
%     for j=1:m
%         if (Area_correct(j,i)>0.0)
%             z=z+1;
%             fprintf(ffile,'%03d %5.2f ',z,Rt_mean_2(1,j));
%             fprintf(ffile,'%10.3E %6.2f%%\n',Area_correct(j,i),Area_correct(j,i)*100.0/sms);
%         end
%     end
%     fclose(ffile);
% end

% export overall retention times 

%Rt_thisfile - globt1; % globt is equal to Rt_thisfile
%RT is correct (checked 23/11/2008 JV)
 rtfile = fopen([project.name,'classic_rt.xls'],'w+t');
 fprintf(rtfile,'nr\tmean\tident');
 for j=1:n 
     [path,name] = fileparts(project.files{j});
     fprintf(rtfile,'\t%s',name);
 end
 fprintf(rtfile,'\n');
 for i=1:size(Rt_thisfile_2,1)
     fprintf(rtfile,'%d\t%8.3f\t%s',i,Rt_mean_2(i),names{i});
     for j=1:size(Rt_thisfile_2,2)
        fprintf(rtfile,'\t%8.3f',Rt_thisfile_2(i,j));
     end
     fprintf(rtfile,'\n');
 end
fclose(rtfile);
