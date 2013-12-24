function y = deco_block_evaluate(block)

%y=deco_retentiontime(block);


%--------------------------------------
% evaluate a block after processing
%--------------------------------------
% global project;
% pw= project.pw;
% f=project.interval;
% nfiles = project.nfiles;
% 
% disp('block evaluate')
% 
% N = size(project.deco{block}.copt,2); % number of compounds in this block
% peaks = zeros(N,nfiles);
% rt    = zeros(N,nfiles);
% pp   = zeros(nfiles,N);
%     
% for file=1:nfiles % for all files 
%      % which peaks lie in the active region of the block 1*pw<x>3*pw
%      [hgt,pos] =max(project.deco{block}.copt((file-1)*4*pw+1:file*pw*4,:));
%      pp(file,:) = pos;
%      [a,b] = find(pos>pw & pos<=3*pw);    % b peaks lie within range
%      peaks(b,file)=1;                     % set these peaks to 1
%      rt(:,file) = ((f * (pos+(block-1)*2*pw)+project.rt_start(file))/60); % get retention times
% end
% zt = rt'; % copy the retention times
% 
% 
% for l=1:N % for all peaks
%     if project.deco{block}.inc(l)==1
%         
%         pvals= pp(:,l); 
%         poss =sort(pvals); % sort the positions of the peaks in the files
%         iid2 = find(poss>=0); % reallly outside (below) the range
%         iid3 = find(poss<=4*pw); % really outside (above) the maximum range
%         iid = intersect(iid2,iid3); % active peaks inside the range
%         if (length(iid)>=1) % if there are any peaks
%         % evaluate based on position in batch points
%             mps = median(poss(iid)); % determine median value
%             stdval  = std(poss(iid));% determine std deviation of remaining peaks           
%             for ll=1:nfiles
%                 if (abs(pp(ll,l)-mps)>2*stdval) % peakposition outside acceptable values
%                     peaks(l,ll)=0; % then remove peak
%                     zt(ll,l) =0;   % remove from retention time array
%                 end  
%             end        
%         else 
%             for ll=1:nfiles 
%                peaks(l,ll)=0;
%                zt(ll,l)=0;
%             end
%         end
%         if block<project.lastblock
%             h=project.deco{block}.sopt(l,:);
%             for j=1:size(project.deco{block+1}.sopt,1)
%                 z=corrcoef(h,project.deco{block+1}.sopt(j,:)) .^ 2;
%                   if (abs(z(1,2))>0.95) 
%                     project.deco{block+1}.inc(j)=0;
%                     project.deco{block}.inc(l) = j+10;
%                         figure;
%                         hold on;
%                         plot(h,'g');
%                         plot(project.deco{block+1}.sopt(j,:),'r');
%                         z
%                         title(['block:' num2str(block) ' peak:' num2str(l) ' b2:' num2str(j) ' r:' num2str(z(1,2))])
%                         hold off;
%                 end
%             end
%         end
%     end   
%     
% end
%        
%     any_info = find(sum(peaks,2)); % which peaks contain any information in one of the spectra
%    
%     %calculate the average retention time of each peak 
%     %(excluding peaks outside hot area)
%     idx = find(zt>0); % find exisiting rts 
%     rtt=sum(zt',2)'; % su mvalues 
%     zt(idx)=1; % set these to 1
%        
%     nn =sum(zt); % found in how many spectra ?
%     idx1 = find(nn==0); nn(idx1)=1; % correct for non existend peaks
%     rt_mean=rtt ./ nn; % average retention time inside hot zone   
%    
%     % find overlap with next block
%     block
%     
%  
%     idx = find(peaks==0);
%     final_area = project.deco{block}.areaopt;
%     final_area(idx)=0;
%          
%     for i=1:length(any_info) % for each component   
%         nr=any_info(i); 
%         sp = sum(peaks(nr,:),2); % number of peak-in-spectra found
%         if (sp ~=project.nfiles && sp>0) % 
%             % if not all components are accounted for 
%             % in current block then look in first discarded peaks for other
%             % missing integrals (if peak quality is OK position difference not too
%             % big and peakshape ok
%             for j=1:nfiles 
%                 [hgt,pos] =max(project.deco{block}.copt((j-1)*4*pw+1:j*pw*4,nr));
%                 pp =((f * (pos+(block-1)*2*pw)+project.rt_start(file))/60);
%                 pd = abs(pp-rt_mean(nr))/rt_mean(nr)*100; % procent diff versus mean position
%                
%                 if (peaks(nr,j)==0 ) % for the missing peaks   
%                     if (pd<=1.0 && pos>=5 && pos<=3*pw+pw/2) % else the difference is too large
%                         final_area(nr,j) = project.deco{block}.areaopt(nr,j);
%                     end
%                 elseif (pd>1.0) % peaks differ too much (more than 1%)
%                         final_area(nr,j) =0.0;
%                 end
%             end
%         end
%     end
%     %rt_mean
%     %final_area 
%     
%     project.deco{block}.rt = rt_mean; % the actual accurate position of peaks in block
%     project.deco{block}.area = final_area;
%    
%     y=1;