function y = deco_predict(block)
global project
%-------------------------------------------------
% function to predict concentrations of a block
% of files contained in the project test set
% use spectra from earlier model calculations
%-------------------------------------------------

pw = project.pw;
if project.ntestfiles>0 
    s = sum(project.deco{block}.sopt,2); % calculate sum of masses in a spectrum 
    totblock = deco_readblock([project.name '.tst'],block,project.pw,project.ntestfiles,project.ntraces,project.nmasses);
    
    thesescans=[1:project.pw*4*project.ntestfiles];
    usemasses = project.deco{block}.usemass;
    sopt = project.deco{block}.sopt(:,usemasses);
    [a(thesescans,:)]= deco_als99(block,totblock(thesescans,usemasses),sopt,project.ntestfiles,1); 
    npeaks = size(a,2); % number of peaks
    areaopt= zeros(npeaks,project.ntestfiles); % area of peak in test spectrum
    pmax   = zeros(npeaks,project.ntestfiles); % position of peak maximum
    %calculate the area of the predicted components
    for peak=1:npeaks,
        for file=1:project.ntestfiles
            areaopt(peak,file)= sum(a((file-1)*4*pw+1:file*4*pw,peak)) * s(peak);
            [hgt,pos] = max(a((file-1)*4*pw+1:file*4*pw,peak));
            pmax(peak,file)   = pos;
        end
    end
    project.deco{block}.tprofile = a; % store the profiles of the test data
    project.deco{block}.pred = areaopt; % store the area's 
    project.deco{block}.tmax = pmax; % store the location of the maxima
end
y=1;


