function  [c1, s1, sdopt, numiter] = deco_fixspectra(block, x, npeaks, actfiles, usemasses)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% S. Krishnan
% Determining component iteratively with appropriate seeding strategy
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
npeaks = 4;
[row, col] = size(x);
c = zeros(row, npeaks);
rnd = 0;
% 
% c = deco_peakpos(x, npeaks, c);
% 
% is1 = c\x;
% [c1, s1, sdopt, numiter] = deco_als990(block, x, is1, actfiles);

% rnd = 0; % Actual implementation
% rnd = 1; % Random initialisation
% rnd = 2; % MCR ALS: sequential % this has no meaning with respect to MCR
% rnd = 3; % MCR ALS: simultaneous % original implementation of the mcr
% rnd = 4; % Actual implementation simultaneous

% rnd = 5; % Random initialisation of ext.nipals simultaneous
% rnd = 6; % Random initialisation of mcr.als
% rnd = 7; % peak exact position of mcr.als

wd = size(x, 1)/actfiles;
x1 = x;
for i=1:npeaks;
    if (rnd == 0)
        %actual implentation
        [mx, idx] = max(sum(x1,2));
        filenameOR = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_original\extN_itr_pp_recon.fig';
        filenameC = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_original\Oc.fig';
        filenameS = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_original\s';
    elseif (rnd == 1)
        % random intialisation
        idx = randi(wd);
        x1d = sum(x1,2);
        mx = x1d(idx);
        filenameOR = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_random_init\extN_itr_rnd_recon.fig';
        filenameC = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_random_init\Oc.fig';
        filenameS = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_random_init\s';
    elseif (rnd == 2)
        %mcr als: sequential
        idx = i*wd/npeaks;
        x1d = sum(x1,2);
        mx = x1d(idx);
        filenameOR = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_sequential\MSqorgvsrecon.fig';
        filenameC = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_sequential\Oc.fig';
        filenameS = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_sequential\s';
    elseif(rnd == 3)
        % mcr als: simultaneous
        idx = wd/npeaks:wd/npeaks:wd;
        x1d = sum(x1,2);
        mx = x1d(idx);
        for j = 1:npeaks
            gs = deco_makegauss(1:row, idx(j), wd*1/15);
            c(:,j) =  mx(j)*gs./max(gs);% create initiation gaussian peaks
%             c(idx,j) = mx(j);
        end
        is1 = c\x;
        [c1, s1, sdopt, numiter] = deco_als990(block, x, is1, actfiles);
        filenameOR = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_simultaneous\MSiorgvsrecon.fig';
        filenameC = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_simultaneous\Oc.fig';
        filenameS = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_simultaneous\s';
        break;
    elseif(rnd == 4)
        % actual implementation simultaneous
        idx = [37 3 11 48];
        x1d = sum(x1,2);
        mx = x1d(idx);
        for j = 1:npeaks
            c(idx,j) = mx(j);
        end
        is1 = c\x;
        [c1, s1, sdopt, numiter] = deco_als990(block, x, is1, actfiles);
        filenameOR = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_simultaneous\extN_sim_pp_recon.fig';
        filenameC = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_simultaneous\Oc.fig';
        filenameS = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_simultaneous\s';
        break;
    elseif (rnd == 5)
        % Random initialisation of ext.nipals simultaneous
        idx = [randi(wd) randi(wd) randi(wd) randi(wd)];
        x1d = sum(x1,2);
        mx = x1d(idx);
        for j = 1:npeaks
            c(idx,j) = mx(j);
        end
        is1 = c\x;
        [c1, s1, sdopt, numiter] = deco_als990(block, x, is1, actfiles);
        filenameOR = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_random_init_simul\extN_sim_rnd_recon.fig';
        filenameC = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_random_init_simul\Oc.fig';
        filenameS = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_random_init_simul\s';
        break;
    elseif (rnd == 6)
        % Random initialisation of mcr.als
        idx = [randi(wd) randi(wd) randi(wd) randi(wd)];
        x1d = sum(x1,2);
        mx = x1d(idx);
        for j = 1:npeaks
            gs = deco_makegauss(1:row, idx(j), wd*1/15);
            c(:,j) =  mx(j)*gs./max(gs);% create initiation gaussian peaks
%             c(idx,j) = mx(j);
        end
        is1 = c\x;
        [c1, s1, sdopt, numiter] = deco_als990(block, x, is1, actfiles);
                filenameOR = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_random_init\mcr_sim_rnd_recon.fig';
                filenameC = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_random_init\Oc.fig';
                filenameS = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_random_init\s';
        break;
    elseif (rnd == 7)
        % peak position of mcr.als
        idx = [37 3 11 48];
        x1d = sum(x1,2);
        mx = x1d(idx);
        for j = 1:npeaks
            gs = deco_makegauss(1:row, idx(j), wd*1/20);
            c(:,j) =  mx(j)*gs./max(gs);% create initiation gaussian peaks
            %             c(idx,j) = mx(j);
        end
        is1 = c\x;
        [c1, s1, sdopt, numiter] = deco_als990(block, x, is1, actfiles);
                filenameOR = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_peak_pos\mcr_sim_pp_recon.fig';
                filenameC = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_peak_pos\Oc.fig';
                filenameS = 'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\mcr_als_peak_pos\s';
        break;
    else
        % nothing
    end
    c(idx,i) = mx;
    is1 = c\x;
    [c1, s1, sdopt, numiter] = deco_als990(block, x, is1, actfiles);
    x1 = x - c1*s1;
end

z = c1*s1;

figure;
hold on;
% plot(sum(z(1:wd, :), 2)./max(sum(z(1:wd, :), 2)));
% plot(sum(x(1:wd, :), 2)./max(sum(x(1:wd, :), 2)));
% saveas(hplot,'C:\Latex\Elsevier\extnipals\Figures\figs28thsept2012forRevisedpaper\extnipals_original\example.fig','fig');
hplot = plot(sum(z(1:wd, :), 2), '--r');
hplot = plot(sum(x(1:wd, :), 2));
diff = sum(z(1:wd, :), 2) - sum(x(1:wd, :), 2);
od = sum(x(1:wd, :), 2);
lfitp = sum(diff.*diff)*100/sum(od.*od);
xlabel('scan number');
ylabel('intensity');
title(['fit% = ' num2str(lfitp)]);
saveas(hplot,filenameOR,'fig');
close;
hold off;

% 
% extnipals_original
% extnipals_random_init
% mcr_als_sequential
% mcr_als_simultaneous
% extnipals_simultaneous

figure;
hplot = plot(c1(1:wd,:));
xlabel('scan number');
ylabel('intensity');
saveas(hplot(end),filenameC,'fig');
close;

for i = 1:npeaks;
    figure;
    s1(i, (s1(i, :)<0.025)) = 0;
    hplot = bar(usemasses, s1(i, :));
    xlabel('m/z');
    ylabel('abundance');
    filename = [filenameS int2str(i) '.fig'];
    saveas(hplot, filename, 'fig');
    close;
end

end


