function y = deco_dentify(block,nr)
global project;

%---------------------------------------------------------
% voor automatisch opstarten van NISTMS library search
%---------------------------------------------------------
% definieer het file AUTOIMP.MSD
%   zet hierin een standaard verwerkingsfile bv: filspec.fil
%defineer dan een file : filespec.fil
% waarin de namen van de te verwerken files staan
% bv  c:\\ontwikkel\\finalv18\\test.msp
% met of de toevoeging overwrite of append om deze aan bestand van nist toe
% te voegen 

% where to locate nist ?

if exist('c:\NISTdemo2\MSSEARCH\nistms$.exe','file') ~=0
    ppos =  project.deco{block}.rt(nr); % peak position in sec
    fid = fopen('compound1.msp','w');
      fprintf(fid,'Name: Rt=%f Block %d\n',ppos,block);
      [a,b] = find(project.deco{block}.sopt(nr,:) > 0.0001);
      n     = size(a,2);
      fprintf(fid,'Num peaks %d\n',n);
      fprintf(fid,'%d %f\n',[project.theMasses(b)+project.minmass-1 ; project.deco{block}.sopt(nr,b)]);
    fclose(fid);
    fid1 = fopen('c:\NISTdemo2\MSSEARCH\filspec.fil','w');
    fprintf(fid1,'%s\\compound1.msp overwrite\n',pwd);
    fclose(fid1);
    system('c:\NISTdemo2\mssearch\nistms$.exe');
else
    msgbox('could not find NIST directory','Error','error')
end

y=1;
