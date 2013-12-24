function deco_save_project()
global project;

if (strcmp(project.name,project.basename)==0) % overwrite
      save ([project.name '.prj'], 'project');
      return;
end;

while (strcmp(project.name,project.basename))
    [a,b] = deco_projectname; % collect answer and new projectname
    if (strcmpi(a,'OK') && strcmpi(b,project.basename)~=1)
        project.name = b;
       save ([project.name '.prj'], 'project');
        return;
    elseif strcmpi(a,'Cancel') 
        return;
    else 
        disp('choose another name')
    end
end


