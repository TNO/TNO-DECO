function name = deco_get_name(i)
global project;
name = lower([project.sdir '\\' project.files{i}]);
[a,b] = fileparts(name);
name = [a b '.dbc'];
