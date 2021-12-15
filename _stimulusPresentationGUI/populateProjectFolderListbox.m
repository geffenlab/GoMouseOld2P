function projects = populateProjectFolderListbox(handles)
global pm
contents = cellstr(get(handles.listbox4,'String'));
pm.person = contents{get(handles.listbox4,'Value')};
% Populate the PROJECT listbox
projects  = dir([pm.stimFolder pm.person '\']);
projects(strcmp({projects.name},'.'))=[]; projects(strcmp({projects.name},'..'))=[];
projects = {projects.name};
set(handles.listbox5,'String',projects);