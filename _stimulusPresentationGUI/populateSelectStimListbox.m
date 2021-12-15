function stimuli = populateSelectStimListbox(handles)
global pm
contents = cellstr(get(handles.listbox5,'String'));
projectSel = contents{get(handles.listbox5,'Value')};
pm.wavDir = [pm.stimFolder pm.person '\' projectSel '\'];
% POPULATE SELECT STIMULI LISTBOX
pm.wavFiles = dir([pm.wavDir '*.wav']);
pm.wavFiles = {pm.wavFiles.name};
set(handles.listbox2,'Value',1);
set(handles.listbox2,'String',pm.wavFiles);
set(handles.listbox2,'Max',length(pm.wavFiles),'Min',1);
stimuli = pm.wavFiles;