function [chanIn,chanOut,nIn,nOut] = getNidaqSettings(handles) %#ok<STOUT>

chIn = [13,15,17,19];
nIn = str2double(get(handles.edit11,'String'));
for ii=1:nIn
    eval(sprintf('chanIn(%d) = str2double(get(handles.edit%d,''String''));',ii,chIn(ii)))
end
if nIn~=sum(~isnan(chanIn))
    warndlg('Your NIDAQ input channels are not right, check number of channels and channel definitions in the NIDAQ settings panel','You need to fix something')
end
chOut = [14,16,18,20];
nOut = str2double(get(handles.edit12,'String'));
for ii=1:nOut
    eval(sprintf('chanOut(%d) = str2double(get(handles.edit%d,''String''));',ii,chOut(ii)))
end
if nOut~=sum(~isnan(chanOut))
     warndlg('Your NIDAQ output channels are not right, check number of channels and channel definitions in the NIDAQ settings panel','You need to fix something')
end

