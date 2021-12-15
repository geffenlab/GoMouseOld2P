function editBlockNumbers(wavs,d,handles)

if ~iscell(d)
    d = num2cell(d);
end
input = cat(2,wavs,d);
f = figure;
set(f,'Name','Edit block numbers then hit ESCAPE or Done to return to GoMouse!')
uitable(f,'Data',input,'ColumnWidth',{360,40,60},'Position',[20 50 500 350],...
    'ColumnName',{'Stimulus','Block','Repeats'},'RowName',{},'ColumnEditable',logical([0,1,1]),...
    'KeyPressFcn',@(src,event)getDataAndClose(src,event,f,handles));
uicontrol('Style', 'pushbutton', 'String', 'Done',...
        'Position', [470 20 50 20],...
        'Callback', @(src,event)gdac(src,event,f,handles));  


function getDataAndClose(src,event,f,handles)
if strcmp('escape',event.Key)
    output = f.Children(2).Data(:,2:3);
    set(handles.uitable2,'Data',output);
    prepPresInfo(handles);
    close(f)
end

function gdac(src,event,f,handles)
    output = f.Children(2).Data(:,2:3);
    set(handles.uitable2,'Data',output);
    prepPresInfo(handles);
    close(f)



