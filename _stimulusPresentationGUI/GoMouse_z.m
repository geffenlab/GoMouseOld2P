function varargout = GoMouse_z(varargin)
% GOMOUSE_Z MATLAB code for GoMouse_z.fig
%      GOMOUSE_Z, by itself, creates a new GOMOUSE_Z or ra-`    ises the existing
%      singleton*.
%
%      H = GOMOUSE_Z returns the handle to a new GOMOUSE_Z or the handle to
%      the existing singleton*.
%                                                                                                                                                                         
%      GOMOUSE_Z('CALLBACK',hObjec  t,eventData,handles,...) calls the local
%      function named CALLBACK in GOMOUSE_Z.M with the given input arguments.
%
%      GOMOUSE_Z('Property','Value',...) creates a new GOMOUSE_Z or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GoMouse_z_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GoMouse_z_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GoMouse_z

% Last Modified by GUIDE v2.5 27-Aug-2018 14:26:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GoMouse_z_OpeningFcn, ...
    'gui_OutputFcn',  @GoMouse_z_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



%% SET VARIABLES
global pm
pm.stimFolder = 'C:\stimuli\';
pm.mouseFolder = 'C:\data\';
pm.filterFolder = 'C:\Users\MMMBop\Documents\GitHub\filters\';
addpath('C:\Users\MMMBop\Documents\GitHub\GoMouse_2P\_matlab_general\')



% --- Executes just before GoMouse_z is made visible.
function GoMouse_z_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GoMouse_z (see VARARGIN)

% Choose default command line output for GoMouse_z
handles.output = hObject;

% Update handles structure
global pm
guidata(hObject, handles);
set(hObject,'Name','GoMouse!')
pm.wavFolders = [];
set(handles.uitable2,'Data',[])

% UIWAIT makes GoMouse_z wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = GoMouse_z_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






%% SELECT YOUR BEAST
function listbox1_Callback(hObject, eventdata, handles) %#ok<*INUSD>
global pm
contents = cellstr(get(hObject,'String'));
pm.mouse = contents{get(hObject,'Value')};
pm.saveFolder = [pm.mouseFolder pm.mouse '\'];
if exist([pm.mouseFolder pm.mouse '\' pm.mouse 'pitchYaw.mat'],'file')
    load([pm.mouseFolder pm.mouse '\' pm.mouse 'pitchYaw.mat'])
    set(handles.edit2,'String',num2str(pitch));
    set(handles.edit1,'String',num2str(yaw));
else
    set(handles.edit2,'String',num2str(0));
    set(handles.edit1,'String',num2str(0));
end

% POPULATE SELECT YOUR BEAST
function listbox1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% Populate the listbox
global pm
f = dir(pm.mouseFolder);
f(strcmp({f.name},'.'))=[]; f(strcmp({f.name},'..'))=[]; f(strcmp({f.name},'zippedEventTraces'))=[];
f = {f.name};
set(hObject,'string',f);



%% SELECT FOLDER
function listbox4_Callback(hObject, eventdata, handles)
global pm
contents = cellstr(get(hObject,'String'));
pm.person = contents{get(hObject,'Value')};

% Populate the PROJECT listbox
populateProjectFolderListbox(handles)



% POPULATE SELECT FOLDER
function listbox4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% Populate the listbox
global pm
f = dir(pm.stimFolder);
f(strcmp({f.name},'.'))=[]; f(strcmp({f.name},'..'))=[];
f = {f.name};
set(hObject,'string',f);




%% SELECT PROJECT FOLDER
function listbox5_Callback(hObject, eventdata, handles)
populateSelectStimListbox(handles);


% POPULATE PROJECT FOLDER
function listbox5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% SELECT STIMULI LISTBOX
function listbox2_Callback(hObject, eventdata, handles)

% POPULATE SELECT STIMULI
function listbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% EDIT BLOCKS NUMBERS
function uitable2_ButtonDownFcn(hObject, eventdata, handles)
wavs = cellstr(get(handles.listbox3,'String'));
d = get(handles.uitable2,'Data');
editBlockNumbers(wavs,d,handles);


%% ADD PUSH BUTTON
function pushbutton2_Callback(hObject, eventdata, handles)
global pm
contents = cellstr(get(handles.listbox3,'String'));
wavSel = get(handles.listbox2,'Value');
wavsToAdd = cellstr(get(handles.listbox2,'String'));
wavsToAdd = wavsToAdd(wavSel);
if strcmp(contents,'Listbox')
    contents=wavsToAdd;
else
    contents(length(contents)+1:length(contents)+length(wavsToAdd)) = wavsToAdd;
end
%  pm.wavsToPresent
set(handles.listbox3,'String',contents)
set(handles.listbox3,'Max',length(contents),'Min',1);
for ii=1:length(wavsToAdd)
    pm.wavFolders{length(pm.wavFolders)+1,1} = pm.wavDir;
end
% edit block numbers
d = get(handles.uitable2,'Data');
a = get(handles.listbox3,'Value');
if ~iscell(d)
    d=num2cell(d);
end
if ~isempty(d)
    nd = cat(1,d,cat(2,repmat(d(a,1),length(wavsToAdd),1),num2cell(ones(length(wavsToAdd),1))));
else
    nd = ones(length(wavsToAdd),2);
end
set(handles.uitable2,'Data',nd)



%% REMOVE PUSH BUTTON
function pushbutton3_Callback(hObject, eventdata, handles)
global pm
contents = cellstr(get(handles.listbox3,'String'));
if length(contents)==1
    contents{1} = 'Listbox';
else
    rm = get(handles.listbox3,'Value');
    if rm(1)-1<1
        set(handles.listbox3,'Value',1);
    else
        set(handles.listbox3,'Value',rm(1)-1);
    end
    contents(rm)=[];
    if isempty(contents)
        contents{1} = 'Listbox';
    end
end
pm.wavFolders(rm)=[];
set(handles.listbox3,'String',contents)
set(handles.listbox3,'Max',length(contents),'Min',1);

d = get(handles.uitable2,'Data');
if ~isempty(d)
    d(rm,:)=[]; nd=d;
end
set(handles.uitable2,'Data',nd)




%% CLEAR ALL PUSH BUTTON
function pushbutton8_Callback(hObject, eventdata, handles)
global pm
set(handles.listbox3,'Value',1)
set(handles.listbox3,'String','Listbox')
pm.wavFolders=[];
set(handles.uitable2,'Data',[]);



%% STIMULI SELECTED LISTBOX
function listbox3_Callback(hObject, eventdata, handles)
% POPULATE STIMULI SELECTED
function listbox3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% MOVE UP PUSH BUTTON
function pushbutton4_Callback(hObject, eventdata, handles)
global pm
contents = cellstr(get(handles.listbox3,'String'));
mu = get(handles.listbox3,'Value');
d = get(handles.uitable2,'Data');
if mu~=1
    newContents = [contents(1:mu(1)-2);contents(mu);contents(mu(1)-1);contents(mu(end)+1:end)];
    pm.wavFolders = cat(1,pm.wavFolders(1:mu(1)-2),pm.wavFolders(mu),pm.wavFolders(mu(1)-1),pm.wavFolders(mu(end)+1:end));
    d = cat(1,d(1:mu(1)-2,:),d(mu,:),d(mu(1)-1,:),d(mu(end)+1:end,:));
    set(handles.listbox3,'String',newContents)
    set(handles.listbox3,'Value',mu-1)
    set(handles.uitable2,'Data',d)
end



%% MOVE DOWN PUSH BUTTON
function pushbutton5_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
global pm
contents = cellstr(get(handles.listbox3,'String'));
mu = get(handles.listbox3,'Value');
d = get(handles.uitable2,'Data');
if mu~=length(contents)
    newContents = [contents(1:mu(1)-1);contents(mu(end)+1);contents(mu);contents(mu(end)+2:end)];
    pm.wavFolders = cat(1,pm.wavFolders(1:mu(1)-1),pm.wavFolders(mu(end)+1),pm.wavFolders(mu),pm.wavFolders(mu(end)+2:end));
    d = cat(1,d(1:mu(1)-1,:),d(mu(end)+1,:),d(mu,:),d(mu(end)+2:end,:));
    set(handles.listbox3,'String',newContents)
    set(handles.listbox3,'Value',mu+1)
    set(handles.uitable2,'Data',d)
end


%% RANDOMISE STIMULI ORDER
function pushbutton12_Callback(hObject, eventdata, handles)
global pm
contents = cellstr(get(handles.listbox3,'String'));
d = get(handles.uitable2,'Data');
r = randperm(length(contents));
contents = contents(r);
pm.wavFolders = pm.wavFolders(r);
set(handles.listbox3,'String',contents)
set(handles.uitable2,'Data',d(r,:));







%% NOISE CLICKS TOGGLE BUTTON
function togglebutton1_Callback(hObject, eventdata, handles)
global nc
% clear -global nc
state = logical(get(hObject,'Value'));
if state
    noiseClickDur = 1; % seconds
    ISI = 5; % seconds
    fs = str2double(get(handles.edit9,'String'));
    filtName = get(handles.edit10,'String');
    load(filtName);
    noise = makeNoiseClick(noiseClickDur,ISI,fs,FILT); % duration, ISI and sample rate
    set(handles.text35,'String','Connecting to NIDAQ card');
    nc.s = connectToNidaq(fs,[],[0,1]);
    set(handles.text35,'String','NIDAQ connected');
    nc.lh = addlistener(nc.s,'DataRequired',@(src,event)continuousNoise(src,event,noiseClickDur,ISI,fs,FILT));
    nc.s.IsContinuous = true;
    nc.s.queueOutputData(noise);
    nc.s.startBackground();
    set(handles.text35,'String','Presenting noise clicks');
    disp('Presenting noise clicks')
else
    nc.s.stop();
    delete(nc.lh);
    nc.s.IsContinuous = false;
    clear -global nc
    set(handles.text35,'String','Stopped noise clicks, now nothing happening');
    disp('Stopped noise clicks')
end











%% BASELINE DURATION
function edit4_Callback(hObject, eventdata, handles)
%
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% FRAME RATE
function edit6_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% SAVE FILENAME
function edit7_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% ORBITAL YAW
function edit1_Callback(hObject, eventdata, handles)


function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% ORBITAL PITCH
function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% ZOOMS
% 1X ZOOM
function pushbutton13_Callback(hObject, eventdata, handles)
set(handles.edit6,'String',num2str(29.543))
set(handles.text34,'String','1X')
prepPresInfo(handles);
% 1.25X ZOOM
function pushbutton14_Callback(hObject, eventdata, handles)
set(handles.edit6,'String',num2str(29.874))
set(handles.text34,'String','1.25X')
prepPresInfo(handles);
% 2X ZOOM
function pushbutton18_Callback(hObject, eventdata, handles)
set(handles.edit6,'String',num2str(30.506))
set(handles.text34,'String','2X')
prepPresInfo(handles);
%
% 1.5X ZOOM
function pushbutton19_Callback(hObject, eventdata, handles)
set(handles.edit6,'String',num2str(30.395))
set(handles.text34,'String','1.5X')
prepPresInfo(handles);
%



%% RESET NIDAQ
function pushbutton9_Callback(hObject, eventdata, handles)
daqreset
delete(instrfindall)



%% UPDATE RECORDING INFO
function [presInfo] = pushbutton10_Callback(hObject, eventdata, handles)
global pm
presInfo = prepPresInfo(handles);
contents = cellstr(get(handles.listbox5,'String'));
projectSel = contents{get(handles.listbox5,'Value')};
set(handles.edit7,'String',[datestr(now,'yymmdd_HHMM') '_' pm.mouse '_' projectSel '_block01.txt'])



%% RECORD
function pushbutton11_Callback(hObject, eventdata, handles)
clear -global nc
global nc pm
nc.blockN = 1;
nc.mouse = pm.mouse;
nc.stimFolder = pm.stimFolder;
playNextBlock(handles)


% %% PRESENT SOUND ONLY
% function pushbutton23_Callback(hObject, eventdata, handles)
% clear -global nc
% global nc pm
% nc.playbackOnly = 1;
% nc.blockN = 1;
% nc.mouse = pm.mouse;
% nc.stimFolder = pm.stimFolder;
% 
% % Connect to the NIDAQ
% [~,chanOut] = getNidaqSettings(handles);
% nc.fs = str2double(get(handles.edit9,'String'));
% set(handles.text35,'String','Connecting to NIDAQ card');
% nc.s = connectToNidaq(nc.fs,[],chanOut);
% set(handles.text35,'String','NIDAQ connected');
% 
% % Add the listeners for continuous playback/acquisition
% nc.lh = addlistener(nc.s,'DataRequired',@(src,event)presentStimContNidaq_stimGUI(src,event,handles));
% nc.la = [];
% nc.s.IsContinuous = true; % set nidaq to continuous mode
% 
% % start counters
% nc.counter=1;
% 
% % Get info about what to present
% presInfo = prepPresInfo(handles);
% nc.firstChunk = 1;
% nc.nBlocks = presInfo.nBlocks;
% nc.nChunks = presInfo.nChunks(nc.blockN);
% nc.stimFiles = presInfo.stimFiles(presInfo.blocks==nc.blockN);
% nc.nFiles = length(nc.stimFiles);
% nc.stimDur = presInfo.stimDur{nc.blockN};
% nc.preStimSil = presInfo.preStimSil;
% 
% queueOutputData(nc.s,presInfo.triggerAcquisition);
% % Initialise the presentation/acquisition (the listeners take over after
% % triggerAcquisition has been presented
% nc.s.startBackground();
% set(handles.text35,'String',['Presenting block ' num2str(nc.blockN) ' of ' num2str(nc.nBlocks)])
% 




%% SAMPLE RATE
function edit9_Callback(hObject, eventdata, handles)

function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% FILTER SELECTION
function edit10_Callback(hObject, eventdata, handles)

function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global pm
set(hObject,'String',[pm.filterFolder '20211213_New2PStim_1k-70k_fs400k_AL.mat'])

function edit10_ButtonDownFcn(hObject, eventdata, handles)
global pm
[FileName,PathName] = uigetfile([pm.filterFolder '*.mat'],'Select filter');
if FileName~=0
    set(handles.edit10,'String',[PathName,FileName]);
end


%% WIDEFIELD RECORDING CHECKBOX
function checkbox1_Callback(hObject, eventdata, handles)
global saveData
a = get(hObject,'Value');
saveData.wfRec = a;
if a==1
    set(handles.text17,'String','Widefield imaging')
else
    set(handles.text17,'String','')
end

%% 2P RECORDING CHECKBOX
function checkbox2_Callback(hObject, eventdata, handles)
global saveData
a = get(hObject,'Value');
saveData.twoPhotonRec = a;
if a==1
    set(handles.text17,'String','2P imaging')
else
    set(handles.text17,'String','')
end

%% EPHYS RECORDING CHECKBOX
function checkbox3_Callback(hObject, eventdata, handles)
global saveData
a = get(hObject,'Value');
saveData.ephysRec = a;
if a==1
    set(handles.text17,'String','Electrophysiology recording')
else
    set(handles.text17,'String','')
end




%% NIDAQ SETTINGS
function edit11_Callback(hObject, eventdata, handles) % NUMBER OF INPUT CHANNELS
% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit12_Callback(hObject, eventdata, handles) % NUMBER OF OUTPUT CHANNELS
function edit12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit13_Callback(hObject, eventdata, handles) % 1ST INPUT CHANNEL
function edit13_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit14_Callback(hObject, eventdata, handles)% 1ST OUTPUT CHANNEL
function edit14_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit15_Callback(hObject, eventdata, handles) % 2ND INPUT CHANNEL
function edit15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit16_Callback(hObject, eventdata, handles) % 2ND OUTPUT CHANNEL
function edit16_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit17_Callback(hObject, eventdata, handles)% 3RD INPUT CHANNEL
function edit17_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit18_Callback(hObject, eventdata, handles)% 3RD OUTPUT CHANNEL
function edit18_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit19_Callback(hObject, eventdata, handles) % 4th INPUT CHANNEL
function edit19_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit20_Callback(hObject, eventdata, handles) % 4TH OUTPUT CHANNEL
function edit20_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% QUIT BUTTON
function pushbutton15_Callback(hObject, eventdata, handles)
fclose('all');
clear
clear global
close all


%% SAVE CONFIG
function pushbutton20_Callback(hObject, eventdata, handles)
global pm
c = cellstr(get(handles.listbox1,'String'));
GUIdata.mouse = c{get(handles.listbox1,'Value')}; % mouse
c = cellstr(get(handles.listbox4,'String'));
GUIdata.person = c{get(handles.listbox4,'Value')}; % person folder
c = cellstr(get(handles.listbox5,'String'));
GUIdata.project = c{get(handles.listbox5,'Value')}; % project folder
c = cellstr(get(handles.listbox3,'String'));
GUIdata.stimuli = c;
GUIdata.stimuliFolders = pm.wavFolders; % stimuli full filename
GUIdata.yaw = get(handles.edit2,'String'); % yaw
GUIdata.pitch = get(handles.edit1,'String'); % pitch
GUIdata.blockReps = get(handles.uitable2,'Data'); % blocks and repeats
[chanIn,chanOut,nIn,nOut] = getNidaqSettings(handles);
GUIdata.nidaqInfo = {chanIn,chanOut,nIn,nOut};
GUIdata.recType(1) = get(handles.checkbox1,'Value'); % widefield
GUIdata.recType(2) = get(handles.checkbox2,'Value'); % 2P
GUIdata.recType(3) = get(handles.checkbox3,'Value'); % electrophys
GUIdata.baselineDur = get(handles.edit4,'String');
GUIdata.opticalZoom = get(handles.text34,'String');
GUIdata.fs = get(handles.edit9,'String');
[filename,pathname] = uiputfile([pm.stimFolder GUIdata.person '\' GUIdata.project '\*.mat'],'save configuration');
if ~isempty(pathname)
    save([pathname filename],'GUIdata');
    set(handles.text35,'String','Configuration saved');
end





%% LOAD CONFIG
function pushbutton21_Callback(hObject, eventdata, handles)
global pm
[filename,pathname]=uigetfile(pm.stimFolder);
load([pathname filename]);
c = cellstr(get(handles.listbox1,'String'));
set(handles.listbox1,'Value',find(strcmp(c,GUIdata.mouse))); % mouse
pm.mouse = GUIdata.mouse;
pm.saveFolder = [pm.mouseFolder pm.mouse '\'];
c = cellstr(get(handles.listbox4,'String'));
set(handles.listbox4,'Value',find(strcmp(c,GUIdata.person))); % person folder
projects = populateProjectFolderListbox(handles);
set(handles.listbox5,'Value',find(strcmp(projects,GUIdata.project))); % project folder
populateSelectStimListbox(handles);
pm.wavFolders = GUIdata.stimuliFolders; % stimuli full filename
set(handles.listbox3,'String',GUIdata.stimuli);
set(handles.edit2,'String',GUIdata.yaw); % yaw
set(handles.edit1,'String',GUIdata.pitch); % pitch
set(handles.uitable2,'Data',GUIdata.blockReps); % blocks and repeats
set(handles.edit11,'String',GUIdata.nidaqInfo{3}); % set n channels in
chIn = [13,15,17];
for ii=1:GUIdata.nidaqInfo{3}
    eval(sprintf('set(handles.edit%d,''String'',%d);',chIn(ii),GUIdata.nidaqInfo{1}(ii)))
end
set(handles.edit12,'String',GUIdata.nidaqInfo{4}); % set n channels out
chOut = [14,16,18];
for ii=1:GUIdata.nidaqInfo{4}
    eval(sprintf('set(handles.edit%d,''String'',%d);',chOut(ii),GUIdata.nidaqInfo{2}(ii)))
end

set(handles.checkbox1,'Value',GUIdata.recType(1)); % widefield
set(handles.checkbox2,'Value',GUIdata.recType(2)); % 2P
set(handles.checkbox3,'Value',GUIdata.recType(3)); % electrophys

set(handles.edit4,'String',GUIdata.baselineDur);
set(handles.text34,'String',GUIdata.opticalZoom);
switch GUIdata.opticalZoom
    case '1X'
        pushbutton13_Callback(hObject, eventdata, handles)
    case '1.25X'
        pushbutton14_Callback(hObject, eventdata, handles)
end
set(handles.edit9,'String',GUIdata.fs);


%% ABORT BUTTON
function pushbutton22_Callback(hObject, eventdata, handles)
global nc
if ~isempty(nc)
    stop(nc.s);
    delete(nc.lh);
    delete(nc.la);
%     if nc.playbackOnly==0
%         fclose(nc.fid);
%         fclose('all');
%         % save everything
%         set(handles.text35,'String',['Saving block ' num2str(nc.blockN)])
%         exptInfo.mouse = nc.mouse;
%         exptInfo.stimFiles = nc.stimFiles;
%         b = unique(exptInfo.stimFiles);
%         for ii=1:length(b)
%             a = load([b{ii}(1:end-4) '_stimInfo.mat']);
%             exptInfo.stimInfo{ii} = a.stimInfo;
%         end
%         exptInfo.preStimSilence = nc.preStimSil;
%         exptInfo.fsStim = nc.fs;
%         exptInfo.yaw = nc.yaw;
%         exptInfo.pitch = nc.pitch;
%         fn = get(handles.edit7,'String');
%         save([fn(1:end-4) '_exptInfo.mat'],'exptInfo')
%     end
    clear -global nc
    fclose('all');
    set(handles.text35,'String','Presentation aborted, nothing happening now')
end










