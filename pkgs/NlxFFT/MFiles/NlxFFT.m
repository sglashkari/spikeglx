function varargout = NlxFFT(varargin)
%GUI M-file for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('Property','Value',...) creates a new GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to NlxFFT_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI('CALLBACK') and GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 08-Aug-2013 14:13:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NlxFFT_OpeningFcn, ...
                   'gui_OutputFcn',  @NlxFFT_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before GUI is made visible.
function NlxFFT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;
UserData = get(handles.output, 'UserData');
%Setup timer
delete(timerfind('tag', 'NlxFFTIntervalT'));
UserData.TimerData.NlxFFTIntervalT = timer;
set(UserData.TimerData.NlxFFTIntervalT, ...
    'ErrorFcn', 'ErrorNlxFFTIntervalT(handles)', ...
    'executionmode', 'fixedRate', ...
    'Name', 'NlxFFTIntervalT', ...
    'Period', 1, ...
    'StartDelay', 0, ...
    'TimerFcn', 'CalcNlxFFT(handles);', ...
    'Tag', 'NlxFFTIntervalT');
% Update handles structure
set(handles.output, 'UserData', UserData);
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.GUI);


% --- Outputs from this function are returned to the command line.
function varargout = NlxFFT_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;



function DataFolder_Callback(hObject, eventdata, handles)
% hObject    handle to DataFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Folder = get(handles.DataFolder, 'string');
if ~strcmp(Folder(end), '\')
    set(handles.DataFolder, 'String', [Folder '\']);
end
    
% Hints: get(hObject,'String') returns contents of DataFolder as text
%        str2double(get(hObject,'String')) returns contents of DataFolder as a double


% --- Executes during object creation, after setting all properties.
function DataFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Connect.
function Connect_Callback(hObject, eventdata, handles)
% hObject    handle to Connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(timerfind('Name', 'NlxFFTIntervalT'));
set(handles.CalcFFT, 'value', 0);
UserData = get(handles.output, 'UserData');
set(handles.TimeInterval, 'enable', 'on');
% Disconnect and reconnect
success = NlxDisconnectFromServer;
if success == -1
GUIMessage(handles, 'Netcom can not be found');
    return; 
end
% All CSDData will be updated from NetComConnect
UserData.ConnectData = [];
set(handles.output, 'UserData', UserData);

% handles are updated in NetComConnect
NetComConnect(handles);


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to LoadInputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(timerfind('Name', 'NlxFFTIntervalT'));
set(handles.CalcFFT, 'value', 0);
UserData = get(handles.output, 'UserData');
%Browse
set(handles.DataFolder, 'string', [uigetdir '\']);
UserData.ConnectData.Path = get(handles.DataFolder, 'string');
Folder = get(handles.DataFolder, 'string');
set(handles.TimeInterval, 'enable', 'off');
% Offline is currently oneshot
GUIMessage(handles, 'Offline files are currently only single shot');
set(handles.OneShot, 'value', 1);
set(handles.TimeInterval, 'string', 0)

if NlxAreWeConnected
    GUIMessage(handles, ['Disconnecting from ' NlxGetServerPCName]);
    NlxDisconnectFromServer;
end
if ~exist(Folder, 'dir')
    GUIMessage(handles, 'Please select folder');
    return
end

UserData.ConnectData.CSCs = CscObj('path', Folder);
GUIMessage(handles, 'Files loaded');

[UserData.ConnectData.cheetahObjects, SortIndex] = sort_nat({UserData.ConnectData.CSCs.File});
UserData.ConnectData.CSCs = UserData.ConnectData.CSCs(SortIndex);
UserData.ConnectData.cheetahTypes = {UserData.ConnectData.CSCs.cheetahType};
UserData.ConnectData.Input = 'file';
set(handles.output, 'UserData', UserData);
UpdateCSCTable_Callback([ ], [ ], handles);



% --- Executes on button press in CommonZero.
function CommonZero_Callback(hObject, eventdata, handles)
% hObject    handle to CommonZero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(timerfind('Name', 'NlxFFTIntervalT'));
set(handles.CalcFFT, 'value', 0);


function TimeInterval_Callback(hObject, eventdata, handles)
% hObject    handle to TimeInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Interval = str2num(get(hObject, 'string'));
if isempty(Interval) || ~isreal(Interval) || Interval <0
    GUIMessage(handles, 'Interval must be >=0');
else
    UserData = get(handles.output, 'UserData'); 
    stop(UserData.TimerData.NlxFFTIntervalT);
    set(handles.CalcFFT, 'value', 0);
end
if Interval ~= 0
    set(handles.OneShot, 'value', 0);
else
    set(handles.OneShot, 'value', 1);
end
% Hints: get(hObject,'String') returns contents of TimeInterval as text
%        str2double(get(hObject,'String')) returns contents of TimeInterval as a double


% --- Executes during object creation, after setting all properties.
function TimeInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CalcFFT.
function CalcFFT_Callback(hObject, eventdata, handles)
% hObject    handle to CalcFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData = get(handles.output, 'UserData');
stop(UserData.TimerData.NlxFFTIntervalT);

%% Check for source
if ~isfield(UserData, 'ConnectData')
    GUIMessage(handles, 'No source specified, Connect or Load', 'blink', 3);
    return
end
%% Check for valid interval and Set Clock for single shot or fixed rate
Period = str2double(get(handles.TimeInterval, 'string'));
if isnan(Period)
    GUIMessage(handles, 'Set Interval >=0', 'blink', 3);
    return
else
    % Determine single shot or repetitive FFT
    if Period
        UserData.TimerData.NlxFFTIntervalT.ExecutionMode = 'fixedRate';
        UserData.TimerData.NlxFFTIntervalT.Period = Period;
    else
        UserData.TimerData.NlxFFTIntervalT.ExecutionMode = 'singleShot';
        UserData.TimerData.NlxFFTIntervalT.Period = 1;
    end
end
%% Check that CSC table is available
UserData.ConnectData.ChosenCSCs = [ ];
UserData.ConnectData.hButton = get(handles.hCSCSelGrp, 'children');
% handles are in reverse order, change to nat order
UserData.ConnectData.hButton = UserData.ConnectData.hButton(end:-1:1);
if isempty(UserData.ConnectData.hButton)
    GUIMessage(handles, 'Select CSCs, ''Press Update CSC Table'' if blank', 'blink', 3);
    return
end
%% Get Chosen CSCs
for csc = 1:numel(UserData.ConnectData.hButton)
    if get(UserData.ConnectData.hButton(csc), 'value')
        UserData.ConnectData.ChosenCSCs = [UserData.ConnectData.ChosenCSCs csc];
    end
end
%% Verify that at least one CSC is chosen
if ~any(UserData.ConnectData.ChosenCSCs)
    GUIMessage(handles, 'Select CSCs, ''Press Update CSC Table'' if blank', 'blink', 3); 
    return
end

%% Get Data per approptiate source and start timer
if get(hObject, 'value')
    UserData.FirstTime = true;
    if ~NlxAreWeConnected && strcmpi(UserData.ConnectData.Input, 'Server') && ...
            (isfieldRecursive(UserData, 'Path') && ~isempty(UserData.ConnectData.Path) || strcmpi(UserData.ConnectData.CSCs(1).cheetahType, 'CscFile'))
        GUIMessage(handles, 'Either Connect or Load', 'blink', 3);
    end
    if strcmpi(UserData.ConnectData.Input, 'file')
        %Source is file(s)
%%%%% This operation was done in Load, It is done here for connect because I don't have a listener set up        
%         UserData.ConnectData.CSCs = CscObj('path', UserData.ConnectData.Path);
    else%Source is Connection
        UserData.ConnectData.CSCs = CscObj('objname', UserData.ConnectData.cheetahObjects, 'objType', UserData.ConnectData.cheetahTypes);
    end
    set(handles.output, 'UserData', UserData);
    start(UserData.TimerData.NlxFFTIntervalT);
else %Calc button is off stop calculations
    stop(UserData.TimerData.NlxFFTIntervalT);
end    



% --- Executes on button press in checkboxCSC.
function checkboxCSC_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCSC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCSC


% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in UpdateCSCTable.
function UpdateCSCTable_Callback(~, ~, handles)
% hObject    handle to UpdateCSCTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData = get(handles.output, 'UserData');
%Stop calculations and reset table
stop(UserData.TimerData.NlxFFTIntervalT);
set(handles.CalcFFT, 'value', 0);
if isfieldRecursive(UserData, 'cheetahObjects')
    if numel(UserData.ConnectData.cheetahObjects) <= 40
        M = 8;
    elseif numel(UserData.ConnectData.cheetahObjects) <= 80
        M = 16;
    else
        M = 24;
    end
    % Remove all but CSCAcqEnt, (logical addressing)
    UserData.ConnectData.cheetahObjects = UserData.ConnectData.cheetahObjects(strcmpi(UserData.ConnectData.cheetahTypes, 'CSCAcqEnt'));
    UserData.ConnectData.hButton = ...
        SetRadioGUI(handles.hCSCSelGrp, M, 5, UserData.ConnectData.cheetahObjects);
else
    GUIMessage(handles, 'Objects are not available. Use Connect or Load to specify source.', 'fontcolor', 'red', 'blink', 3);
end

set(handles.output, 'UserData', UserData);


% --- Executes on button press in OneShot.
function OneShot_Callback(hObject, eventdata, handles)
% hObject    handle to OneShot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OneShot
stop(timerfind('Name', 'NlxFFTIntervalT'));
set(handles.CalcFFT, 'value', 0);
if get(hObject, 'value')
    set(handles.TimeInterval, 'string', 0)
end


% --- Executes on button press in dB.
function dB_Callback(hObject, eventdata, handles)
% hObject    handle to dB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dB


% --- Executes when user attempts to close GUI.
function GUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
diary off
% Hint: delete(hObject) closes the figure
delete(hObject);
