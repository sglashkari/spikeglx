function varargout = NetComConnect(varargin)
%   NetComConnect is a GUI that connects to a computer. It also gets a
%   configuration file for Cheetah commands
% Output: Stored in object handles.output 'UserData' as structure ConnectData
%   .cheetahObjects
%   .cheetahTypes
%   .status          - Connection status
%   .Server          - IP or computer name
%   .ConfigFile      - Cheetah Command file

% If NetComConnect is called by another GUI, This GUI is called the Parent.
% All connection data is stored in (handles.output, 'UserData').
% If called by a parent then connection data is in
% (handles.output, 'UserData')
% (ConnectData.ParentsHandles.output, 'UserData'), otherwise
% Connection data is stored in (handles.output, 'UserData'), and is only good for
% the life of NetComConnect

% NETCOMCONNECT M-file for NetComConnect.fig
%      NETCOMCONNECT, by itself, creates a new NETCOMCONNECT or raises the existing
%      singleton*.
%
%      H = NETCOMCONNECT returns the handle to a new NETCOMCONNECT or the handle to
%      the existing singleton*.
%
%      NETCOMCONNECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NETCOMCONNECT.M with the given input arguments.
%
%      NETCOMCONNECT('Property','Value',...) creates a new NETCOMCONNECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NetComConnect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NetComConnect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIhandles

% Edit the above text to modify the response to help NetComConnect

% Last Modified by GUIDE v2.5 17-Sep-2012 11:19:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NetComConnect_OpeningFcn, ...
                   'gui_OutputFcn',  @NetComConnect_OutputFcn, ...
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


% --- Executes just before NetComConnect is made visible.
function CD = NetComConnect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NetComConnect (see VARARGIN)


% Choose default command line output for NetComConnect
handles.output = hObject;
% set(handles.output,'CloseRequestFcn', 'CloseNetComConnect'); % use this
% if you have an external shut down program

%% Opening Sequence
if ~isempty(varargin)
    %get handles structure and data from (handles.output, 'userData')
    % Assumes that UserData is avail if not then do nothing
    % This will allow to retrieve all the data from this GUI when 
    % returning to calling function, (parent GUI). ConnectData will be updated
    % in the parent's handle structure when 'connect' is pressed and the GUI is deleted.pdated  
    ParentsHandles = varargin{1};
    ParentsUserData = get(ParentsHandles.output, 'UserData');
    ConnectData = ParentsUserData.ConnectData;
    ConnectData.ParentsHandles = ParentsHandles;
else %Create a holder for field that is not a handle
    ConnectData.ParentsHandles.output = {'None'};
end
ConnectData.HasParent = ishandle(ConnectData.ParentsHandles.output);

%% Callback body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get Fresh connection, populate fields with default data from ParentGUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Verify Connection Data is valid
if ConnectData.HasParent 
    % Get Config from Parent's Data
    if isfieldRecursive(ConnectData, 'ConfigFile')
        set(handles.PathFile, 'String', ConnectData.ConfigFile);
    else
        ConnectData.ConfigFile = [ ];
    end
    % Create ServerName, IP if ~ exist
    if ~isfieldRecursive(ConnectData, 'ServerName')
        ConnectData.ServerName = [ ];
    end
    if ~isfieldRecursive(ConnectData, 'ServerIP')
        ConnectData.ServerIP = [ ];
    end
    % Set Name or IP on Connect GUI
    if isfieldRecursive(ConnectData, 'NameT_OR_IPF')
        if ConnectData.NameT_OR_IPF
            set(handles.RemoteComputer, 'String', ConnectData.ServerName);
        else
            set(handles.RemoteComputer, 'String', ConnectData.ServerIP);
        end
    end
end
        
% Disconnect, reconnect for fresh connection
if NlxAreWeConnected
    if NlxDisconnectFromServer;
        GUIMessage(handles, 'Disconnected');
    else%Did not succeed in disconnect
        GUIMessage(handles, 'Unable to Disconnect for fresh connection', 'ForegroundColor', 'Red', 'blink',2);
    end
end

%% Closing Sequence
if ConnectData.HasParent
    ParentsUserData.ConnectData = ConnectData;
    set(ConnectData.ParentsHandles.output, 'UserData', ParentsUserData);
end
set(handles.output, 'UserData', ConnectData);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NetComConnect wait for user response (see UIRESUME)
% uiwait(handles.NetComConnect);


% --- Outputs from this function are returned to the command line.
function varargout = NetComConnect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;



% --- Executes during object creation, after setting all properties.
function RemoteComputer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RemoteComputer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RemoteComputer_Callback(hObject, eventdata, handles)
% hObject    handle to RemoteComputer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RemoteComputer as text
%        str2double(get(hObject,'String')) returns contents of RemoteComputer as a double

% --- Executes on button press in RemoteComputerConnect.
function RemoteComputerConnect_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to RemoteComputerConnect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Opening Sequence
ConnectData = get(handles.output, 'UserData');
ConnectData.HasParent = ishandle(ConnectData.ParentsHandles.output);
if ConnectData.HasParent
    ParentsUserData = get(ConnectData.ParentsHandles.output, 'UserData');
    ConnectData = ParentsUserData.ConnectData;
end


%% Callback Body
% All the data on this GUI will be sent to the handles structure of the
% calling GUI and then this GUI will be deleted

DispMessages = true;
%Attempt Connection
NlxDisconnectFromServer;
GUIMessage(handles, ['Attempting to connect to ' get(handles.RemoteComputer, 'String')], 'pause', 1);
if NlxConnectToServer(get(handles.RemoteComputer, 'String'))
    ConnectData.ServerName = NlxGetServerPCName;
    ConnectData.ServerIP = NlxGetServerIPAddress;
    ConnectData.Input = 'Server';
    [~, Fs] = NlxSendCommand('-GetSampleFrequency CSC1');
    ConnectData.Fs = str2double(Fs{1});
    [S ConnectData.IOBoards] = NlxSendCommand('-GetDigitalIOBoardList');
    GUIMessage(handles, ['Connected to ' ConnectData.ServerName ' (' ConnectData.ServerIP ')' ]);
    if ConnectData.HasParent
        GUIMessage(ConnectData.ParentsHandles, ['Connected to ' ConnectData.ServerName ' (' ConnectData.ServerIP ')' ]);
        pause(2)
    end
else
    GUIMessage(handles, ['Unable to connect to ' get(handles.RemoteComputer, 'String') ], 'blink',2);
    return
end

% Send Cheetah config commands if file exists
ConfigFile = get(handles.PathFile, 'String');
if exist(ConfigFile, 'file') == 2;
    ConnectData.ConfigFile = ConfigFile;
    [ConnectData.NlxMatCommands ConnectData.ReturnMessages ConnectData.ReturnMessagesCounter] = ...
        NlxSetupCheetah(ConfigFile);
else
    GUIMessage(handles, 'Connecting without Config file');
end

%Get Cheetah objects and types
%%%%%%%%%%%%%%%%%%%%%% ARG check for ParentsHandles.GUITitle %%%%%%%%%%%%%%
if ConnectData.HasParent 
    try
        AppName  =  get(ConnectData.ParentsHandles.GUI, 'Name');
    catch
        AppName  =  get(ConnectData.ParentsHandles.figure1, 'Name');
    end
else
    AppName =  'Untitled';
end

[cheetahObjects, cheetahTypes, status] = CheetahConnect(ConnectData.ServerIP, AppName, DispMessages);
   
ConnectData.cheetahObjects = cheetahObjects;
ConnectData.cheetahTypes = cheetahTypes;
ConnectData.status = status;

% Update Message board of results
if status
   GUIMessage(handles, ['Connected to ' ConnectData.ServerName], 'pause', 1);
   ConnectData.Message = ['Connected to ' ConnectData.ServerName];
else
   GUIMessage(handles, 'Connection Not Successful', 'pause', 1);
   ConnectData.Message = 'Connection Not Successful';
end


%% Closing Sequence
if ConnectData.HasParent
    ParentsUserData.ConnectData = ConnectData;
    set(ConnectData.ParentsHandles.output, 'UserData', ParentsUserData);
end
set(handles.output, 'UserData', ConnectData);

%% Delete Connect GUI or minimize it
if NlxAreWeConnected
    pause(1)
    v = ver('Matlab');
    if str2double(v(1).Version) <= 8.0
        warning off
        %Warning: figure JavaFrame property will be obsoleted in a future release. For more information see the JavaFrame resource on the MathWorks web site. 
        jFrame = get(handles.NetComConnect, 'JavaFrame');
        jFrame.setMinimized(true)
        warning on
        % These are only necessary if you aren't going to delete NetComConnect GUI
        set(handles.output, 'Userdata', ConnectData);
        guidata(hObject, handles);
    else
        % Delete NetComConnect GUI, Data is saved in parent
        delete(get(hObject, 'parent'))
    end
end



% --- Executes on button press in RemoteComputerDisconnect.
function RemoteComputerDisconnect_Callback(hObject, eventdata, handles)
% hObject    handle to RemoteComputerDisconnect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Opening Sequence
ConnectData = get(handles.output, 'UserData');
ConnectData.HasParent = ishandle(ConnectData.ParentsHandles.output);
if ConnectData.HasParent
    ParentsUserData = get(ConnectData.ParentsHandles.output, 'UserData');
    ConnectData = ParentsUserData.ConnectData;
end


%% Callback Body
if NlxAreWeConnected;
    %Retrieve parents (Connection data)
    if ~isempty(ConnectData.cheetahObjects)
        succeeded = CheetahClose(ConnectData.cheetahObjects);
    end
    succeeded = NlxDisconnectFromServer();
    if succeeded || ~NlxAreWeConnected()
        set(handles.Message, 'String', 'Disconnected');
    end
    guidata(hObject, handles);
else
    GUIMessage(handles, 'Not Connected');
end

%% Closing Seqence
if ConnectData.HasParent
    ParentsUserData.ConnectData = ConnectData;
    set(ConnectData.ParentsHandles.output, 'UserData', ParentsUserData);
end
set(handles.output, 'UserData', ConnectData);
% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Opening Sequence
ConnectData = get(handles.output, 'UserData');
ConnectData.HasParent = ishandle(ConnectData.ParentsHandles.output);
if ConnectData.HasParent
    ParentsUserData = get(ConnectData.ParentsHandles.output, 'UserData');
    ConnectData = ParentsUserData.ConnectData;
end

%% Callback Body
[FileName,PathName,FilterIndex] = uigetfile('.txt');
if ischar(FileName) && ischar(PathName)
    set(handles.PathFile,'String',sprintf('%s%s',PathName,FileName));
    ConnectData.ConfigFile = sprintf('%s%s',PathName,FileName);
    set(handles.Message, 'string', '')
else
    set(handles.Message, 'string', 'Invalid File')
end 

%% Closing Sequence
if ConnectData.HasParent
    ParentsUserData.ConnectData = ConnectData;
    set(ConnectData.ParentsHandles.output, 'UserData', ParentsUserData);
end
set(handles.output, 'UserData', ConnectData);
guidata(hObject, handles);


% --- Executes on button press in VerifyConnection.
function VerifyConnection_Callback(hObject, eventdata, handles)
% hObject    handle to VerifyConnection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if NlxAreWeConnected()
    set(handles.Message, 'string', 'Valid Connection');
else
    set(handles.Message, 'string', 'Invalid Connection. Disconnect and then Reconnect');
end
guidata(hObject, handles)


% --- Executes when user attempts to close NetComConnect.
function NetComConnect_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to NetComConnect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Opening Sequence
ConnectData = get(handles.output, 'UserData');
ConnectData.HasParent = ishandle(ConnectData.ParentsHandles.output);
if ConnectData.HasParent
    ParentsUserData = get(ConnectData.ParentsHandles.output, 'UserData');
    ConnectData = ParentsUserData.ConnectData;
end

%% Callback body
button = questdlg('This will remove access to Cheetah Objects through the handle structure, You will still be connected to Cheetah. Continue','Close NetComConnect','Yes','No','No'); 
if strcmp(button, 'Yes')
%     RemoteComputerDisconnect_Callback(hObject, eventdata, handles)
    delete(hObject);
end

%% Closing Sequence





% --- Executes on key press with focus on RemoteComputer and none of its controls.
function RemoteComputer_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to RemoteComputer (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmpi(eventdata.Key, 'return')
    RemoteComputerConnect_Callback(handles.RemoteComputer, [], handles);
end




