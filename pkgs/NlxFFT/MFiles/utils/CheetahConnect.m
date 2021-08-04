function [cheetahObjects cheetahTypes Connected] = CheetahConnect(serverName, AppName, DispMess)
% [CheetahNames] = CheetahConnect(ConfigFileName, serverName)

% CheetahConnect connects to Cheetah on servername using the parameters in ConfigFileName.
% 1)Connects
% 2)Sets application name
% 3)OPens streams for each object
% ConfigFileName:     Sets the channels to collect CSC(1:n) and the input range.
%                     Sets the operations to be done
% Returns null if fail on connect, name or retrieving objects
if nargin < 2
    AppName = 'Matlab Routine';
end

if nargin < 1
    serverName = 'localhost';
end

if nargin < 3
    DispMess = false;
end

% Connect to Cheetah using Nlx NetCom
if NlxAreWeConnected() 
    Connected = 2; %already connected
elseif NlxConnectToServer(serverName);
    Connected = 1; %Connection successful
else 
    Connected = 0; %Connection not successful
    cheetahObjects = [ ];
    cheetahTypes = [ ];
    return
end

if DispMess
    if Connected 
        disp(sprintf('Connected to %s.', serverName));
    else 
        disp(sprintf('FAILED connect to %s. Exiting script.', serverName));
        return
    end
end

%Identify this program to the server we're connected to.
AppNameSucceeded = NlxSetApplicationName(AppName);
if DispMess
    if AppNameSucceeded
        disp 'PASSED set the application name';
    else
        disp 'FAILED set the application name';
    end
end

%get a list of all objects in Cheetah, along with their types.
[succeed, cheetahObjects, cheetahTypes] = NlxGetCheetahObjectsAndTypes;

%open up a stream for all objects
S = zeros(1,numel(cheetahObjects));
for index = 1:length(cheetahObjects)
    S(index) = NlxOpenStream(cheetahObjects(index));
    if ~succeed && DispMess
        disp(sprintf('FAILED to open stream for %s', char(cheetahObjects(index))));
    end
end
cheetahObjects = cheetahObjects(logical(S));
cheetahTypes = cheetahTypes(logical(S));
Connected = true;
if all(S)
    disp 'PASSED open stream for all current objects'
else
    disp 'Failed to open stream for one or more current objects.'
end
