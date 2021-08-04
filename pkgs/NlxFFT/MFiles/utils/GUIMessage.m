function handles = GUIMessage(handles, String, varargin)
% handles = GUIWarning(handles, String, 'parmeter1', value1, 'parmeter2', value2 ...)
% Parmeters, Values:
%%%%% Defined GUIMessage.m properties
% 'FontColor', [R G B]  --- Warning: this will changs all messages to the
% new color. Recommend usage is as alert then return to black
% 'Blink', N
% 'Pause', N, pause for N seconds
% 
%%%%% Automatic: These are properties of a GUI Text object and are therefore
% passed directly to the object without checking 
% 'Fontweight',  'light'|'normal'|'demi'|'bold' 
% 'FontSize', n pixels
% 'FontName', legal font name
% see documenation for Matlab GUI for complete list
% Used alone
% persistent Messages 
% FontSize, FontWeight, FontColor
% 'fontsize', FontSize, 'fontweight', FontWeight,  'foregroundcolor', FontColor
NMessages = 5;

%% arg chk
if nargin > 2
    % verify that there are an even number of optional arguements
    if rem(length(varargin),2)
        error('There must be an even number of optional arguments for GUIWarning')
    end
    if ishandle(handles.Message)
        hMessage = handles.Message;
    elseif ishandle(handles)
        hMessage= handles;
    else
        error('First Arg to GUIMessage must be a handle')
    end
end
%% FontColor
% if one of the parmeternames is 'fontcolor' then change it to
% 'ForegroundColor:
% 1)   Find index of parameter fontcolor if it exists
%       a.  looking at every other arg ->varargin(1:2:end)
%       b.  See if 'fontcolor' is one of them -> ismember(lower( ---  ), 'fontcolor')
%       c.   ismember returns a vector index of the position, find the position of the
%            first '1' -> find( ---  ,1)
%       d.  Then multiply be 2 and subtract 1 to get the index of 'fontcolor' in varargin ->
%             FCIndex =  --- *2-1;
FCIndex = find(ismember(lower(varargin(1:2:end)), 'fontcolor'),1)*2-1;
% 2)    Now change that arg to the actual property name in a Matlab GUI
if ~isempty(FCIndex); varargin(FCIndex) = {'ForegroundColor'}; end;

%% Blink
% Find if Blink is set, set to blink and remove Blink parameter from
% varagin
% Find position of Blink parameter if it exists
PropertyIndex = find(ismember(lower(varargin(1:2:end)), 'blink'),1)*2-1;
if ~isempty(PropertyIndex)
    % retreive the number of times to blink
    NBlink = varargin{PropertyIndex+1};
    %%% Remove Blink parameters from varargin
    %Create index vector with 'Blink' and 'Blink_val' set to 1, all others are 0
    VararginRemoveSet = zeros(1,length(varargin));
    VararginRemoveSet([PropertyIndex PropertyIndex+1]) = 1;
    varargin = varargin(not(VararginRemoveSet));
else 
    NBlink = false;
end

%% Pause Setup
PropertyIndex = find(ismember(lower(varargin(1:2:end)), 'pause'),1)*2-1;
if ~isempty(PropertyIndex)
    % retreive the time to pause
    PauseTime = varargin{PropertyIndex+1};
    %%% Remove parameters from varargin
    %Create index vector to remove parameter set
    VararginRemoveSet = zeros(1,length(varargin));
    VararginRemoveSet([PropertyIndex PropertyIndex+1]) = 1;
    varargin = varargin(not(VararginRemoveSet));
else
    PauseTime = 0;
end
%% Message
Messages = get(handles.Message, 'string');
if isempty(Messages) || size(Messages,1) < NMessages
     Messages = cell(NMessages,1);
end

% Truncate messages to NMessages
set(handles.Message, 'string', Messages(1:NMessages));


% Display Set
if ~isempty(String) %if string is null then don't add message
    Messages = [{String}; Messages];
end
if NBlink
    BlinkTime = .25;
    BlinkColor = {'red' 'black'};
    for n = 1:NBlink
        VarArgIn = {varargin{:}  'foregroundcolor' BlinkColor{1}};
        set(handles.Message, 'String', Messages, VarArgIn{:});    
        pause(BlinkTime)
        VarArgIn = {varargin{:}  'ForegroundColor' BlinkColor{2}};
        set(handles.Message, 'String', Messages, VarArgIn{:});    
        pause(BlinkTime)
    end
else
    set(handles.Message, 'String', Messages, varargin{:} );
end

%% Pause
pause(PauseTime);
guidata(handles.Message, handles);

