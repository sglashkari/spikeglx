function succeeded = CheetahClose(cheetahObjects)
% CheetahClose(CheetahObjects)
% Close NetCom connection 

DispMess = true;

%close all open streams before disconnecting
Success = zeros(1,numel(cheetahObjects));
for index = 1:numel(cheetahObjects)
    Success(index) = NlxCloseStream(cheetahObjects(index));
    if ~Success(index)
        disp(sprintf('FAILED to close stream for %s', char(cheetahObjects(index))));
    end
end;
if all(Success) && DispMess
    disp 'PASSED close stream for all current objects'
else
    disp 'Failed to close one or more current objects'
end


%Disconnects from the server and shuts down NetCom
succeeded = NlxDisconnectFromServer();
if succeeded ~= 1 && DispMess
    disp 'FAILED disconnect from server'
elseif DispMess
    disp 'PASSED disconnect from server'
end


    