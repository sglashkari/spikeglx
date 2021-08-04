function Stream_Channel_Events_LuTom()

%Enter Server Name of NetCom Server.    Class: String
Server_Name = 'localhost';

%Enter Event Name.                        Class: String.
Channel_Name = 'Events';

%Enter Update Time in seconds.          Class: Double
Update_Time = .5;

NumStim=1;
Duration=0.00001;
StayTime=0;
step=0.1;
% Delay=0;
%Connect to Server_Name using NlxConnectToServer
NlxDisconnectFromServer();
succeeded = NlxConnectToServer(Server_Name);

if succeeded == 1
    %Connection succeeded
    disp(['Connection to ' Server_Name ' established.']);
    StateV=1;
    [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 0']);
    while succeeded ==1 
    %Open a stream to Channel_Name using NlxOpenStream
    NlxOpenStream(Channel_Name);    
    
%     [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 0']) % Setting Port 0 as Output back once the subject visited

        %Pause Update_Time to allow the buffer for Channel_Name to fill.
        pause(Update_Time-0.2);
        
              
        %Get new data for Channel_Name using NlxGetNewEventData
        [succeeded, timeStampArray, eventIDArray, ttlValueArray, eventStringArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewEventData('Events');
        %Last_Zone = char(eventStringArray(14:19));
        if eventIDArray ==6
           Zone = char(eventStringArray{end});
        end
         if numRecordsReturned > 0          
            Present_Zone = str2num(Zone(1,18));
            
            
                        switch StateV
                case 1
                     if Present_Zone==1
                         StimDone=0;
                         while StimDone<StayTime
                                 [succeeded, timeStampArray, eventIDArray, ttlValueArray, eventStringArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewEventData('Events');
                                 
                                 if eventIDArray ==6
                                    
                                    StateV=0;

                                 else
                                     pause(step);
                                     StimDone=StimDone+step;
                                     
                                 end
                                 
                         end
                         

%                     [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 1']) % Settining Port 0 as Input will prevent TTL to go out when subject visit same area in two successive entries.
%                     [succeeded, cheetahReply] = NlxSendCommand(['-SetDelay 200']) % Setting Port 0 as Output back once the subject visited
%                     [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 0']) % Setting Port 0 as Output back once the subject visited
                           if StateV==0                     
                             StateV=1;
                           else
                             BurstStim(NumStim,Duration);
                             StateV=2;
                           end
                     end

                case 2

                     if Present_Zone==2
                         StimDone=0;
                         while StimDone<StayTime
                                 [succeeded, timeStampArray, eventIDArray, ttlValueArray, eventStringArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewEventData('Events');
                                 
                                 if eventIDArray ==6
                                    
                                    StateV=0;

                                 else
                                     pause(step);
                                     StimDone=StimDone+step;
                                     
                                 end
                                 
                         end
                         

%                     [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 1']) % Settining Port 0 as Input will prevent TTL to go out when subject visit same area in two successive entries.
%                     [succeeded, cheetahReply] = NlxSendCommand(['-SetDelay 200']) % Setting Port 0 as Output back once the subject visited
%                     [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 0']) % Setting Port 0 as Output back once the subject visited
                           if StateV==0                     
                             StateV=2;
                           else
                             BurstStim(NumStim,Duration);
                             StateV=1;
                           end
                     end


            end

            
            
         end
%         if numRecordsReturned > 0
%            
%             Present_Zone = Zone(1,[14:18])
%             findstr('')
    end   
    %LIST;
    %Close the stream to Channel_Name using NlxCloseStream
    NlxCloseStream(Channel_Name);    
    
    %Update Title
    %title(['Stream for ' Channel_Name ' closed.']); 
  
else
    %Connection failed.
    disp(['Unable to connect to ' Server_Name '.']);
end
 
%end


function BurstStim(NumStim,Duration)

for i=1:NumStim

    [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 1']); % Settining Port 0 as Input will prevent TTL to go out when subject visit same area in two successive entries.
    [succeeded, cheetahReply] = NlxSendCommand(['-SetDelay 7']) % Setting Port 0 as Output back once the subject visited
    [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 0']) % Setting Port 0 as Output back once the subject visited

%     pause(Duration);
%     [succeeded, cheetahReply] = NlxSendCommand(['-SetDigitalIOPortValue AcqSystem1_0 0 0']); % Setting Port 0 as Output back once the subject visited
%      pause(Delay);

end

