function Stream_Channel()

%Enter Server Name of NetCom Server.    Class: String
Server_Name = 'localhost';

%Enter CSC Name.                        Class: String.
Channel_Name = 'CSC1';

%Enter Plot Height in microvolts.       Class: Double
Plot_Heigth = 50000;

%Enter Plot Length in seconds.          Class: Double
Plot_Length = 30;

%Enter Update Time in seconds.          Class: Double
Update_Time = .5;


%Connect to Server_Name using NlxConnectToServer
NlxDisconnectFromServer();
succeeded = NlxConnectToServer(Server_Name);

if succeeded == 1
    %Connection succeeded
    disp(['Connection to ' Server_Name ' established.']);
    
    %Get the sample frequency using NlxSendCommand
    [~, cheetahReply] = NlxSendCommand(['-GetSampleFrequency ' Channel_Name]);
    
    %Preallocate Channel_Data_Array
    Channel_Data_Array = zeros(1, (Plot_Length * str2double(cheetahReply{1})));
    Channel_Time_Array = 1/str2double(cheetahReply{1}):1/str2double(cheetahReply{1}):Plot_Length;
    
    %Open a stream to Channel_Name using NlxOpenStream
    NlxOpenStream(Channel_Name);    
    
    %Configure the plot
    set(gcf,'CurrentCharacter','a')
    plot(Channel_Time_Array,Channel_Data_Array);
    
    while(double(get(gcf,'CurrentCharacter'))~=27)
        %Pause Update_Time to allow the buffer for Channel_Name to fill.
        pause(Update_Time-0.2);
        
        %Get the sample frequency using NlxSendCommand
        [~, cheetahReply] = NlxSendCommand(['-GetInputRange ' Channel_Name]);
        
        %Get new data for Channel_Name using NlxGetNewCSCData
        [~,dataArray, ~, ~, ~, ~, ~, ~ ] = NlxGetNewCSCData(Channel_Name);

        %Convert dataArray to uV
        dataArray = double(dataArray)*str2double(cheetahReply{1})/32767;
        
        %Concatenate Channel_Data_Array and dataArray
        Channel_Data_Array(1,1:(length(Channel_Data_Array)-length(dataArray))) = Channel_Data_Array(1,(length(dataArray)+1):length(Channel_Data_Array));
        Channel_Data_Array(1,(length(Channel_Data_Array)-length(dataArray)+1):length(Channel_Data_Array)) = dataArray;
        
        %Plot the data
        %figure(1)
        plot(Channel_Time_Array,Channel_Data_Array);
        xlabel('Time [s]');
        ylabel('Microvolts [uV]');
        title(['Streaming ' Channel_Name '.  Press Esc to stop.']);   
        ylim([-Plot_Heigth Plot_Heigth]);
    end   
    
    %Close the stream to Channel_Name using NlxCloseStream
    NlxCloseStream(Channel_Name);    
    
    %Update Title
    title(['Stream for ' Channel_Name ' closed.']); 
  
else
    %Connection failed.
    disp(['Unable to connect to ' Server_Name '.']);
end

end
