function Data = CheetahRead(cheetahObjects, varargin)
if isempty(varargin)
    cheetahTypes ={'CscAcqEnt'};
else
    cheetahTypes = varargin{1};
end
DispMess = false;
% Data(objectIndex) = CheetahRead(cheetahObjects)
for objectIndex = 1:length(cheetahObjects)
    ObjectToRetrieve = char(cheetahObjects(objectIndex));
    %determine the type of acquisition entity we are currently indexed
    %to and call the appropriate function for that type
    if strcmp('CscAcqEnt', char(cheetahTypes(objectIndex))) == 1
%         Data = CSC;
        [succeeded, Data(objectIndex).DataArray, Data(objectIndex).TimeStampArray, Data(objectIndex).ChannelNumberArray, ...
            Data(objectIndex).SamplingFreqArray, Data(objectIndex).NumValidSamplesArray, Data(objectIndex).NumRecordsReturned, ...
            Data(objectIndex).NumRecordsDropped ] = NlxGetNewCSCData(ObjectToRetrieve);
        Data(objectIndex).cheetahType = cheetahTypes(objectIndex);
            
        if succeeded == 0 && DispMess
             disp(sprintf('FAILED to get new Data(objectIndex) for CSC stream %s on pass %d', ObjectToRetrieve, pass));
            return;
        elseif DispMess
            disp(sprintf('Retrieved %d CSC records for %s with %d dropped.', ...
                Data(objectIndex).NumRecordsReturned, ObjectToRetrieve, Data(objectIndex).NumRecordsDropped));

            %Here is where you'll perform some calculation on any of 
            %the returned values. Make sure any calculations done here
            %don't take too much time, otherwise NetCom will back up 
            %and you'' have dropped records
            plot(Data(objectIndex).DataArray);
        end
    %The test and actions are repeated for each acquisition entity type
    elseif strcmp('SEScAcqEnt', char(cheetahTypes(objectIndex))) == 1
        [succeeded, Data(objectIndex).DataArray, Data(objectIndex).TimeStampArray, Data(objectIndex).SpikeChannelNumberArray, ...
            Data(objectIndex).CellNumberArray, Data(objectIndex).FeatureArray, Data(objectIndex).NumRecordsReturned, ...
            Data(objectIndex).NumRecordsDropped ] = NlxGetNewSEData(ObjectToRetrieve);
        Data(objectIndex).cheetahType = cheetahTypes(objectIndex);
        if succeeded == 0 && DispMess
             disp(sprintf('FAILED to get new Data(objectIndex) for SE stream %s on pass %d', ObjectToRetrieve, pass));
            return;
        elseif DispMess
            disp(sprintf('Retrieved %d SE records for %s with %d dropped.', ...
                Data(objectIndex).NumRecordsReturned, ObjectToRetrieve, Data(objectIndex).NumRecordsDropped));
                Data(objectIndex).cheetahType = 'SES';
            plot(Data(objectIndex).DataArray);
        end
    elseif strcmp('STScAcqEnt', char(cheetahTypes(objectIndex))) == 1
        [succeeded, Data(objectIndex).DataArray, Data(objectIndex).TimeStampArray, Data(objectIndex).SpikeChannelNumberArray, ...
            Data(objectIndex).CellNumberArray, Data(objectIndex).FeatureArray, Data(objectIndex).NumRecordsReturned, ...
            Data(objectIndex).NumRecordsDropped ] = NlxGetNewSEData(ObjectToRetrieve);
        Data(objectIndex).cheetahType = cheetahTypes(objectIndex);
        if succeeded == 0 && DispMess
             disp(sprintf('FAILED to get new Data(objectIndex) for ST stream %s on pass %d', ObjectToRetrieve, pass));
            return;
        elseif DispMess 
            disp(sprintf('Retrieved %d ST records for %s with %d dropped.', ...
                Data(objectIndex).NumRecordsReturned, ObjectToRetrieve, Data(objectIndex).NumRecordsDropped));
            plot(Data(objectIndex).DataArray);
        end
     elseif strcmp('TTScAcqEnt', char(cheetahTypes(objectIndex))) == 1
        [succeeded, Data(objectIndex).DataArray, Data(objectIndex).TimeStampArray, Data(objectIndex).SpikeChannelNumberArray, ...
            Data(objectIndex).CellNumberArray, Data(objectIndex).FeatureArray, Data(objectIndex).NumRecordsReturned, ...
            Data(objectIndex).NumRecordsDropped ] = NlxGetNewSEData(ObjectToRetrieve);
        Data(objectIndex).cheetahType = cheetahTypes(objectIndex);
        if succeeded == 0 && DispMess
             disp(sprintf('FAILED to get new Data(objectIndex) for TT stream %s on pass %d', ObjectToRetrieve, pass));
            return;
        elseif DispMess
            disp(sprintf('Retrieved %d TT records for %s with %d dropped.', ...
                Data(objectIndex).NumRecordsReturned, ObjectToRetrieve, Data(objectIndex).NumRecordsDropped));
            plot(Data(objectIndex).DataArray);
        end
     elseif strcmp('EventAcqEnt', char(cheetahTypes(objectIndex))) == 1
        [succeeded, Data(objectIndex).TimeStampArray, Data(objectIndex).EventIDArray, Data(objectIndex).TTLValueArray, ...
            Data(objectIndex).EventStringArray, Data(objectIndex).NumRecordsReturned, Data(objectIndex).NumRecordsDropped ] = ...
            NlxGetNewEventData(ObjectToRetrieve);
        Data(objectIndex).cheetahType = cheetahTypes(objectIndex);
        if succeeded == 0 && DispMess
             disp(sprintf('FAILED to get new Data for event stream %s on pass %d', ObjectToRetrieve, pass));
            return;
        elseif DispMess
            disp(sprintf('Retrieved %d event records for %s with %d dropped.', ...
                Data(objectIndex).NumRecordsReturned, ObjectToRetrieve, Data(objectIndex).NumRecordsDropped));
            for recordIndex=1:Data(objectIndex).NumRecordsReturned                
                disp(sprintf('Event String: %s Event ID: %d TTL Value: %d', ...
                    char(Data(objectIndex).EventStringArray(recordIndex)), Data(objectIndex).EventIDArray(recordIndex), ...
                    Data(objectIndex).TTLValueArray(recordIndex)));
            end
        end
     elseif strcmp('VTAcqEnt', char(cheetahTypes(objectIndex))) == 1
        [succeeded,  Data(objectIndex).TimeStampArray, Data(objectIndex).ExtractedLocationArray, Data(objectIndex).ExtractedAngleArray, ...
            Data(objectIndex).NumRecordsReturned, Data(objectIndex).NumRecordsDropped ] = NlxGetNewVTData(ObjectToRetrieve);
        Data(objectIndex).cheetahType = cheetahTypes(objectIndex);
        if succeeded == 0 && DispMess
             disp(sprintf('FAILED to get new Data for VT stream %s on pass %d', ...
                 ObjectToRetrieve, pass));
            return;
        elseif DispMess
            disp(sprintf('Retrieved %d VT records for %s with %d dropped.', ...
                Data(objectIndex).NumRecordsReturned, ObjectToRetrieve, Data(objectIndex).NumRecordsDropped));
           plot(Data(objectIndex).ExtractedLocationArray);
        end   
    end
end

if succeeded == 0
    return;
end

if succeeded == 0 && DispMess
    disp 'FAILED to get Data consistently for all open streams'
elseif DispMess
    disp 'PASSED get Data consistently for all open streams'
end

