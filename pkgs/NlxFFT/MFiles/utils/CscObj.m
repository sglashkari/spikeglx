classdef CscObj 
% Instantiation: 
% CSCs = CscObj('File', {file(s)} | 'Path', path | 'objname', {name(s)}), 
% CscObj - creates an empty object
% 'File' looks for specific file. If path is not included, (ex)
%   CscObj('file', ['T:\Data\Filename '.ncs']) then current ML directory is
%   assumed. A warning will result if the file is not found.
% 'Path' will load all .ncs files on the provided path. 
% 'objname' assumes a valid cheetah connection   
% All CscObj must be followed with Read(CSCobjname) - for files, or Stream(CSCobjname) - for online
%
% The CscObj is a container for Nlx CSC data. As such calculations on the CSCs 
% such as fft are done separately with other functions such as Myffthelp 
% 
    
%     had CscObj < handle. idea was to save memory for large data sets
%     since handle passes references rather than values. But if a method is
%     called and assigned to another variable, the origninal data is still
%     changed so the results are not expected
    %object of CSC or Cheetah streamed data
    %cell list of object names or files optional on Instantiation 
    %obj = CscObj('propname', 'value')
    %          propname = 'File', 'ObjName'
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Issues to resolve
%     should be able to read muliple path, files, streams
% should be able to repetitively read in the above

% Has
% NlxGetNewCSCData: CscAcqEnt
% % Include streaming for 
% NlxGetNewSEData: TTScAcqEnt SEScAcqEnt TTScAcqEnt
% NlxGetNewEventData: EventAcqEnt
% NlxGetNewVTData: VTAcqEnt

    properties
        DataArray 
        TimeStampArray 
        ChannelNumberArray 
        SamplingFrequency 
        InputRange
        NumValidSamplesArray 
        NumRecordsReturned 
        NumRecordsDropped
        Header
        Fs
        File
        Path
        Ext
        Name
        Units = 'ADC';
        ADBitVolts = 2^(-16); %+/- 32768/32767
        cheetahType = 'CscAcqEnt';
        ExtractMode = 1;
        ModeArray = [];
%         input % 'CSC', 'NRD', 'Server'
%         inputsource % Input- CSC: Folder, NRD: file, server: ipaddress
    end
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Instantiate
        function Obj = CscObj(varargin)
            % Obj = CscObj(Names), where Source is 
            %        'Name' cell list of CSC Names
            %        'File'     cell list of path|files
            if mod(numel(varargin), 2)
                error('input arguments must be (''property'', value) pairs');
             else
                Props = reshape(varargin, 2, [])';
            end
            %Test for simulataneous use of Path and File
            if strfind(lower([Props{:,1}]), 'file') & strfind(lower([Props{:,1}]), 'path')
                error('Use either ''file'' or ''path''. Use of both is ambiguous')
            end
            for prop_i = 1:size(Props, 1)
                switch lower(Props{prop_i,1})
                    case 'file'
                        if ischar(Props{prop_i,2})
                            %File specifed as single string
                            Files = Props(prop_i,2);
                        else %Files specified as string(s) in cell
                            Files = Props{prop_i,2};
                        end
                        for csc_i = 1 : numel(Files)
                            %Verify arg is a file
                            if strcmp(Files{csc_i}(end), '\') 
                                error('file was specified, but path was given')
                            else
                                %Verify file exists, if so get full file.ext
                                file = dir([Files{csc_i} '*']);
                                if ~isempty(file)
                                    [Obj(csc_i).Path, ~, ~] = fileparts(Files{csc_i});
                                    if isempty(Obj(csc_i).Path)
                                        %File was found without path,therefore it must be working
                                        %directory
                                        Obj(csc_i).Path = pwd;
                                    end
                                    [~, Obj(csc_i).File, Obj(csc_i).Ext] = fileparts(file.name);
                                    Obj(csc_i).Name = Obj(csc_i).File;
                                    %fileparts returns path without '\', which is official syn for path 
                                    Obj(csc_i).Path(end+1) = '\';
                                else %file is not on disk, 
                                    [Obj(csc_i).Path, Obj(csc_i).File, Obj(csc_i).Ext] = fileparts(Files{csc_i});
                                    Obj(csc_i).Path(end+1) = '\';
                                    warning([Obj(csc_i).Path Obj(csc_i).File Obj(csc_i).Ext ' Not found on Disk'])
                                end
                            end
                        end
                        Obj(csc_i).cheetahType = 'CSCAcqEnt';
%                         Obj(csc_i).cheetahType = 'CscFile';
                    case 'path'
                        %Verify arg is a directory
                        if ~strcmp(Props{prop_i, 2}(end), '\') 
                            error('Path was specified, but either file was given or path specified as cell')
                        else
                            Files = dir(Props{prop_i, 2});
                            %Sort in natural alphnumeric order
                            [~, SortInd] = sort_nat({Files.name});
                            Files = Files(SortInd);
                            % Make objects
                            Objs = [];
                            for file_i = 1:numel(Files)
                                [path, filename, ex] = fileparts(Files(file_i).name);
                                if strcmpi(ex, '.ncs')
                                    csc = CscObj('file', [Props{prop_i,2} filename ex]);
                                    csc.Path = Props{prop_i,2};
                                    Objs = [Objs csc];
                                end
                            end
                            if isempty(Objs)
                                %No csc files found on path
                                fprintf('No .ncs files found on %s\n', Props{prop_i, 2})
                                return
                            else
                            Obj = Objs;
                            end
                        end
                    case 'objname'
                        for csc_i = 1 : numel(Props{prop_i,2})
                            Obj(csc_i).Name = Props{prop_i,2}{csc_i};
                        end
                        if NlxAreWeConnected
                            for csc_i = 1 :numel(Obj)
                                [~, BitVolts] =  NlxSendCommand(['-GetADBitVolts ' Obj(csc_i).Name]);
                                [~, SamplingFreq] = NlxSendCommand(['-GetSampleFrequency ' Obj(csc_i).Name]);
                                Obj(csc_i).ADBitVolts = str2double(BitVolts);
                                Obj(csc_i).SamplingFrequency = str2double(SamplingFreq);
                                Obj(csc_i).Fs = Obj(csc_i).SamplingFrequency;
                            end
                        else
                            warning('Not connected to Cheetah, BitVolt and Fs not set')
                        end %NlxAreWeConnected
                    case 'objtype'
                        for csc_i = 1 : numel(Props{prop_i,2})
                            Obj(csc_i).cheetahType = Props{prop_i,2}{csc_i};
                        end
                    otherwise
                        error(['arg set ' num2str(prop_i) ' was not recognised'...
                            'Valid args are: ''File '' ''Path'' ''ObjName'' ''ObjType'' ' ]);
                end %switch
            end %for 
        end%Instantiation 
        
        % % The output args are optional, since CSC is a handle the data is passed as a reference. 
            % Use output args for clarity as appropriate.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Stream Data from Cheetah
        function Obj = Stream(Obj)
            % [Success, obj] = Stream(objList)
            % Stream an array CSC objects from Cheetah. 
            % Assumes: 
            % Connection is established
            % obj(i).Name is CSCAcqEnt obj string returned from NlxGetCheetahObjects, verification should be done prior to call
            if ~NlxAreWeConnected
                error('Not connected to Cheetah')
            else
                ValidCSCs = [];
                SuccessData = [];
                for obj_i = 1:numel(Obj)
                    if strcmpi(Obj(obj_i).cheetahType, 'CscAcqEnt') && ischar(Obj(obj_i).Name)
                        [Success(obj_i), Obj(obj_i).DataArray, Obj(obj_i).TimeStampArray, ChNum, ...
                            ~, Obj(obj_i).NumValidSamplesArray, Obj(obj_i).NumRecordsReturned, ...
                            Obj(obj_i).NumRecordsDropped ] = NlxGetNewCSCData(Obj(obj_i).Name);
                        ValidCSCs = [ValidCSCs obj_i];
                        if isempty(Obj(obj_i).Fs)
                            [~, Fs] = NlxSendCommand(['-GetSampleFrequency ' Obj(obj_i).Name]);
                            Obj(obj_i).Fs = str2double(Fs);
                            Obj(obj_i).SamplingFrequency = Obj(obj_i).Fs;
                        end
                        if isempty(Obj(obj_i).ADBitVolts)
                            [~, ADBitVolts] = NlxSendCommand(['-GetADBitVolts ' Obj(obj_i).Name]);
                            Obj(obj_i).ADBitVolts = str2double(ADBitVolts);
                        end
                        Obj(obj_i).Units = 'ADC';
                    else
                        warning('Acquisition type or name not recognised');
                        continue
                    end
                    SuccessData = [SuccessData Obj(obj_i).NumRecordsReturned];
                end
                if ~all(SuccessData)
%                     GUIMessage(handles, 'One or more cscs did not return data', 'blink', 3);
                    warning('One or more cscs did not return data');
                end
                Obj(obj_i).ChannelNumberArray = ChNum;
                if ~any(Success)
                    warning('One or more CSCs were not successfully streamed from Cheetah');
                end
            end
            Obj = Obj(ValidCSCs);
        end %stream     
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function [Obj Rec2Drop] = AlignRec(Obj)
        % [Data Rec2Drop] = AlignRec(Data)
        % Assumes Data is structure of CSC!!!!
        % Align all datarecords with Most lagging stream. Return with Aligned
        % records, Corrected Returned, Dropped, and vector of adjustment, Rec2Drop
        % Rev 4/3/12 Rec2Drop could be non Integer 
        % Rec2Drop = (m-Time-Time)/... -> Rec2Drop = round((m-Time-Time)/...
            Time = zeros(1,numel(Obj)-1);
            try
            for stream_i = 1:numel(Obj)
%                 if isfield(Obj, 'eventIDArray')
                if strcmpi(Obj(stream_i).cheetahType, 'eventIDArray')
                    warning('AlignRec assumes CSC only')
                    continue
                end
                Time(stream_i) = Obj(stream_i).TimeStampArray(1);
            end
            catch
                1;
            end
            [m i]= max(Time);
            RecLen = length(Obj(1).DataArray)/length(Obj(1).TimeStampArray);
            Rec2Drop = round((m-Time)/(RecLen/double(Obj(1).SamplingFrequency)*1e6));
            for stream_i = 1:numel(Obj)
%                 if isfield(Obj, 'eventIDArray')
                if strcmpi(Obj(stream_i).cheetahType, 'eventIDArray')
                    warning('AlignRec assumes CSC only')
                    continue
                end
                RecBegin = min([Rec2Drop(stream_i)+1 length(Obj(stream_i).TimeStampArray)]);
                Obj(stream_i).DataArray = Obj(stream_i).DataArray(Rec2Drop(stream_i)*RecLen+1:end);
                Obj(stream_i).TimeStampArray = Obj(stream_i).TimeStampArray(RecBegin:end);
                Obj(stream_i).ChannelNumberArray = Obj(stream_i).ChannelNumberArray(RecBegin:end);
                Obj(stream_i).NumRecordsReturned = Obj(stream_i).NumRecordsReturned - (Rec2Drop(stream_i));
                Obj(stream_i).NumRecordsDropped = Obj(stream_i).NumRecordsDropped + Rec2Drop(stream_i);
            end
        end %function
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Get from file    
        function varargout = Read(Obj, varargin)
        %[Success, obj] = Read(Obj, prop1, VAL1, prop2, VAL2)
        % Assumes: 
        % obj(i).File is CSC file
        % (prop, VAL) pairs: 
        %   'ExtractMode', Modes listed below,
        %   'ExtractModeVector', numeric vector
        %   'file', FileName.ncs
        %   'path', Path\
        %'ModeArray', Vector or 'ExtractModeVector', integer 1-5,
        % Extract Mode/Array:
        % 1 OR all, default     /   value is ignored
        % 2 OR Index Range      /   two indices range of records 
        % 3 OR Index List       /   indice list of records
        % 4 OR Timestamp Range  /   two indices range of timestamps
        % 5 OR Timestamp List   /   indice list of timestamps
        % see also Nlx2MatCSC documentation
        
        
        %Verify that varargin consists of property, VAL pairsCSC
        Obj = SetParameterPairs(Obj, varargin{:});
        
               
        TimeStamps      =   true;
        ChNum           =   true;
        SampleFreq      =   true;
        ValSamples      =   true;
        Samples         =   1;
        FieldSelection = [TimeStamps ChNum SampleFreq ValSamples Samples];

        ExtractHeader   = 1;
         for obj_i = 1:numel(Obj)
             if isempty(dir([Obj(obj_i).Path Obj(obj_i).File Obj(obj_i).Ext]))
                 error([Obj(obj_i).Path Obj(obj_i).File Obj(obj_i).Ext ' Not found']);
             end
            if strcmpi(Obj(obj_i).cheetahType, 'CSCAcqEnt')
               [Obj(obj_i).TimeStampArray, Obj(obj_i).ChannelNumberArray, Obj(obj_i).SamplingFrequency, Obj(obj_i).NumValidSamplesArray, Obj(obj_i).DataArray, Obj(obj_i).Header] = ...
                    Nlx2MatCSC([Obj(obj_i).Path Obj(obj_i).File Obj(obj_i).Ext], FieldSelection, ExtractHeader, Obj(obj_i).ExtractMode, Obj(obj_i).ModeArray );
%                 Obj(obj_i).DataArray = reshape(Obj(obj_i).DataArray, [ ], 1);
                Obj(obj_i).Fs = Obj(obj_i).SamplingFrequency(1);
                %Default for name is file name. In future there will be a
                %method that allows the user to customize the name.
                if isempty(Obj(obj_i).Name)
                    Obj(obj_i).Name = [Obj(obj_i).File num2str(obj_i)];
                end
                try
                    Hdr = ExtractNlxHeader(Obj(obj_i).Header);
                    Obj(obj_i).ADBitVolts = Hdr.ADBitVolts;
                catch
                    warning('Error reading file')
                end
            end
%             fprintf('.');
        end % for file = obj
%         fprintf('\n');
        switch nargout
            case 0
                varargout{1} = Obj;
            case 1
                varargout{1} = Obj;
            case 2 
                varargout{1} = Obj;
                varargout{2} = Obj;
        end % switch
    end %Read
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Write to file    
    function Write(Obj, varargin)
    %[Success, obj] = Write(CSCFile, prop1, VAL1, prop2, VAL2)
    % Assumes: 
    % prop, VAL is 'ModeArray', Vector or 'ExtractMode', integer 1-5
    warning('Writing header causes Massive ML error. Therefore, it is disabled')
    %[Success, obj] = Read(Obj, prop1, VAL1, prop2, VAL2)
    % Assumes: 
    % obj(i).File is CSC file
    % (prop, VAL) pairs: 
    %   'ExtractMode', Modes listed below,
    %   'ExtractModeVector', numeric vector
    %   'file', FileName.ncs
    %   'path', Path\
    %'ModeArray', Vector or 'ExtractModeVector', integer 1-5,
    % Extract Mode/Array:
    % 1 OR all, default     /   value is ignored
    % 2 OR Index Range      /   two indices range of records 
    % 3 OR Index List       /   indice list of records
    % 4 OR Timestamp Range  /   two indices range of timestamps
    % 5 OR Timestamp List   /   indice list of timestamps
    % see also Nlx2MatCSC documentation


    Obj = SetParameterPairs(Obj, varargin{:});
    
    
    AppendToFileFlag = 0; % false
    ExportMode = 1; % All
    ExportModeVector = 1; % ignore

    TimeStamps  = true;
    ChNum                  = true;
    SampleFreq        = true;
    ValSamples            = 1;
    Samples               = 1;
    WriteHeader = 0;
    FieldSelectionFlags = [TimeStamps ChNum SampleFreq ValSamples Samples];
    Obj = Volts2AD(Obj);

    for obj_i = 1:numel(Obj)
        Obj(obj_i) = VerifyCSC(Obj(obj_i));
        NRec = numel(Obj(obj_i).DataArray)/512;
        if logical(NRec-round(NRec))
            warning('Number of Samples is not a mult of record size, 512')
        end
        if isempty(Obj(obj_i).SamplingFrequency) || length(Obj(obj_i).SamplingFrequency) ~= NRec
            Obj(obj_i).SamplingFrequency = repmat(Obj(obj_i).Fs, 1, NRec);
        end
        if strcmpi(Obj(obj_i).cheetahType, 'CSCAcqEnt')
           Mat2NlxCSC([Obj(obj_i).Path Obj(obj_i).File '.ncs'], ...
               AppendToFileFlag, ExportMode, ExportModeVector, FieldSelectionFlags, ...
               Obj(obj_i).TimeStampArray, Obj(obj_i).ChannelNumberArray, Obj(obj_i).SamplingFrequency, Obj(obj_i).NumValidSamplesArray, ...
               reshape(Obj(obj_i).DataArray, 512, []), ...
               Obj(obj_i).Header);
        end
%             fprintf('.');
    end % for file = obj
%         fprintf('\n');
    end %Write

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Obj = VerifyCSC(Obj)
    % Check all required parameters, create if necessay
        for obj_i = 1: length(Obj(:)')
            % Verify if Obj is writeable
            if isempty(Obj(obj_i).cheetahType)
                error('cheetahType is not defined');
            end
            if isempty(Obj(obj_i).DataArray)
                error('obj has no data')
            end
            if isempty(Obj(obj_i).Fs)
                error('obj must have Fs, sampling frequency, length > 0')
            end
            if strcmpi(Obj(obj_i).Units, 'ADC') && abs(Obj(obj_i).DataArray(1)) < 1-eps
                error(['Obj.' Obj(obj_i).Units ' is inconsistent with Obj.DataArray'])
            elseif strcmpi(Obj(obj_i).Units, 'Volts') && Obj(obj_i).DataArray(1) > 1
                fprintf('\n')
                warning(['obj.' Obj(obj_i).Units ' appears inconsistent with obj.DataArray'])
                fprintf('\n')
            end
            if strcmpi(Obj(obj_i).cheetahType, 'CSCAcqEnt')
                NRec = length(Obj(obj_i).DataArray(:))/512;
                %test for non 512 rec size
                if logical(round(NRec) - NRec)
                    warning('Padding data to multiple of 512');
                    NRec = ceil(NRec);
                    Obj(obj_i).DataArray = [Obj(obj_i).DataArray(:)' zeros(1, NRec*512-length(Obj(obj_i).DataArray))];
                else
                    Obj(obj_i).DataArray = Obj(obj_i).DataArray;
                end
                Obj(obj_i).DataArray  = double(reshape(Obj(obj_i).DataArray, 512, []));
                % Verify TimeStampArray
                if length(Obj(obj_i).TimeStampArray) ~= NRec
                    TSA = CSCTimeStampArray(0, Obj(obj_i).Fs(1), length(Obj(obj_i).DataArray(:)));
                    Obj(obj_i).TimeStampArray = double(TSA(1:512:end));
                else
                    Obj(obj_i).TimeStampArray  = double(Obj(obj_i).TimeStampArray );
                end
                
                % Verify ChannelNumberArray
                if length(Obj(obj_i).ChannelNumberArray) ~= NRec
                    if isempty(Obj(obj_i).ChannelNumberArray)
                        Obj(obj_i).ChannelNumberArray = ones(1, NRec);
                    else
                        Obj(obj_i).ChannelNumberArray = repmat(Obj(obj_i).ChannelNumberArray(1), 1, NRec);
                    end
                end
                Obj(obj_i).ChannelNumberArray = double(Obj(obj_i).ChannelNumberArray);
                % Verify Sampling Frequency Array
                if length(Obj(obj_i).Fs) ~= NRec
                    Obj(obj_i).Fs = Obj(obj_i).Fs(1);
                    Obj(obj_i).SamplingFrequency = repmat(Obj(obj_i).Fs(1), 1, NRec);
                end
                % Verify NumValidSamplesArray
                if numel(Obj(obj_i).NumValidSamplesArray) ~= NRec
                    Obj(obj_i).NumValidSamplesArray = repmat(512, 1, NRec);
                end 
                Obj(obj_i).NumValidSamplesArray = double(Obj(obj_i).NumValidSamplesArray);
                % Make sure Header has minimum entries
                Obj(obj_i) = MakeHeader(Obj(obj_i));
            end %for obj_i = 1:numel(Obj)
    %             fprintf('.');
        end % for file = obj
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %Convert to Volts
        function varargout = AD2Volts(obj)
            for i = 1:numel(obj)
                if strcmpi(obj(i).Units, 'ADC')
                    obj(i).DataArray = double(obj(i).DataArray)*obj(i).ADBitVolts;
                    obj(i).Units = 'Volts';
                end
            end
            if nargout == 1
               varargout{1} = obj;
            end 
        end %AD2Volts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %Convert To  ADInt64
        function varargout = Volts2AD(obj)
            for i = 1:numel(obj)
                if strcmpi(obj(i).Units, 'Volts')
                    obj(i).DataArray = (obj(i).DataArray/obj(i).ADBitVolts);
                    obj(i).Units = 'ADC';
                    if abs(max(obj(i).DataArray(:)) > 32768)
                        warning(['Obj ' num2str(i) ' is out of AD Range'])
                    end
                end
            end
            if nargout == 1
               varargout{1} = obj;
            end  
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %Create Audio object of Cscs
        function SndObj = CSC2Audio(obj, Scale, AudioID)
            obj = AD2Volts(obj);
            %validate Scale arg
            if nargin < 2 || isempty(Scale)
                Scale = 132e-3;
            elseif ~isscalar(Scale) && abs(Scale) > 132e-3 
                    error('Second arg to CSC2Wav must be scalar < 132mV')
            end
            %Create audio object                       
            for csc = 1:numel(obj)
                %Scale for audio and create object
                if nargin < 3
                    %HaveID = false
                    SndObj(csc) = audioplayer(obj(csc).DataArray(:)/Scale, obj(csc).Fs, 16);
                else
                    %HaveID = true
                    SndObj(csc) = audioplayer(obj(csc).DataArray(:)/Scale, obj(csc).Fs, 16, AudioID);
                end
            end
        end %End CSC2Wav

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %threshold CSCs with defined function
        function obj = Threshold(obj, varargin)
            for csc = 1:numel(obj)
                if numel(varargin) == 0
                    obj(csc).DataArray = Threshold(obj(csc).DataArray, @DblDiode);
                else
                    obj(csc).DataArray = Threshold(obj(csc).DataArray, varargin{:});
                end % if numel
            end %end for obj
        end %end Threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %length
        function [L, T]  = DataLength(obj)
            %Returns Length in samples and if Fs is defined, length in Time
            for csc = 1:size(obj(:),1)
                %Since we are redefining length for CscObj, we must use size
                L(csc) = size(obj(csc).DataArray(:),1);
                if ~isempty(obj(csc).Fs)
                    T(csc) = L(csc)/obj(csc).Fs(1);
                end
            end
        end %length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %Get input Range
        function obj = GetInputRange(obj)
            if NlxAreWeConnected
               for csc = 1:length(obj(:))
                   if ~isempty(obj(csc).Name)
                       [~, IR] = NlxSendCommand(sprintf('-GetInputRange %s',obj(csc).Name));
                       obj(csc).InputRange = str2double(IR);
                   else
                       warning(['obj Name is not specified. Skipping InputRange for ' num2str(csc) ' object.']);
                   end
               end
            end
        end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Create Minimal Header
        function Obj = MakeHeader(Obj)
            SP = ' ';
            for obj_i = 1:length(Obj(:))
                 HdrStr = reshape(lower(char(Obj(obj_i).Header))', 1,[ ]);
                if isempty(HdrStr)
                    Obj(obj_i).Header = {'######## Neuralynx Data File Header'; ['## Minimum header created on ' date]};
                    HdrStr = reshape(lower(char(Obj(obj_i).Header))', 1,[ ]);
                end
                if isempty(strfind(HdrStr, '-filetype'))
                    Obj(obj_i).Header(end+1, 1) = {'	-FileType CSC'};
                end
                if isempty(strfind(HdrStr, '-samplingfrequency'))
                    Obj(obj_i).Header(end+1, 1) = {['	-SamplingFrequency' SP sprintf('%6.0f',Obj(obj_i).Fs(1)) ]};
                end
                if isempty(strfind(HdrStr, '-adbitvolts'))
                    Obj(obj_i).Header(end+1, 1) = {['	-ADBitVolts' SP sprintf('%2.10f',(Obj(obj_i).ADBitVolts))]};
                end
                if isempty(strfind(HdrStr, '-admaxvalue'))
                    Obj(obj_i) = Volts2AD(Obj(obj_i));
                    MaxAD = 32767;%max(o.DataArray(:));
                    Obj(obj_i).Header(end+1, 1) = {['	-ADMaxValue' SP sprintf('%6.0f', MaxAD)]};
                end
                if isempty(strfind(HdrStr, '-acqentname'))
                    Obj(obj_i).Header(end+1, 1) = {['	-AcqEntName' SP	Obj(obj_i).Name]};
                end
                if isempty(strfind(HdrStr, '-numadchannels'))
                    Obj(obj_i).Header(end+1, 1) = {['	-NumADChannels' SP sprintf('%3.0f', numel(Obj))]};
                end
                if isempty(strfind(HdrStr, '-inputrange'))
                    Obj(obj_i).Header(end+1, 1) = {['	-InputRange 1000']};% SP sprintf('%3.0f', numel(Obj)])};
                end
                if isempty(strfind(HdrStr, '-recordsize'))
                    Obj(obj_i).Header(end+1, 1) = {'	-RecordSize 1044'};
                end
            end
        end% MakeHeader
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    function Obj = SetParameterPairs(Obj, varargin)
            %Verify that varargin consists of property, VAL pairsCSC
    if mod(numel(varargin), 2)
        error('input arguments must be (''property'', value) pairs');
    else
        %Rearrange so that there are n rows of [property, VAL]
        Props = reshape(varargin, 2, [])';
    end 
    for csc_i = 1:length(Obj)
    % Update all properties, do extract mode first
        for prop_i = 1:size(Props, 1)
            % execute code based on property, leave in switch syntax to
            % allow for easy addition of properties,Val sets for Read
            switch lower(Props{prop_i,1})
                case 'extractmode'
                   % Val = 1, 2, 3, 4, or 5 see Nlx2MatCSC documentation
                   if ischar(Props{prop_i, 2})
                        [ValidMode, ModeN] = ismember(lower(Props{prop_i, 2}), ...
                                {'all', 'index range', 'index list', ...
                                'timestamp range', 'timestamp list'});
                    else
                        [ValidMode, ModeN] = ismember(Props{prop_i, 2}, 1:5);
                    end
                    if  ValidMode
                        for csc_i = 1 : numel(Obj)
                            Obj(csc_i).ExtractMode = ModeN;
                        end
                    else
                        error('Extract Mode = scalar 1:5 or ''all'',''index range'', ''index list'', ''timestamp range'', ''timestamp list'' ')
                    end
                case 'extractmodevector'
                    if ~isnumeric(Props{prop_i,2})
                        error('Mode Vector must be numeric array')
                    end
                    for csc_i = 1 : numel(Obj)
                        %Store as 2 number vector or list
                        if Obj(csc_i).ExtractMode == 2 ||  Obj(csc_i).ExtractMode == 4
                            if length(Props{prop_i, 2}) < 2
                                error('Record Range mode (mode 2) requires 2 parameters.  1 were specified.')
                            end
                            Obj(csc_i).ModeArray = [Props{prop_i, 2}(1) Props{prop_i, 2}(end)];
                        else
                            Obj(csc_i).ModeArray = Props{prop_i, 2};
                        end
                    end
                case 'file'
                    Obj(csc_i).File = [Props{prop_i, 2} num2str(csc_i)];
                case 'path'
                    Obj(csc_i).Path = Props{prop_i, 2};
                otherwise
                    error(['CscObj Read(CSC, ' num2str(Props{prop_i,1}) ' val not recognised']);
            end %switch
        end %for Prop
    end %csc
    end %SetParameterPairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function Hdr = ExtractNlxHeader(NlxHdr)
    % Extracts parameters from csc and nrd Headers and creates structure for access 
    %% Intitalize structure
    Hdr.AD2Volts = [ ];     %double
    Hdr.ADChan = [ ];       %int singlar for CSC, vector of int for nrd
    Hdr.ADGain   = [ ];     %double
    Hdr.Gain = [ ];         %double
    Hdr.Date = [ ];         %Date/Time Vector [yyyy mm dd hh mm ss.ddd]
    Hdr.DateOpen = [ ];     %Date/Time Vector [yyyy mm dd hh mm ss.ddd]
    Hdr.HPF = [ ];          %double
    Hdr.Threshold = [ ];    %double
    Hdr.LPF  = [ ];         %double
    Hdr.MaxAD = [];         %int
    Hdr.NChan = [];         %int singular for CSC, length of Hdr.ADChan for nrd
    Hdr.SamFreq = [];       %double
    Hdr.Ver = [];           %String
    Hdr.ADBitVolts = [];    %double

    % parse header and put in structure
    %%
    for row = 1:length(NlxHdr)
    % AD Bit Volts, ADSample to volt conversion
        if ~isempty(strfind(lower(NlxHdr{row}), 'bitvolts'));    
            Hdr.AD2Volts = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
        end

    % CSC Channel
        if ~isempty(strfind(lower(NlxHdr{row}), 'adchan'));
            Hdr.ADChan = str2double(regexprep(NlxHdr{row},'0-9',''));
            Hdr.ADChan = cast(Hdr.ADChan, 'uint16');
        end    
    % AD Gain
        if ~isempty(strfind(lower(NlxHdr{row}), 'adgain'));    
            Hdr.ADGain = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
        end

    % Amp Gain    
        if ~isempty(strfind(lower(NlxHdr{row}), 'ampgain'));
            Hdr.Gain = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
        end

    % Date Time Closed
        if ~isempty(strfind(lower(NlxHdr{row}), 'time closed'))
            Date = regexp(NlxHdr{row},'(\d|\d\d)/(\d|\d\d)/\d\d\d\d','match');
            Time = regexp(NlxHdr{row},'(\d|\d\d):(\d|\d\d):(\d|\d\d).\d\d\d','match');
            Hdr.Date = datevec([Date{:} ' ' Time{:}]);
        end

    % Date Time Open
        if ~isempty(strfind(lower(NlxHdr{row}), 'time open'))
            Date = regexp(NlxHdr{row},'(\d|\d\d)/(\d|\d\d)/\d\d\d\d','match');
            Time = regexp(NlxHdr{row},'(\d|\d\d):(\d|\d\d):(\d|\d\d).\d\d\d','match');
            Hdr.DateOpen = datevec([Date{:} ' ' Time{:}]);
        end

    % File Name   
        if ~isempty(strfind(lower(NlxHdr{row}), 'file name'));
            %to be done later. cleary the user knows the file name
        end

    % High Pass Filter
        if ~isempty(strfind(lower(NlxHdr{row}), 'hpf'));
            Hdr.HPF = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
        end

    % Threshold Value
        if ~isempty(strfind(lower(NlxHdr{row}), 'threshold'));  
            Hdr.Threshold = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
        end

    % Low Pass Filter
        if ~isempty(strfind(lower(NlxHdr{row}), 'lpf'));
            Hdr.LPF =  str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
        end

    % Maximum AD Value
        if ~isempty(strfind(lower(NlxHdr{row}), 'maxadval'));
            Hdr.MaxAD = str2double(regexprep(NlxHdr{row},'\D',''));
            Hdr.MaxAD = cast(Hdr.MaxAD, 'uint32');
        end

    % Number of Channels
        if ~isempty(strfind(lower(NlxHdr{row}), 'channels'));
            Hdr.NChan = str2double(regexprep(NlxHdr{row},'\D',''));
            Hdr.NChan = cast(Hdr.NChan, 'uint16');
        end

    % Sampling Frequency, Fs
        if ~isempty(strfind(lower(NlxHdr{row}), 'frequency'));
            Hdr.SamFreq = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
        end

    % Cheetah Version
        if ~isempty(strfind(lower(NlxHdr{row}), 'cheetah'));
            Hdr.Ver = regexprep(NlxHdr{row},'[^0-9^\.]','');
        end

    % BitVolts
        if ~isempty(strfind(lower(NlxHdr{row}), 'adbitvolts'))
            Hdr.ADBitVolts = regexprep(NlxHdr{row},'[^0-9^\.]','');
        end
    end
    end %Extract Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end %methods
    
end %classdef