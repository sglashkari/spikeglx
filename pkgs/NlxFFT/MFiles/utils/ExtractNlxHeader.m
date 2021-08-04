function Hdr = ExtractNlxHeader(NlxHdr)
% Extracts parameters from csc and nrd Headers and creates structure for access 
% Updated 2013 08 07 to make backwards compatable with ReadNlxHeader
% Needs updating for spike files (.nse, .nst, .ntt) and video files, (.nvt)
% Need to double check CheetahCmds with Reference Guide
%% Intitalize structure

Hdr.FileName                = []; %Legacy ReadNlxHeader
Hdr.FileType                = [];       
Hdr.FileVersion             = [];      
Hdr.CheetahRev              = [];
Hdr.TimeOpen                = [];
Hdr.TimeClosed              = [];
Hdr.HardwareSubSystemType   = [];
Hdr.HardwareSubSystemName   = [];
Hdr.NumADChannels           = [];
Hdr.ADBitVolts              = [];
Hdr.AD2Volts                = []; %Legacy ReadNlxHeader
Hdr.MaxAD                   = []; %Legacy ReadNlxHeader
Hdr.ADMaxValue              = [];
Hdr.SamplingFrequency       = [];
Hdr.RecordSize              = [];
Hdr.AcqEntName              = [];
Hdr.ADChan                  = [];
Hdr.ADGain                  = [];
Hdr.Gain                    = [];
Hdr.InputRange              = [];
Hdr.Threshold               = [];
Hdr.InputInverted           = [];
Hdr.DspDelayCompensation    = [];
Hdr.DspFilterDelay_us       = [];
Hdr.DspLowCutFilterType     = [];
Hdr.DspLowCutFilterEnabled  = [];
Hdr.DspLowCutNumTaps        = [];
Hdr.DSPHighCutFilterType    = [];
Hdr.DSPHighCutFilterEnabled = [];
Hdr.DspHighCutNumTaps       = [];
Comments                    = []; %Cell of String
MatlabCmds                  = []; %Cell of String
CheetahCmds                 = []; %Cell of String
Unknown                     = []; %Unknown commands

%% parse header and put in structure

for row = 1:length(NlxHdr)
    CurrentLine = strtrim(NlxHdr{row});
    if ~isempty(CurrentLine)
    % Comments
        if strcmp(CurrentLine(1), '#')
            Comments = [Comments; {CurrentLine}]; %#ok<AGROW> %
        end
    % Matlab
        if strcmp(CurrentLine(1), '~')
            MatlabCmds = [MatlabCmds; {CurrentLine}]; %#ok<AGROW> %#
            continue;
        end
    % CheetahCmds List
        if strcmp(CurrentLine(1), '-')
            CheetahCmds = [CheetahCmds; {CurrentLine}];
        end
    %% Cheetah Cmds and Queries  
        Strs = strsplit(CurrentLine);
        Cmd = Strs{1}(2:end);
        CmdArg = Strs{2};
        %% General information
        %Date
        % Date Time Open
        if ~isempty(strfind(lower(CurrentLine), 'time open'))
            Date = regexp(NlxHdr{row},'(\d|\d\d)/(\d|\d\d)/\d\d\d\d','match');
            Time = regexp(NlxHdr{row},'(\d|\d\d):(\d|\d\d):(\d|\d\d).\d\d\d','match');
            Hdr.TimeOpen = datevec([Date{:} ' ' Time{:}]);
            continue;
        end    
        % Date Time Closed
        if ~isempty(strfind(lower(CurrentLine), 'time closed'))
            Date = regexp(NlxHdr{row},'(\d|\d\d)/(\d|\d\d)/\d\d\d\d','match');
            Time = regexp(NlxHdr{row},'(\d|\d\d):(\d|\d\d):(\d|\d\d).\d\d\d','match');
            Hdr.TimeClosed = datevec([Date{:} ' ' Time{:}]);
            continue;
        end
        if strcmpi(Strs{2}, 'file')
            Hdr.FileName = CurrentLine(strfind(CurrentLine, 'File Name')+10:end);
            continue;
        end
        switch lower(Cmd)
            case 'file name'
            %to be done later. cleary the user knows the file name
            case 'filetype'
    %             Strs = strsplit(NlxHdr{row});
                Hdr.FileType = Strs{2};
            case 'fileversion'
                Hdr.FileVersion = CmdArg;%regexprep(NlxHdr{row},'[^0-9^\.]','');
            case 'cheetahrev'
                Hdr.CheetahRev = regexprep(NlxHdr{row},'[^0-9^\.]','');
            case 'time open'
                Date = regexp(NlxHdr{row},'(\d|\d\d)/(\d|\d\d)/\d\d\d\d','match');
                Time = regexp(NlxHdr{row},'(\d|\d\d):(\d|\d\d):(\d|\d\d).\d\d\d','match');
                Hdr.DateOpen = datevec([Date{:} ' ' Time{:}]);
            case'time closed'
                Date = regexp(NlxHdr{row},'(\d|\d\d)/(\d|\d\d)/\d\d\d\d','match');
                Time = regexp(NlxHdr{row},'(\d|\d\d):(\d|\d\d):(\d|\d\d).\d\d\d','match');
                Hdr.DateClosed = datevec([Date{:} ' ' Time{:}]);
            case 'hardwaresubsystemtype'
    %             Strs = strsplit(NlxHdr{row});
                Hdr.HardwareSubSystemType = Strs{2};
            case 'hardwaresubsystemname'
    %             Strs = strsplit(NlxHdr{row});
                Hdr.HardwareSubSystemName = Strs{2};
            case 'numadchannels'
                Hdr.NumADChannels = str2double(CmdArg);%regexprep(NlxHdr{row},'\D',''));
                Hdr.NumADChannels = cast(Hdr.NumADChannels, 'uint16');
            case 'adbitvolts'
                Hdr.ADBitVolts = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'bitvolts' %legacy
                Hdr.AD2Volts = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'maxadval' %legacy
                Hdr.MaxAD = str2double(CmdArg);%regexprep(NlxHdr{row},'\D',''));
                Hdr.MaxAD = cast(Hdr.MaxAD, 'uint32');
            case 'admaxvalue'
                Hdr.ADMaxValue = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'samplingfrequency'
                Hdr.SamplingFrequency = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'recordsize'
                Hdr.RecordSize = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]','')));

            %% Channel Specific
            case 'acqentname'
    %             Strs = strsplit(NlxHdr{row});
                Hdr.AcqEntName = Strs{2};
            case 'adchannel'
                Hdr.ADChannel = str2double(CmdArg);%regexprep(NlxHdr{row},'0-9',''));
                Hdr.ADChannel = cast(Hdr.ADChan, 'uint16');
            case 'adgain'
                Hdr.ADGain = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]','')));
            case 'ampgain'
                Hdr.Gain = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]','')));
            case 'inputrange'
                Hdr.InputRange = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'threshold'
                Hdr.Threshold = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]','')));
            case 'inputinverted'
                if ~isempty(CmdArg)
                    Hdr.InputInverted = logical(str2num(lower(CmdArg)));%regexprep(NlxHdr{row},'truefalse','');
                else 
                    Hdr.InputInverted = [];
                end
            %% Filter Settings

            %General DSP
            case 'dspdelaycompensation'
                Hdr.DspDelayCompensation = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'dspfilterdelay_µs'
                Hdr.DspFilterDelay_us = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]',''));
            %Low
            case 'dsplowcutfrequency'
                Hdr.DspLowCutFrequency = str2double(CmdArg);%regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'dsplowcutfiltertype'
                Hdr.DspLowCutFilterType = CmdArg;%regexprep(NlxHdr{row},'[^0-9^\.]','');
            case 'dsplowcutfilterenabled'
                if ~isempty(CmdArg)
                    Hdr.DSPLowCutFilterEnabled = logical(str2num(lower(CmdArg)));%#ok<ST2NM> %regexprep(NlxHdr{row},'[^0-9^\.]','');
                else
                    Hdr.DSPLowCutFilterEnabled =  [];
                end
            case 'dsplowcutnumtaps'
                Hdr.DspLowCutNumTaps = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
            %High
            case 'dsphighcutfrequency'
                Hdr.DspHighCutFrequency = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'dsphighcutfiltertype'
                Hdr.DspHighCutFilterType = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
            case 'dsphighcutfilterenabled'
                if ~isempty(CmdArg)
                    Hdr.DSPHighCutFilterEnabled = logical(str2num(lower(CmdArg)));%#ok<ST2NM> %regexprep(NlxHdr{row},'[^0-9^\.]','');
                else
                    Hdr.DSPHighCutFilterEnabled = [];
                end
            case 'dsphighcutnumtaps'
                Hdr.DspHighCutNumTaps = str2double(regexprep(NlxHdr{row},'[^0-9^\.]',''));
            otherwise
                if ~strcmp(NlxHdr{row}(1), '#')
                    Unknown = [Unknown NlxHdr(row)]; %#ok<AGROW>
                end

        end%switch
    end %if isempty
end % for
%Legacy See ReadNlxHeader which was derived from Phil's code
Hdr.MatlabCmds = MatlabCmds;
Hdr.CheetahCmds = CheetahCmds;
Hdr.Comments = Comments;
Hdr.Unknown = Unknown;
if ~isempty(Unknown)
    warning('Unknown entries found in header')
end

 











