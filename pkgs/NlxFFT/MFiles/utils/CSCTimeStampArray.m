function TS_Array = CSCTimeStampArray(arg1, arg2, nSam) %#ok<*FNDEF>
% TS_Array = CSCTimeStampArray(arg1 = Obj or Ti, arg2 = Fs or Ti, arg = num elements)
% Generates a timeIndex, (int64 data type), for a CSC Obj
% Ti - intitial time (us)
% nSam - num of samples
if ~isa(arg1, 'CscObj') 
    % Legacy mode
    % TS_Array = CSCTimeStampArray(Ti, Fs, lenDataArray)
    % Generates a timeIndex, (int64 data type), for a CSC entity
    % Ti - intitial time (us)
    % Fs - sampling Rate
    % nSam - num of samples
    Ti = arg1;
    Fs = arg2;
    Tf = double(arg1) + 1e6*(nSam-1)/double(Fs);
    TS_Array = int64(double(Ti) : 1e6/double(Fs): Tf);
else %Used as method
    % TS_Array = CSCTimeStampArray(CscObj, Ti, nSam)
    for obj_i = 1:length(arg1)
        if isempty(arg1(obj_i).TimeStampArray)
            warning(['CscObj ' num2str(obj_i) 'has an empty TimeStampArray. Result = null.']);
            TS_Array{obj_i} = [];
            nSam = 0;
        end
        if isempty(arg1(obj_i).Fs)
            warning(['CscObj ' num2str(obj_i) 'has an empty samplinging freq, Fs, property. Result = null.']);
            TS_Array{obj_i} = [];
            nSam = 0;
        end
        if arg2 > arg1(obj_i).TimeStampArray(end)
            warning(['Initial Timestamp for CscObj ' num2str(obj_i) 'is greater than last timestamp. Result = null.']);
            TS_Array{obj_i} = [];
            nSam = 0;
        end        
        if nSam>length(arg1(obj_i).DataArray(:))
            warning(['CscObj ' num2str(obj_i) 'has fewer than nSam samples, truncating']);
            nSam = length(arg1(obj_i).DataArray(:));
        end

        if nSam>0
            Fs = arg1(obj_i).Fs;
            Ti = arg2;
            Tf = double(Ti) + 1e6*(nSam-1)/double(Fs);
            TS_Array{obj_i} = int64(double(Ti) : 1e6/double(Fs): Tf);
        end
    end
end