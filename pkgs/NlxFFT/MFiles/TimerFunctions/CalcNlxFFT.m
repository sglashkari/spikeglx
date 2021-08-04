function CalcNlxFFT(handles)
% Timer function to calculate and plot ffts for either a streaming object
% or a CSC read from a file.
UserData = get(handles.output, 'UserData');
%Verify that CSCs objects are available to process
if isfieldRecursive(UserData, 'CSCs')
    ChosenCSCs = UserData.ConnectData.ChosenCSCs;
    %% Get Data from source
    if strcmpi(UserData.ConnectData.Input, 'file') 
        % Input is file
        %put time extraction in CscObj
        try
            UserData.ConnectData.CSCs(ChosenCSCs) = Read(UserData.ConnectData.CSCs(ChosenCSCs));
        catch
            GUIMessage(handles, 'Problem reading CSCs from file');
            stop(UserData.TimerData.NlxFFTIntervalT);
            return
        end
    else %Input must be server
        try
            UserData.ConnectData.CSCs(ChosenCSCs) = Stream(UserData.ConnectData.CSCs(ChosenCSCs));
        catch
            GUIMessage(handles, 'Problem streaming CSCs from connection');
            stop(UserData.TimerData.NlxFFTIntervalT);
            return
        end
    end
    %% Convert data to uvolts
    UserData.ConnectData.CSCs(ChosenCSCs) = AD2Volts(UserData.ConnectData.CSCs(ChosenCSCs));
    
    %% One time Setup
    %   Align CSCs, set up figures
    if UserData.FirstTime %Reset when Calc Button is pushed 
        % Align all CSCs with most lagging CSC, (only when online)
        if strcmpi(UserData.ConnectData, 'server')
             UserData.ConnectData.CSCs(ChosenCSCs) = AlignRec(UserData.ConnectData.CSCs(ChosenCSCs));
        end
        %Set up Figures
        if get(handles.CommonZero, 'Value')
            %Single figure with multiple traces
            UserData.hFigs = figure;
        else 
            %Mult plots, 2 options: Multiple figures, one figure multiple
            %subplots
%             for csc_i = 1:numel(ChosenCSCs(:))
%                 % Multiple figures
%                 UserData.hFigs(csc_i) = figure;
%%%%% Alternate code for subplots %%%%%%%%
            % create subplot grid
            NSubPltRow = ceil(sqrt(numel(ChosenCSCs)));
            NSubPltCol = NSubPltRow;
            if numel(ChosenCSCs) <= NSubPltCol * (NSubPltRow-1)
                NSubPltCol = NSubPltCol-1;
            end
            % Generate subplot handles
            UserData.hFigs = figure;
            for csc_i = 1:numel(ChosenCSCs(:))
                UserData.hSubPlts(csc_i) = subplot(NSubPltRow, NSubPltCol, csc_i);
            end
            UserData.COLORORDER = get(UserData.hSubPlts(csc_i), 'ColorOrder');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    %% Get Fig current and find fft
    for csc_i = 1:numel(ChosenCSCs(:))
        if get(handles.CommonZero, 'Value')
            %Allow first plot to set figure properties, (hold off), then
            %hold for subquent plots, 
            h = figure(UserData.hFigs);
            if csc_i == 1
                hold off
            else
                hold all
            end
        else
%             %Mult Figures
%             figure(UserData.hFigs(csc_i))
%%%%% Alternate code for subplots %%%%%%%%
            h = UserData.hSubPlts(csc_i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
        end
        
        %This allows user to change scale
        if ~UserData.FirstTime
            XLIM = get(h, 'xlim');
%             YLIM = get(h, 'ylim');
        end
        dB = get(handles.dB, 'value');
        if dB
            yUnits = 'dBVolt';
        else
            yUnits = 'Volt';
        end
        % Calculate and plot fft on current axis, (gca)
        Myfft(UserData.ConnectData.CSCs(ChosenCSCs(csc_i)).DataArray(:), ...
            UserData.ConnectData.CSCs(ChosenCSCs(csc_i)).Fs, ...
            h, true, dB, 2^15);
        if get(handles.CommonZero, 'Value')
            figure(h)
            ylabel(yUnits);
            grid('on')
        else
            ColorIndex = mod(csc_i-1, size(UserData.COLORORDER, 1))+1;
            set(get(h, 'children'), 'color', UserData.COLORORDER(ColorIndex, :))
            ylabel(h, yUnits);
            grid(h, 'on')
        end
        

        
        %This allows user to change scale
        if ~UserData.FirstTime
            set(h, 'xlim', XLIM)
%             set(h, 'ylim', YLIM)
        end
        
        if ~get(handles.CommonZero, 'Value')
            % Title subplots or multiple figures
            title(h, UserData.ConnectData.CSCs(ChosenCSCs(csc_i)).Name)      
        else
            title('ffts')
        end
    end
    
    % Legend for common Zero plot
    if get(handles.CommonZero, 'Value')
        legend({UserData.ConnectData.CSCs(ChosenCSCs).Name})
    end
    
    % Reset Calc toggle button
    if strcmpi(UserData.TimerData.NlxFFTIntervalT.ExecutionMode, 'singleShot')
        set(handles.CalcFFT, 'value', 0);
    end
    
    % Reset First time flag
    UserData.FirstTime = false;
    
    set(handles.output, 'UserData', UserData);
    
else
    GUIMessage(handles, 'Objects are not available', 'fontcolor', 'red', 'blink', 3)
    stop(UserData.TimerData.NlxFFTIntervalT);
end
pause(0.01)% Allow callbacks to be seen
