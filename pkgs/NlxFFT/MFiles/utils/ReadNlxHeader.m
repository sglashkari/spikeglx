function Header = ReadNlxHeader(HeaderIn)

Header.Comments = [];
Header.CheetahCmds = [];
Header.Unknown = [];
Header.MatlabCmds = [];
for CurrLine = 1:length(HeaderIn)
    CurrentLine = strtrim(HeaderIn{CurrLine});
    %Determine the first character in the command line.
    if isempty(CurrentLine)
        continue
    else
        CommandTypeCharacter = CurrentLine(1);
    end
    %Determine the command type and make the corresponding action
    if strcmp(CommandTypeCharacter,'#') == 1
        %This line has been commented out to provide information to the
        %user. Skip this line.
        Header.Comments = [Header.Comments HeaderIn(CurrLine)];
    elseif strcmp(CommandTypeCharacter,'-') == 1
        %This corresponds to a Cheetah command. 
        Cmd = HeaderIn{CurrLine};
        spaces = strfind(Cmd, ' ');
        Field = Cmd(2:spaces(1)-1);
        CmdVal = Cmd(spaces(1):end);
        if ~isempty(str2num(CmdVal))
            CmdVal = str2num(CmdVal);
        end
        Header.CheetahCmds.(genvarname(Field)) = CmdVal;
    elseif strcmp(CommandTypeCharacter,'~') == 1
        %This corresponds to a Neuralynx Matlab Command. 
        Header.MatlabCmds = [Header.MatlabCmds HeaderIn(CurrLine)];
    else 
        Header.Unknown = [Header.Unknown HeaderIn(CurrLine)];
    end
end
