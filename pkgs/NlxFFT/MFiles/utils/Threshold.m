function x = Threshold(x, varargin)
% x = Threshold(x, Fcn)
% 2013 03 21
switch numel(varargin)
    case 0
        %Use default function
        Fcn = @DblDiode;
        k = 1;
        x = Fcn(x, k);
    case 1
        Fcn = varargin{1};
        x = Fcn(x);
    otherwise
        Fcn = varargin{1};
        if nargin(varargin{1}) ~= numel([varargin{2:end}])+1
            warning(['It appears that the number of args supplied do not '...
                'match what is required for ' char(varargin{1})])
        end
        x = Fcn(x, varargin{2:end});
end

end%end function

