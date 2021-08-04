function hButton = SetRadioGUI(h, M, N, Names)
% SetRadioGUI(h, M, N, {Names})
% Setup matrix of upto MxN Named radio buttons on figure with handle h 
% if number of Names < M*N then the matrix will not be full
% if number of Names > M*N then the first M*N names are used
% Reccommend N = [ ] so default is used.
if ~ishandle(h)
    %Requires a defined button group
    return
else
    %Clear group in preparation for new set
    delete(get(h, 'children'));
end
if M > 20
    warning('Rows > 20 may not scale well')
end
if nargin< 4
    nButtons = M*N;
    NoNames = true;
else
    nButtons = min(numel(Names), M*N);
    NoNames = false;
end


ButtonXSize = 0.2; 
ButtonYSize = 0.05;
[Row, Col] = meshgrid(linspace(0+ButtonYSize, 1-ButtonYSize, M), linspace(0, 1-ButtonXSize, N));
Row = Row'; Col = Col'; 
Row = Row(end:-1:1,:);
button_i = 1;
while button_i <= nButtons
    if NoNames
        name = ['CSC' num2str(button_i-1)];
    else
        name = Names{button_i};
    end
    hButton(button_i) = uicontrol('Style','checkbox','String', name, 'parent',h, ...
         'HandleVisibility','on', 'units', 'normal', ... 
         'pos', [Col(button_i), Row(button_i), ButtonXSize, ButtonYSize], ...
         'fontsize', 6, 'visible', 'on'); 
     button_i = button_i+1;
end


