function ErrorNlxFFTIntervalT(handles)
UserData = get(handles.output, 'UserData');
UserData.error = lasterror;
GUIMessage(handles, 'Error in timer fcn. See UserData.error for details', 'fontcolor', 'red', 'blink', 3);