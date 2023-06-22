function msgboxw(message,fontSize)
try
	CreateStruct.Interpreter = 'tex';
	CreateStruct.WindowStyle = 'modal';
% 	CreateStruct.Title = 'MATLAB Message';
	
% 	fontSize = 14;
	% Embed the required tex code in before the string.
	latexMessage = sprintf('\\fontsize{%d} %s', fontSize, message);
	uiwait(msgbox(latexMessage, 'MATLAB message', CreateStruct));
catch ME
	errorMessage = sprintf('Error in msgboxw():\n%s', ME.message);
	fprintf('%s\n', errorMessage);
	uiwait(warndlg(errorMessage));
end
return;