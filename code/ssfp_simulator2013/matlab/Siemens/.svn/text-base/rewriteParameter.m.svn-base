function output = rewriteParameter(protocolString, parameterName, newval)
%   
%   AUTHOR: Danny Park

indx = strfind(protocolString,parameterName);
indx = indx+length(parameterName);

output = protocolString;
for n=1:length(indx),
    opnBrkt = strfind(protocolString(indx(n)+1:end),'{');
    opnBrkt = indx(n)+opnBrkt(1);
    clsBrkt = strfind(protocolString(opnBrkt+1:end),'}');
    clsBrkt = opnBrkt + clsBrkt(1);

%     if (indx(n) < 1) || (opnBrkt < 1) || (clsBrkt < 1),
%         error('not found');
%     end

    output = [output(1:opnBrkt) ' ' newval ' ' output(clsBrkt:end)];
end

end

function val = readParameter(protocol, parameter)
%READPARAMETER Reads a given parameter from the text based protocol.

indx = strfind(protocol,parameter) + length(parameter);
valstr = sscanf(protocol(indx:indx+35)','%s');
val = str2num(valstr(isstrprop(valstr,'digit')));

end