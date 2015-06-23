function [ outstr ] = freadstr( file_fid )
%FREADSTR Reads a null terminated string from the given file, consuming the
%null character
%   Detailed explanation goes here
%   
%   AUTHOR: Danny Park


%% Read string
%tempstr = char(max_len);
nc = 0;
while ~feof(file_fid),
    nc = nc + 1;
    tempstr(nc) = fread(file_fid,1,'*char');
    if tempstr(nc) == char(0),
        break;
    end
end

outstr = tempstr(1:nc-1);

end

