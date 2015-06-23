function [ image3d ] = readIceIma( folder, num2read, dim )
%function [ image3d ] = readIceIma( folder, num2read, dim )
%READICEIMA Reads a series of .ima files output by an ice simulation and
%puts them into a 3D matrix for viewing.
%   folder - the folder to find the files in with the name of
%       WriteToFile_####.ima where # represents a single digit (0-9)
%   num2read - the number of .ima to read (starting with 0001)
%   dim - the x and y dimension in pixels of the picture (only square supported)
%   
%   AUTHOR: Danny Park

fileprefix = 'WriteToFile_';
filesuffix = '.ima';

if length(dim) == 1,
    dim(2) = dim(1);
end

image3d = zeros(dim(1),dim(2), num2read, 'int16');

for n=1:num2read,
    fid = fopen([folder ...
                 fileprefix  ...
                 num2str(mod(floor(n/1000),10)) ...
                 num2str(mod(floor(n/100),10)) ...
                 num2str(mod(floor(n/10),10)) ...
                 num2str(mod(n,10)) ...
                 filesuffix], 'rb', 'ieee-le');
    image3d(:,:,n) = reshape(fread(fid, inf, 'int16'), dim);
    fclose(fid);
end

end

