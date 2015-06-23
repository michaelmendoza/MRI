function [ fileHeaders ] = readProtHeadFromFiles( structFields )
%READPROTHEADFROMFILES Reads all fields of the structFields into a struct
%from individual files.
%   Detailed explanation goes here
%   
%   AUTHOR: Danny Park
%
%   See also readDatFile, readDatHeaderProtocol, readMdh,
%   readProtHeadFromFiles, writeDatFile, writeDatHeaderProtocol, writeMdh,
%   writeProtHead2Files

for n=1:length(structFields),
    fidin = fopen(structFields{n}, 'rt');
    fileHeaders.(structFields{n}) = fread(fidin,inf,'*char');
    fclose(fidin);
end

end

