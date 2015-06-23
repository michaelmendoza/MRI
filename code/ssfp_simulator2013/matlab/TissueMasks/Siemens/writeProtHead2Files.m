function [  ] = writeProtHead2Files( fileHeaders, prefix )
%WRITEPROTHEAD2FILES Writes all fields of the fileHeaders structs to
%individual files.
%   Detailed explanation goes here
%   
%   AUTHOR: Danny Park
%
%   See also readDatFile, readDatHeaderProtocol, readMdh,
%   readProtHeadFromFiles, writeDatFile, writeDatHeaderProtocol, writeMdh,
%   writeProtHead2Files

if nargin < 2,
    prefix = '';
end

headerFields = fieldnames(fileHeaders);

for n=1:length(headerFields),
    fidout = fopen([prefix headerFields{n}], 'wt');
    fprintf(fidout,'%c',fileHeaders.(headerFields{n}));
    fclose(fidout);
end

end

