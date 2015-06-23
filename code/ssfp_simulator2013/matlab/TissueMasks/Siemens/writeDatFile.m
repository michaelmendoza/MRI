function [  ] = writeDatFile( rawData, fileHeaders, mdhs, filename )
%function [  ] = writeDatFile( rawData, fileHeaders, mdhs, filename )
%WRITEDATFILE Writes a Siemens raw data file (.dat) like it came from an MRI scanner
%   Detailed explanation goes here
%   
%   AUTHOR: Danny Park
%
%   See also readDatFile, readDatHeaderProtocol, readMdh,
%   readProtHeadFromFiles, writeDatFile, writeDatHeaderProtocol, writeMdh,
%   writeProtHead2Files

%% Constants
MDH_ACQEND = 0 + 1;%0 indexed in C, 1 indexed in MATLAB

%% argument checking
if nargin < 4
  [temp path] = uiputfile('*.dat','Select File to Read');
  filename = [path temp];
end

if length(mdhs) ~= length(rawData)
    error('Number of lines of raw data must be equal to the number of Mdhs');
end

%% using data structure defined on page 150 of the ice manual

dat_fid = fopen(filename,'w','ieee-le');
% dataFieldLoc = fread(dat_fid,1,'int32');
fwrite(dat_fid, 0, 'int32'); % Placeholder for the position of the data

headerFields = fieldnames(fileHeaders);
numHeaderStructs = length(headerFields);
fwrite(dat_fid,numHeaderStructs,'int32');

for n = 1:numHeaderStructs,
    writeDatHeaderProtocol(dat_fid, headerFields{n}, fileHeaders.(headerFields{n}));
end

%% Byte align the data
tell1 = ftell(dat_fid);
fwrite(dat_fid, zeros((32-mod(tell1,32)),1,'uint8'), 'uint8');
dataFieldLoc = ftell(dat_fid);

%% Write dataFieldLoc in place holder
fseek(dat_fid, 0, 'bof');
fwrite(dat_fid, dataFieldLoc, 'int32');
fseek(dat_fid, dataFieldLoc, 'bof');


%% Write all the data and mdh information

for n=1:length(mdhs),
    writeMdh(dat_fid, mdhs(n));
    temp = [real(rawData{n}); imag(rawData{n})];
    fwrite(dat_fid, temp, 'single');
end


%% Clean up
fclose(dat_fid);

end
