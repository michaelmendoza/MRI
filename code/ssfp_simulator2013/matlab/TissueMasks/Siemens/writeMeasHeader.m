function [] = writeMeasHeader( writefile, outstruct, datafile )
%WRITEMEASHEADER Summary of this function goes here
%   Detailed explanation goes here
%   
%   AUTHOR: Danny Park


%% using data structure defined on page 150 of the ice manual

fid = fopen(writefile,'w+','ieee-le');
fwrite(fid,0,'int32');
fwrite(fid,length(outstruct.dataStruct),'int32');

for n = 1:length(outstruct.dataStruct),
    %write null terminated string
    fwrite(fid, outstruct.dataStruct(n).name, '*char');
    fwrite(fid, char(0), '*char');
    
    %write an int (length of the data structure) 
    fwrite(fid,length(outstruct.dataStruct(n).data),'int32');
    
    %write the data structure of length the last variable
    fwrite(fid, outstruct.dataStruct(n).data, '*char');
end

%32 bit padding -- isn't really needed?
% fwrite(fid, outstruct.padding, 'uint8');
temptell = ftell(fid);

%pad to 32 byte boundary
if mod(temptell, 32) ~= 0,
    fwrite(fid, zeros(32-mod(temptell, 32),1,'uint8'),'uint8');
    temptell = ftell(fid);
end

%write the length to the beginning of the file
fseek(fid,0,'bof');
fwrite(fid,temptell,'int32');
fseek(fid,temptell,'bof');

%read data from file
dfid = fopen(datafile, 'r');
alldata = fread(dfid, inf, 'uint8');
fclose(dfid);

%write data to a new file
fwrite(fid,alldata,'uint8');
fclose(fid);

end

