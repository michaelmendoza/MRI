function [ rawData, fileHeaders, mdhs ] = readDatFile( filename )
%function [ rawData, fileHeaders, mdhs ] = readDatFile( filename )
%READDATFILE Reads a Siemens raw data file (.dat) from an MRI scanner
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
if nargin < 1
  [temp path] = uigetfile('*.dat','Select File to Read');
  filename = [path temp];
end

%% using data structure defined on page 150 of the ice manual

dat_fid = fopen(filename,'r','ieee-le');
dataFieldLoc = fread(dat_fid,1,'int32');
numHeaderStructs = fread(dat_fid,1,'int32');

for n = 1:numHeaderStructs,
    [cname cprotocol] = readDatHeaderProtocol(dat_fid);
    fileHeaders.(cname) = cprotocol;
end

%% Byte align the data
tell1 = ftell(dat_fid);
if mod(tell1,32) ~= 0,
    bytealign = fread(dat_fid, (32-mod(tell1,32)), 'uint8');
else
    bytealign = [];
end
tell2 = ftell(dat_fid);

fseek(dat_fid,0,'eof');
file_end = ftell(dat_fid);
fseek(dat_fid,tell2,'bof');

if dataFieldLoc ~= tell2
    error('Headers do not match up to where the data should start');
end


%% Read all the data and mdh information
parameter_prot = 'Meas';
cols = readParameter(fileHeaders.(parameter_prot)', 'iNoOfFourierColumns');
lins = readParameter(fileHeaders.(parameter_prot)', 'iNoOfFourierLines');
pars = readParameter(fileHeaders.(parameter_prot)', 'iNoOfFourierPartitions');
aves = readParameter(fileHeaders.(parameter_prot)', 'lAverages');
reps = readParameter(fileHeaders.(parameter_prot)', 'lRepetitions')+1;
epif = readParameter(fileHeaders.(parameter_prot)', 'lEPIFactor');
coilsinfo = readParameterArray(fileHeaders.('MeasYaps')', 'asList');
coils = length(coilsinfo);
sliceinfo = readParameterArray(fileHeaders.('MeasYaps')', 'asSlice');
slices = length(sliceinfo);

line_counter = 1;
% rawData = zeros(731*2*20,820,'single');
cur_size = lins*pars;
if epif > 1,
    cur_size = cur_size + 3;
end
if reps > 0,
    cur_size = cur_size*reps;
end
cur_size = (cur_size*aves*slices+1)*coils;
rawData = cell(1, cur_size);
fields = {'ulFlagsAndDMALength','lMeasUID','ulScanCounter','ulTimeStamp','ulPMUTimeStamp','aulEvalInfoMask','ushSamplesInScan','ushUsedChannels','sLC','sCutOff','ushKSpaceCentreColumn','ushCoilSelect','fReadOutOffcentre','ulTimeSinceLastRF','ushKSpaceCentreLineNo','ushKSpaceCentrePartitionNo','aushIceProgramPara','aushFreePara','sSD','ushChannelId','ushPTABPosNeg'};
mdhs = cell2struct(cell(cur_size,length(fields)),fields,2);
while(ftell(dat_fid)~=file_end)
    if(line_counter > cur_size)
        rawData(cur_size+1:cur_size*2) = cell(1, cur_size);
        mdhs(cur_size+1:cur_size*2) = cell2struct(cell(cur_size,length(fields)),fields,2);
        cur_size = cur_size*2;
        disp('Initial size calculation wrong, too small');
    end
    mdhs(line_counter) = readMdh(dat_fid);
    if(bitget(mdhs(line_counter).aulEvalInfoMask(1),MDH_ACQEND)),
        disp('MDH_ACQEND found');
%         break;
    end
    temp = fread(dat_fid, [2 double(mdhs(line_counter).ushSamplesInScan)], 'single=>single');
    rawData{line_counter} = complex(temp(1,:),temp(2,:));
    line_counter = line_counter + 1;
end
line_counter = line_counter - 1;
cur_size
line_counter
if(line_counter ~= cur_size)
    rawData = rawData(1:line_counter);
    mdhs = mdhs(1:line_counter);
    disp('Initial size calculation wrong, too big');
end


%% Clean up
fclose(dat_fid);

end


function val = readParameter(protocol, parameter)
%READPARAMETER Reads a given parameter from the text based protocol.

indx = strfind(protocol,parameter) + length(parameter);
valstr = sscanf(protocol(indx:indx+35)','%s');
val = str2num(valstr(isstrprop(valstr,'digit')));

end

function val = readParameterArray(protocol, parameter)
%READPARAMETER Reads a given parameter from the text based protocol.

indxa = strfind(protocol,parameter)+ length(parameter);
cnt = 0;

for n=1:length(indxa),
	indx = indxa(n);
    valindx = sscanf(protocol(indx:indx+35)','[%d]',1)+1;
    indx = indx + length(['[' num2str(valindx) '].']);
    valfieldend = strfind(protocol(indx:indx+35),' ');
    valfield = protocol(indx:(indx+valfieldend(1)-2));
    valfield = strrep(valfield,'.','_');
    indx = indx + valfieldend(1) -2 +3;
    valstr = sscanf(protocol(indx:indx+35),'%s',1);
    val(valindx).(valfield) = str2num(valstr(isstrprop(valstr,'digit')));
end

end
