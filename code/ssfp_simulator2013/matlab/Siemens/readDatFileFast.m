function [ rawData ] = readDatFileFast( filename )
%function [ rawData ] = readDatFileFast( filename )
%READDATFILE Reads a Siemens raw data file (.dat) from an MRI scanner
%quickly without extra information.
%   This returns a matrix rather than a cell structure with each line.
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

clear fileHeaders;

line_counter = 1;
% rawData = zeros(731*2*20,820,'single');
cur_size = lins*pars;
if epif > 1,
    cur_size = cur_size + 3;
end
if reps > 0,
    cur_size = cur_size*reps;
end
% cur_size = (cur_size*aves*slices+1)*coils; % the +1 accounts for ACQEND signals which are ignored in this implementation
cur_size = (cur_size*aves*slices)*coils;
rawData = complex(zeros(cols, cur_size/coils,coils,'single'));%complex(zeros(cols, lins, pars, aves, reps, coils,'single'));
% fields = {'ulFlagsAndDMALength','lMeasUID','ulScanCounter','ulTimeStamp','ulPMUTimeStamp','aulEvalInfoMask','ushSamplesInScan','ushUsedChannels','sLC','sCutOff','ushKSpaceCentreColumn','ushCoilSelect','fReadOutOffcentre','ulTimeSinceLastRF','ushKSpaceCentreLineNo','ushKSpaceCentrePartitionNo','aushIceProgramPara','aushFreePara','sSD','ushChannelId','ushPTABPosNeg'};
while(ftell(dat_fid)~=file_end)
    mdh = readMdh(dat_fid);
    temp = fread(dat_fid, [2 double(mdh.ushSamplesInScan)], 'single=>single');
    if(bitget(mdh.aulEvalInfoMask(1),MDH_ACQEND)),
        disp('MDH_ACQEND found');
        continue;
    end
    if(line_counter > cur_size)
%         rawData(cur_size+1:cur_size*2) = cell(1, cur_size);
%         mdhs(cur_size+1:cur_size*2) = cell2struct(cell(cur_size,length(fields)),fields,2);
%         cur_size = cur_size*2;
        warning('Initial size calculation wrong, too small.  There are more lines than expected.  This may mean some data is lost.');
    end
    if mdh.ushSamplesInScan == cols,
%         temp = fread(dat_fid, [2 double(mdh.ushSamplesInScan)], 'single=>single');
        if mdh.sLC.ushLine+1 > lins, 
            error('too many lines');
        end
%         rawData(:,mdh.sLC.ushLine+1,mdh.sLC.ushPartition+1,mdh.sLC.ushAcquisition+1,mdh.sLC.ushRepetition+1,mdh.ushChannelId+1) = complex(temp(1,:),temp(2,:));
        rawData(:,cur_size,mdh.ushChannelId+1) = complex(temp(1,:),temp(2,:));
        if mod(line_counter,10000) == 0,
            fprintf(1,'line=%d, part=%d, line_counter=%d\n', mdh.sLC.ushLine,mdh.sLC.ushPartition,line_counter);
        end
    else
        warning(['The number of samples in scan(' num2str(mdh.ushSamplesInScan) ') does not match what is expected(' num2str(cols) ').']);
    end
    line_counter = line_counter + 1;
end
line_counter = line_counter - 1;
if(line_counter ~= cur_size)
    warning('Initial size calculation wrong, too big');
end


%% Clean up
fclose(dat_fid);

end


function val = readParameter(protocol, parameter)
%READPARAMETER Reads a given parameter from the text based protocol.

indx = strfind(protocol,parameter) + length(parameter);
valstr = sscanf(protocol(indx:indx+35)','%s');
val = str2num(valstr(isstrprop(valstr,'digit')));
if isempty(val),
    val = 0;
end

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
