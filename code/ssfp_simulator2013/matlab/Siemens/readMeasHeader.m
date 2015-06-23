function [ out ] = readMeasHeader( infile, outfile, writedata, setevalparams, lin_per_part, num_part )
%READMEASHEADER Summary of this function goes here
%   Detailed explanation goes here
%   
%   AUTHOR: Danny Park

if nargin < 6, num_part = 0; end
if nargin < 4, setevalparams = false; lin_per_part = 0; end
if nargin < 3, writedata = false; end

% variables to find and store the index of the yaps struct
yapsname = 'Meas';
yapsindex = 0;

% Name of the parameter that tells the number of data points aquired per
% readout
NxName = 'iNoOfFourierColumns';

% variables for navigating the data
mdhlength = 128; %in bytes
evalinfomaskindex = 21; %21st byte
infomasklength = 4*2;
MDH_ACQEND = 0 + 1;%0 indexed in C, 1 indexed in MATLAB
MDH_ONLINE = 3 + 1;%0 indexed in C, 1 indexed in MATLAB
% MDH_OFFLINE = 4 + 1;%0 indexed in C, 1 indexed in MATLAB -- not sure what
% this does
addinfomask1 = uint32(0);
addinfomask1 = bitset(addinfomask1, MDH_ONLINE);
addinfomask2 = uint32(0);

%% using data structure defined on page 150 of the ice manual

fid = fopen(infile,'r','ieee-le');
out.dataField = fread(fid,1,'int32');
out.numStructs = fread(fid,1,'int32');

for n = 1:out.numStructs,
    %read null terminated string
    curname = char(50);
    nc = 0;
    while true,
        nc = nc + 1;
        curname(nc) = fread(fid,1,'*char');
        if curname(nc) == char(0), break; end
    end
    out.dataStruct(n).name = curname(1:nc-1);
    if strcmp(out.dataStruct(n).name, yapsname),
        yapsindex = n;
    end
%     out.dataStruct(n).name = fscanf(fid, '%c');
%     garbage = fread(fid,1,'*char');
    
    %read an int (length of the data structure) 
    curlength = fread(fid,1,'int32');
    
    %read the data structure of length the last variable
    out.dataStruct(n).data = fread(fid, curlength, '*char');
    
    fidout = fopen([outfile '.' out.dataStruct(n).name],'w+');
    fprintf(fidout,'%c',out.dataStruct(n).data);
    fclose(fidout);
end

% out.padding = fread(fid, 4, 'uint8'); -- isn't really needed?
out.tell1 = ftell(fid);
out.bytealign = fread(fid, (32-mod(out.tell1,32)), 'uint8');
out.tell2 = ftell(fid);
alldata = fread(fid, inf, 'uint8');
fclose(fid);

if writedata,
    fidout = fopen([outfile '.data'],'w+');
    if setevalparams,
        Nxindex = strfind(out.dataStruct(yapsindex).data',NxName) + length(NxName);
        Nxstr = sscanf(out.dataStruct(yapsindex).data(Nxindex:Nxindex+35)','%s');
        Nx = str2num(Nxstr(isstrprop(Nxstr,'digit')));
        
        % Open file and position it correctly for reading
        fid = fopen(infile,'r','ieee-le');
        fseek(fid,out.dataField,'bof');
        
        %loop through for each line
        counter = 0;
        part_counter = 0;
        acq_counter = 0;
        while(~feof(fid))
            %first get through mdh
            writecount = 0;
            %ignore before evalinfomask
%             beforeread = ftell(fid);
%             beforewrite = ftell(fidout);
            ignore = fread(fid, (evalinfomaskindex - 1), 'uint8');
            writecount = writecount + fwrite(fidout, ignore, 'uint8');
            
            evalinfomask1 = fread(fid, 1, 'uint32');
            evalinfomask2 = fread(fid, 1, 'uint32');
            % Don't change the ACQEND MDH Header
            if bitget(evalinfomask1, MDH_ACQEND) == 0,
                evalinfomask1 = bitor(evalinfomask1, addinfomask1);
                evalinfomask2 = bitor(evalinfomask2, addinfomask2);
            end
%             disp(dec2base(evalinfomask1,2,32));
%             disp(num2str(counter));
            writecount = writecount + 4*fwrite(fidout, evalinfomask1, 'uint32');
            writecount = writecount + 4*fwrite(fidout, evalinfomask2, 'uint32');
            ROpoints = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, ROpoints, 'uint16');
            %ignore the number of channels
            ignore = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, ignore, 'uint16');
            %work with the sloopcounter
            curlin = fread(fid, 1, 'uint16');
%             writecount = writecount + 2*fwrite(fidout, curlin, 'uint16');
            %temporary fix
%             writecount = writecount + 2*fwrite(fidout, mod(counter-1,lin_per_part), 'uint16');
            writecount = writecount + 2*fwrite(fidout,counter, 'uint16');
            curacq = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curacq, 'uint16');
            curslc = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curslc, 'uint16');
            curpart = fread(fid, 1, 'uint16');
%             writecount = writecount + 2*fwrite(fidout, curpart, 'uint16');
            %temporary fix
%             writecount = writecount + 2*fwrite(fidout, floor((counter-1)/lin_per_part), 'uint16');
            writecount = writecount + 2*fwrite(fidout, part_counter, 'uint16');
            curecho = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curecho, 'uint16');
            curphs = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curphs, 'uint16');
            currep = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, currep, 'uint16');
            curset = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curset, 'uint16');
            curseg = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curseg, 'uint16');
            curida = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curida, 'uint16');
            curidb = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curidb, 'uint16');
            curidc = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curidc, 'uint16');
            curidd = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curidd, 'uint16');
            curide = fread(fid, 1, 'uint16');
            writecount = writecount + 2*fwrite(fidout, curide, 'uint16');
            %ignore the rest
            ignore = fread(fid, (mdhlength-writecount), 'uint8');
            writecount = writecount + fwrite(fidout, ignore, 'uint8');
%             afterread = ftell(fid);
%             afterwrite = ftell(fidout);
%             disp(['Read: ' num2str(afterread-beforeread) ', Write: ' num2str(afterwrite-beforewrite)]);
            
            %then read the data
            datatemp = fread(fid, Nx*2, 'float');
            fwrite(fidout, datatemp, 'float');
            counter = counter + 1;
            if(counter == lin_per_part)
                counter = 0;
                part_counter = part_counter + 1;
                if(part_counter == num_part)
                    counter = 0;
                    part_counter = 0;
                    acq_counter = acq_counter + 1;
                end
            end
        end
        
        fclose(fid);
    else
        fwrite(fidout, alldata, 'uint8');
    end
    fclose(fidout);
end

end

