function [  ] = writeDatHeaderProtocol( dat_fid, name, protocol )
%WRITEDATHEADERPROTOCOL Writes an individual header protocol to a Siemens
%raw data file (.dat) like it was output from an MRI machine.
%   A header protocol in a Siemens raw data file (.dat) contains
%   information about the measurement as a whole.  Each raw data file
%   contains many protocols for various reasons.  It is stored with the
%   name of the protocol (null terminated string), length of the protocol,
%   then text for the given length.  Most protocols are in the Siemens
%   proprietary XProtocol format.
%   File structure defined on page 150 of the ICE manual
%   
%   AUTHOR: Danny Park
%
%   See also readDatFile, readDatHeaderProtocol, readMdh,
%   readProtHeadFromFiles, writeDatFile, writeDatHeaderProtocol, writeMdh,
%   writeProtHead2Files

%write the name as a null terminated string
fwrite(dat_fid, [name char(0)], '*char');

%write place holder for an int (length of the data structure) 
prot_length_pos = ftell(dat_fid);
fwrite(dat_fid,0,'int32');

%read the protocol of length prot_length
prot_length = fprintf(dat_fid, '%c',protocol);

%write the length of the protocol at the place holder
cpos = ftell(dat_fid);
fseek(dat_fid, prot_length_pos, 'bof');
fwrite(dat_fid,prot_length,'int32');
fseek(dat_fid, cpos, 'bof');

end

