function [ name, protocol ] = readDatHeaderProtocol( dat_fid )
%READDATHEADERPROTOCOL Reads an individual header protocol from a Siemens
%raw data file (.dat) output from an MRI machine.
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

%read the name as a null terminated string
name = freadstr(dat_fid);

%read an int (length of the data structure) 
prot_length = fread(dat_fid,1,'int32');

%read the protocol of length prot_length
protocol = fread(dat_fid, prot_length, '*char');

end

