function [  ] = anonymizeFromFile( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Constants
MDH_ACQEND = 0 + 1;%0 indexed in C, 1 indexed in MATLAB

%% argument checking
if nargin < 1
  [temp path] = uigetfile('*.dat','Select File to Read');
  filename = [path temp];
end

%% using data structure defined on page 150 of the ice manual

dat_fid = fopen(filename,'r+','ieee-le');
dataFieldLoc = fread(dat_fid,1,'int32');
numHeaderStructs = fread(dat_fid,1,'int32');
curPos = ftell(dat_fid);


%%Regular expressions to delete any kind of format of the parameter
parameters2erase = {'<ParamString."tPatientName">[\r\n \t]*{[^\{\}]*}',...
    '<ParamDouble."flPatientHeight">[\r\n \t]*{[^\{\}]*}',...
    '<ParamString."PatientsName">[\r\n \t]*{[^\{\}]*}',...
    '<ParamString."PatientID">[\r\n \t]*{[^\{\}]*}',...
    '<ParamDouble."flPatientAge">[\r\n \t]*{[^\{\}]*}',...
    '<ParamString."PatientBirthDay">[\r\n \t]*{[^\{\}]*}',...
    '<ParamLong."PatientSex">[\r\n \t]*{[^\{\}]*}',...
    '<ParamLong."lPatientSex">[\r\n \t]*{[^\{\}]*}',...
    '<ParamDouble."flUsedPatientWeight">[\r\n \t]*{[^\{\}]*}',...
    '<ParamString."PatientLOID">[\r\n \t]*{[^\{\}]*}',...
    '<ParamString."Patient">[\r\n \t]*{[^\{\}]*}',...
    '<ParamString."PatientLoid">[\r\n \t]*{[^\{\}]*}',...
    '<ParamString."PatientName">[\r\n \t]*{[^\{\}]*}'};



%% Get all the file as characters and then find the parameters and replace them by an empty string them with the same size

orig_file = char(fread(dat_fid,dataFieldLoc-curPos,'uint8'))';

for i=1:length(parameters2erase)
    patt = parameters2erase{i};
    parameter = regexp(orig_file,patt,'match');
    for j=1:length(parameter)
        orig_file = strrep(orig_file, parameter{j}, char(ones(size(parameter{j}))*32));
    end
end

fseek(dat_fid,curPos,'bof');
fprintf(dat_fid,'%c',orig_file);

end


%%
function val = readParameter(protocol, parameter)
%READPARAMETER Reads a given parameter from the text based protocol.

indx = strfind(protocol,parameter) + length(parameter);
valstr = sscanf(protocol(indx:indx+35)','%s');
val = str2num(valstr(isstrprop(valstr,'digit')));

end

%%
function val = readStringParameter(protocol, parameter)
%READPARAMETER Reads a given parameter from the text based protocol.

indx = strfind(protocol,parameter) + length(parameter) + 1;
valstr = sscanf(protocol(indx:indx+200)','%c');
indx2 = find(valstr == '{'); %find the quotation marks to get the string between
indx3 = find(valstr == '}'); %find the quotation marks to get the string between

if(isempty(indx2)) % in case it is empty and parameter doesn't exist
    val = '';
else
val = valstr((indx2(1)+1): (indx3(1)-1)); %get the text between the indexes
%val = val(isstrprop(val,'alphanum'));
val = strtrim(val);
end


end

%% old Code

% for n = 1:numHeaderStructs,
%     [cname cprotocol] = readDatHeaderProtocol(dat_fid);
%      fileHeaders.(cname) = cprotocol;
%     % create blank string size of cprotocol
%     blank_prot.(cname) = char(ones(size(cprotocol))*32);
%     % find relevant information in cprotocol
%     % index cell with {} to return inner data rather than 1x1 cell
%      name = readStringParameter(fileHeaders.(cname)', strings2keep{1});
%      birthday = readStringParameter(fileHeaders.(cname)',strings2keep{2});
%         
%     % take out all the relevant information and then copy into blank string
%      if(~isempty(name))
%         new_prot = strrep(cprotocol',name,char(ones(size(name))*32));
%      end
%      if(~isempty(birthday))
%          new_prot = strrep(new_prot,birthday,'');
%      end
%      blank_prot.(cname) = new_prot;
%         
% end
% 
% % for in place, go back to beginning of file (leaving position of data)
% % loop over headers writing all blank except relavent info
% headerFields = fieldnames(blank_prot);
% numHeaderStructs = length(headerFields);
% fseek(dat_fid,4,'bof');
% fid_new = fopen(destinationFile,'w','ieee-le');
% 
% for n = 1:numHeaderStructs,
%     writeDatHeaderProtocol(fid_new, headerFields{n}, blank_prot.(headerFields{n})');
% end
% 
% % %% Byte align the data
% tell1 = ftell(dat_fid);
% if mod(tell1,32) ~= 0,
%     bytealign = fread(dat_fid, (32-mod(tell1,32)), 'uint8');
% else
%     bytealign = [];
% end
% tell2 = ftell(dat_fid);
% 
% fseek(dat_fid,0,'eof');
% file_end = ftell(dat_fid);
% fseek(dat_fid,dataFieldLoc,'bof');
% 
% if dataFieldLoc ~= tell2
%      error('Headers do not match up to where the data should start');
% end
% 
% 
% % Read all the data and mdh information
% 
% 
% while(~feof(dat_fid))
%     t = fread(dat_fid , 'int32');
%     fwrite(fid_new,t,'int32');
% end
% 
%  fclose(fid_new);




% Finding information to erase
% for i=1:length(strings2keep)
% parameter2erase = readStringParameter(orig_file, strings2keep{i});
% orig_file = strrep(orig_file, parameter2erase,char(ones(size(parameter2erase))*32));
% end
