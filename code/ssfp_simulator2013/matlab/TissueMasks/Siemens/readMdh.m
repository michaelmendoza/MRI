function [ mdh ] = readMdh( dat_fid )
%READMDH Reads a Siemens measurement data header for raw data aquired with
%an MRI
%   Raw data files from a Siemens MRI scanner contain a structure called a
%   "measurement data header".  This mdh contains information about a
%   single aquisitions line (which would be a line in kspace in a normal
%   cartesian grid).
%   
%   Structure of mdh found in
%   C:\MIDEA\N4_VB15A_LATEST_20070519\n4\pkg\MrServers\MrMeasSrv\SeqIF\MDH\
%   mdh.h on an IDEA development machine
%   
%   AUTHOR: Danny Park
%
%   See also readDatFile, readDatHeaderProtocol, readMdh,
%   readProtHeadFromFiles, writeDatFile, writeDatHeaderProtocol, writeMdh,
%   writeProtHead2Files

MDH_NUMBEROFEVALINFOMASK = 2;
MDH_NUMBEROFICEPROGRAMPARA = 4;
MDH_FREEHDRPARA = 4;

mdh = struct();
mdh.ulFlagsAndDMALength        = fread(dat_fid,                          1, 'uint32=>uint32'); % bit  0..24: DMA length [bytes]
                                                                                               % bit     25: pack bit
                                                                                               % bit 26..31: pci_rx enable flags                   4 byte
mdh.lMeasUID                   = fread(dat_fid,                          1, 'int32=>int32'  ); % measurement user ID                               4
mdh.ulScanCounter              = fread(dat_fid,                          1, 'uint32=>uint32'); % scan counter [1...]                               4
mdh.ulTimeStamp                = fread(dat_fid,                          1, 'uint32=>uint32'); % time stamp [2.5 ms ticks since 00:00]             4
mdh.ulPMUTimeStamp             = fread(dat_fid,                          1, 'uint32=>uint32'); % PMU time stamp [2.5 ms ticks since last trigger]  4
mdh.aulEvalInfoMask            = fread(dat_fid,   MDH_NUMBEROFEVALINFOMASK, 'uint32=>uint32'); % evaluation info mask field                        8
mdh.ushSamplesInScan           = fread(dat_fid,                          1, 'uint16=>uint16'); % # of samples acquired in scan                     2
mdh.ushUsedChannels            = fread(dat_fid,                          1, 'uint16=>uint16'); % # of channels used in scan                        2   =32
mdh.sLC                        = readLoopCounter(dat_fid);                                     % loop counters                                    28   =60
mdh.sCutOff                    = readCutOffData(dat_fid);                                      % cut-off values                                    4
mdh.ushKSpaceCentreColumn      = fread(dat_fid,                          1, 'uint16=>uint16'); % centre of echo                                    2
mdh.ushCoilSelect              = fread(dat_fid,                          1, 'uint16=>uint16'); % Bit 0..3: CoilSelect                              2
mdh.fReadOutOffcentre          = fread(dat_fid,                          1, 'single=>single'); % ReadOut offcenter value                           4
mdh.ulTimeSinceLastRF          = fread(dat_fid,                          1, 'uint32=>uint32'); % Sequence time stamp since last RF pulse           4
mdh.ushKSpaceCentreLineNo      = fread(dat_fid,                          1, 'uint16=>uint16'); % number of K-space centre line                     2
mdh.ushKSpaceCentrePartitionNo = fread(dat_fid,                          1, 'uint16=>uint16'); % number of K-space centre partition                2
mdh.aushIceProgramPara         = fread(dat_fid, MDH_NUMBEROFICEPROGRAMPARA, 'uint16=>uint16'); % free parameter for IceProgram                     8   =88
mdh.aushFreePara               = fread(dat_fid,            MDH_FREEHDRPARA, 'uint16=>uint16'); % free parameter                          4 * 2 =   8
mdh.sSD                        = readSliceData(dat_fid);                                       % Slice Data                                       28   =124
mdh.ushChannelId               = fread(dat_fid,                          1, 'uint16=>uint16'); % channel Id must be the last parameter             2
mdh.ushPTABPosNeg              = fread(dat_fid,                          1, 'uint16=>uint16'); % negative,

end

function [ sLoopCounter ] = readLoopCounter( dat_fid )

sLoopCounter = struct();
sLoopCounter.ushLine        = fread(dat_fid, 1, 'uint16=>uint16'); % line index
sLoopCounter.ushAcquisition = fread(dat_fid, 1, 'uint16=>uint16'); % acquisition index
sLoopCounter.ushSlice       = fread(dat_fid, 1, 'uint16=>uint16'); % slice index
sLoopCounter.ushPartition   = fread(dat_fid, 1, 'uint16=>uint16'); % partition index
sLoopCounter.ushEcho        = fread(dat_fid, 1, 'uint16=>uint16'); % echo index
sLoopCounter.ushPhase       = fread(dat_fid, 1, 'uint16=>uint16'); % phase index
sLoopCounter.ushRepetition  = fread(dat_fid, 1, 'uint16=>uint16'); % measurement repeat index
sLoopCounter.ushSet         = fread(dat_fid, 1, 'uint16=>uint16'); % set index
sLoopCounter.ushSeg         = fread(dat_fid, 1, 'uint16=>uint16'); % segment index  (for TSE)
sLoopCounter.ushIda         = fread(dat_fid, 1, 'uint16=>uint16'); % IceDimension a index
sLoopCounter.ushIdb         = fread(dat_fid, 1, 'uint16=>uint16'); % IceDimension b index
sLoopCounter.ushIdc         = fread(dat_fid, 1, 'uint16=>uint16'); % IceDimension c index
sLoopCounter.ushIdd         = fread(dat_fid, 1, 'uint16=>uint16'); % IceDimension d index
sLoopCounter.ushIde         = fread(dat_fid, 1, 'uint16=>uint16'); % IceDimension e index

end

function [ sCutOffData ] = readCutOffData( dat_fid )

sCutOffData = struct();
sCutOffData.ushPre  = fread(dat_fid, 1, 'uint16=>uint16'); % write ushPre zeros at line start
sCutOffData.ushPost = fread(dat_fid, 1, 'uint16=>uint16'); % write ushPost zeros at line end

end

function [ sSliceData ] = readSliceData( dat_fid )

sSliceData = struct();
sSliceData.sVector = readVector( dat_fid );                     % slice position vector
sSliceData.aflQuaternion = fread(dat_fid, 4, 'single=>single'); % rotation matrix as quaternion

end

function [ sVector ] = readVector( dat_fid )

sVector = struct();
sVector.flSag = fread(dat_fid, 1, 'single=>single');
sVector.flCor = fread(dat_fid, 1, 'single=>single');
sVector.flTra = fread(dat_fid, 1, 'single=>single');

end

