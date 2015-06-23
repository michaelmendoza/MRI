function [  ] = writeMdh( dat_fid, mdh )
%READMDH Writes a Siemens measurement data header for raw data aquired with
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

% MDH_NUMBEROFEVALINFOMASK = 2;
% MDH_NUMBEROFICEPROGRAMPARA = 4;
% MDH_FREEHDRPARA = 4;

fwrite(dat_fid, mdh.ulFlagsAndDMALength,        'uint32'); % bit  0..24: DMA length [bytes]
                                                                   % bit     25: pack bit
                                                                   % bit 26..31: pci_rx enable flags                   4 byte
fwrite(dat_fid, mdh.lMeasUID,                   'int32'  ); % measurement user ID                               4
fwrite(dat_fid, mdh.ulScanCounter,              'uint32'); % scan counter [1...]                               4
fwrite(dat_fid, mdh.ulTimeStamp,                'uint32'); % time stamp [2.5 ms ticks since 00:00]             4
fwrite(dat_fid, mdh.ulPMUTimeStamp,             'uint32'); % PMU time stamp [2.5 ms ticks since last trigger]  4
fwrite(dat_fid, mdh.aulEvalInfoMask,            'uint32'); % evaluation info mask field                        8
fwrite(dat_fid, mdh.ushSamplesInScan,           'uint16'); % # of samples acquired in scan                     2
fwrite(dat_fid, mdh.ushUsedChannels,            'uint16'); % # of channels used in scan                        2   =32
writeLoopCounter(dat_fid, mdh.sLC);                                % loop counters                                    28   =60
writeCutOffData(dat_fid, mdh.sCutOff);                             % cut-off values                                    4
fwrite(dat_fid, mdh.ushKSpaceCentreColumn,      'uint16'); % centre of echo                                    2
fwrite(dat_fid, mdh.ushCoilSelect,              'uint16'); % Bit 0..3: CoilSelect                              2
fwrite(dat_fid, mdh.fReadOutOffcentre,          'single'); % ReadOut offcenter value                           4
fwrite(dat_fid, mdh.ulTimeSinceLastRF,          'uint32'); % Sequence time stamp since last RF pulse           4
fwrite(dat_fid, mdh.ushKSpaceCentreLineNo,      'uint16'); % number of K-space centre line                     2
fwrite(dat_fid, mdh.ushKSpaceCentrePartitionNo, 'uint16'); % number of K-space centre partition                2
fwrite(dat_fid, mdh.aushIceProgramPara,         'uint16'); % free parameter for IceProgram                     8   =88
fwrite(dat_fid, mdh.aushFreePara,               'uint16'); % free parameter                          4 * 2 =   8
writeSliceData(dat_fid, mdh.sSD);                                  % Slice Data                                       28   =124
fwrite(dat_fid, mdh.ushChannelId,               'uint16'); % channel Id must be the last parameter             2
fwrite(dat_fid, mdh.ushPTABPosNeg,              'uint16'); % negative,

end

function [  ] = writeLoopCounter( dat_fid, sLoopCounter )

fwrite(dat_fid, sLoopCounter.ushLine,        'uint16'); % line index
fwrite(dat_fid, sLoopCounter.ushAcquisition, 'uint16'); % acquisition index
fwrite(dat_fid, sLoopCounter.ushSlice,       'uint16'); % slice index
fwrite(dat_fid, sLoopCounter.ushPartition,   'uint16'); % partition index
fwrite(dat_fid, sLoopCounter.ushEcho,        'uint16'); % echo index
fwrite(dat_fid, sLoopCounter.ushPhase,       'uint16'); % phase index
fwrite(dat_fid, sLoopCounter.ushRepetition,  'uint16'); % measurement repeat index
fwrite(dat_fid, sLoopCounter.ushSet,         'uint16'); % set index
fwrite(dat_fid, sLoopCounter.ushSeg,         'uint16'); % segment index  (for TSE)
fwrite(dat_fid, sLoopCounter.ushIda,         'uint16'); % IceDimension a index
fwrite(dat_fid, sLoopCounter.ushIdb,         'uint16'); % IceDimension b index
fwrite(dat_fid, sLoopCounter.ushIdc,         'uint16'); % IceDimension c index
fwrite(dat_fid, sLoopCounter.ushIdd,         'uint16'); % IceDimension d index
fwrite(dat_fid, sLoopCounter.ushIde,         'uint16'); % IceDimension e index

end

function [  ] = writeCutOffData( dat_fid, sCutOffData )

fwrite(dat_fid, sCutOffData.ushPre,  'uint16'); % write ushPre zeros at line start
fwrite(dat_fid, sCutOffData.ushPost, 'uint16'); % write ushPost zeros at line end

end

function [  ] = writeSliceData( dat_fid, sSliceData )

writeVector( dat_fid, sSliceData.sVector );                  % slice position vector
fwrite(dat_fid, sSliceData.aflQuaternion, 'single'); % rotation matrix as quaternion

end

function [  ] = writeVector( dat_fid, sVector )

fwrite(dat_fid, sVector.flSag, 'single');
fwrite(dat_fid, sVector.flCor, 'single');
fwrite(dat_fid, sVector.flTra, 'single');

end

