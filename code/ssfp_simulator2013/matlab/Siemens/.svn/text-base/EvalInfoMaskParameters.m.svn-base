EvalInfoMask.MDH_ACQEND            = 0 + 1;
EvalInfoMask.MDH_RTFEEDBACK        = 1 + 1;
EvalInfoMask.MDH_HPFEEDBACK        = 2 + 1;
EvalInfoMask.MDH_ONLINE            = 3 + 1;
EvalInfoMask.MDH_OFFLINE           = 4 + 1;
EvalInfoMask.MDH_SYNCDATA          = 5 + 1;       % readout contains synchroneous data
EvalInfoMask.MDH_LASTSCANINCONCAT  = 8 + 1;       % Flag for last scan in concatination

EvalInfoMask.MDH_RAWDATACORRECTION = 10 + 1;      % Correct the rawadata with the rawdata correction factor
EvalInfoMask.MDH_LASTSCANINMEAS    = 11 + 1;      % Flag for last scan in measurement
EvalInfoMask.MDH_SCANSCALEFACTOR   = 12 + 1;      % Flag for scan specific additional scale factor
EvalInfoMask.MDH_2NDHADAMARPULSE   = 13 + 1;      % 2nd RF exitation of HADAMAR
EvalInfoMask.MDH_REFPHASESTABSCAN  = 14 + 1;      % reference phase stabilization scan         
EvalInfoMask.MDH_PHASESTABSCAN     = 15 + 1;      % phase stabilization scan
EvalInfoMask.MDH_D3FFT             = 16 + 1;      % execute 3D FFT         
EvalInfoMask.MDH_SIGNREV           = 17 + 1;      % sign reversal
EvalInfoMask.MDH_PHASEFFT          = 18 + 1;      % execute phase fft     
EvalInfoMask.MDH_SWAPPED           = 19 + 1;      % swapped phase/readout direction
EvalInfoMask.MDH_POSTSHAREDLINE    = 20 + 1;      % shared line               
EvalInfoMask.MDH_PHASCOR           = 21 + 1;      % phase correction data    
EvalInfoMask.MDH_PATREFSCAN        = 22 + 1;      % additonal scan for PAT reference line/partition
EvalInfoMask.MDH_PATREFANDIMASCAN  = 23 + 1;      % additonal scan for PAT reference line/partition that is also used as image scan
EvalInfoMask.MDH_REFLECT           = 24 + 1;      % reflect line              
EvalInfoMask.MDH_NOISEADJSCAN      = 25 + 1;      % noise adjust scan --> Not used in NUM4        
EvalInfoMask.MDH_SHARENOW          = 26 + 1;      % all lines are acquired from the actual and previous e.g. phases
EvalInfoMask.MDH_LASTMEASUREDLINE  = 27 + 1;      % indicates that the current line is the last measured line of all succeeding e.g. phases
EvalInfoMask.MDH_FIRSTSCANINSLICE  = 28 + 1;      % indicates first scan in slice = needed for time stamps + 1
EvalInfoMask.MDH_LASTSCANINSLICE   = 29 + 1;      % indicates  last scan in slice = needed for time stamps + 1
EvalInfoMask.MDH_TREFFECTIVEBEGIN  = 30 + 1;      % indicates the begin time stamp for TReff = triggered measurement + 1
EvalInfoMask.MDH_TREFFECTIVEEND    = 31 + 1;      % indicates the   end time stamp for TReff = triggered measurement + 1
EvalInfoMask.MDH_MDS_REF_POSITION  = 32 + 1;      % indicates the reference position for move during scan images = must be set once per slice/partition in MDS mode + 1
EvalInfoMask.MDH_SLC_AVERAGED      = 33 + 1;      % indicates avveraged slice for slice partial averaging scheme

EvalInfoMask.MDH_FIRST_SCAN_IN_BLADE       = 40 + 1;  % Marks the first line of a blade
EvalInfoMask.MDH_LAST_SCAN_IN_BLADE        = 41 + 1;  % Marks the last line of a blade
EvalInfoMask.MDH_LAST_BLADE_IN_TR          = 42 + 1;  % Set for all lines of the last BLADE in each TR interval

EvalInfoMask.MDH_RETRO_LASTPHASE           = 45 + 1;  % Marks the last phase in a heartbeat
EvalInfoMask.MDH_RETRO_ENDOFMEAS           = 46 + 1;  % Marks an ADC at the end of the measurement
EvalInfoMask.MDH_RETRO_REPEATTHISHEARTBEAT = 47 + 1;  % Repeat the current heartbeat when this bit is found
EvalInfoMask.MDH_RETRO_REPEATPREVHEARTBEAT = 48 + 1;  % Repeat the previous heartbeat when this bit is found
EvalInfoMask.MDH_RETRO_ABORTSCANNOW        = 49 + 1;  % Just abort everything
EvalInfoMask.MDH_RETRO_LASTHEARTBEAT       = 50 + 1;  % This adc is from the last heartbeat = a dummy + 1
EvalInfoMask.MDH_RETRO_DUMMYSCAN           = 51 + 1;  % This adc is just a dummy scan, throw it away
EvalInfoMask.MDH_RETRO_ARRDETDISABLED      = 52 + 1;  % Disable all arrhythmia detection when this bit is found