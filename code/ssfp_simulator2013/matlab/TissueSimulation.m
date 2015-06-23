function TissueSimulation
    %clear; clf; close all;
    
    if 0
        GenerateSpectraPlots()
%         GeneratePhaseDiffPlots()
        
    else
        
        clear; clf; close all;
        
        if 0
            GenerateNewData()
            %SimpleSeparation();
        else
            %GenerateNewDataForDixonSeparation()
            [S, S2, S3] = LoadDixonData()
            %DixonFatWaterSeparation(S, S2, S3);
        end
        
    end
    
end

function [S, S2, S3] = LoadDixonData()
    %load('DixonWithNoise.mat')
    %load('DixonWithNoisewithGradients.mat')
    %load('Data/SSPFsignalwithdPhiNoGradMod2')
    %load('6SSPFsignalNoNoise.mat')
    load('6SSPFsignalWithNoise.mat');
end

function SimpleSeparation

    % Load data
    %load('SSPFsignals.mat')
    %load('SSPFsignalswithGrad')
    %load('SSPFsignalwithdPhi')
    %load('SSPFsignalwithdPhiNoGrad')
    %load('SSPFsignalwithdPhiNoGradMod')
    load('SSPFsignalwithdPhiNoGradMod2')
    %load('6SSPFsignalwithdPhiNoGradMod.mat');

    % Get tissue data
    tissues = GenerateTissueProperties(); 
    p = [Pm2, P, P2, P3];
    s = {Sm2, S, S2, S3};
    for n = 1:length(p)
        Mxy = GenerateMxySignal(tissues.WhiteMatter, p(n));
        Mc = Mxy(1) + 1i * Mxy(2); 
        Mtheta(n) = angle(Mc);
    end
    
    % Calculate Constant Phase Difference for WhiteMatter
    startIndex = 3;
    for n = 1:length(p)
        dTheta(n) = Mtheta(n) - Mtheta(startIndex);
    end
    
    % Adjust Signal Phase
    for n = 1:length(p) 
        s{n} = s{n} * exp(-1i*dTheta(n));
    end
    
    % Select Top Signal Magnitudes
    PlotMagnitiudePhaseFor4Images(s{1}, s{2}, s{3}, s{4});
    [s, indices] = SelectTop3Pixels(s);
    ps = CalculatePhaseState(indices);
    PlotMagnitiudePhaseFor3Images(s{1}, s{2}, s{3});
    
    % Normalize Signal Magnitude
    for n = 1:3
        mag{n} = abs(s{n});
        phase{n} = angle(s{n});
    end
    for r = 1:P.Nx
        for c = 1:P.Ny
            norm = mag{1}(r,c) * mag{2}(r,c) * mag{3}(r,c) / 3;
            s{1}(r,c) = norm * exp(1i*phase{1}(r,c));
            s{2}(r,c) = norm * exp(1i*phase{2}(r,c));
            s{3}(r,c) = norm * exp(1i*phase{3}(r,c));
        end
    end

    PlotMagnitiudePhaseFor3Images(s{1}, s{2}, s{3});

    % Decompose Fat/Water Images - Using 2-Point Dixon
    [Water, Fat] = TwoPointDixon(S, S2);
    [AlgoWater, AlgoFat] = TwoPointDixonWithSort(s{1}, s{2}, s{3}, ps{1}, ps{2}, ps{3});
    figure();
    subplot(2,2,1);
    imshow(abs(Water),[]);
    title('Water Image');
    subplot(2,2,2);
    imshow(abs(Fat),[]);
    title('Fat Image');
    subplot(2,2,3);
    imshow(abs(AlgoWater),[]);
    title('Water Image');
    subplot(2,2,4);
    imshow(abs(AlgoFat),[]);
    title('Fat Image');

    % Decompose Fat/Water Images - Using 3-Point Dixon
    [Water2, Fat2] = ThreePointDixon(Sm2, S, S2);
    [AlgoWater2, AlgoFat2] = ThreePointDixonWithSort(s{1}, s{2}, s{3}, ps{1}, ps{2}, ps{3});
    figure();
    subplot(2,2,1);
    imshow(abs(Water2),[]);
    title('Water Image');
    subplot(2,2,2);
    imshow(abs(Fat2),[]);
    title('Fat Image');
    subplot(2,2,3);
    imshow(abs(AlgoWater2),[]);
    title('Water Image');
    subplot(2,2,4);
    imshow(abs(AlgoFat2),[]);
    title('Fat Image');
    
%     %Calcalute SNR and CNR
%     SNRW = CalculateSNR(Water, tissue, tissues.WhiteMatter);
%     SNRF = CalculateSNR(Fat, tissue, tissues.Fat);
%     SNRW2 = CalculateSNR(Water2, tissue, tissues.WhiteMatter);
%     SNRF2 = CalculateSNR(Fat2, tissue, tissues.Fat);
%     [SNRW SNRF SNRW2 SNRF2]
%     
%     CNRW = CalculateCNR(Water, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRF = CalculateCNR(Fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRW2 = CalculateCNR(Water2, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRF2 = CalculateCNR(Fat2, tissue, tissues.WhiteMatter, tissues.Fat);
%     [CNRW CNRF CNRW2 CNRF2]
        
end

function DixonFatWaterSeparation(S, S2, S3)
  
    % Select Top Signal Magnitudes
    PlotMagnitudePhaseImage(S);
    PlotMagnitudePhaseImage(S2);
    PlotMagnitiudePhaseFor3Images(S, S2, S3);

    % Decompose Fat/Water Images - Using 2-Point Dixon
    [Water, Fat] = TwoPointDixon(S, S2);
    figure();
    subplot(1,2,1);
    imshow(abs(Water),[]);
    title('Water Image');
    subplot(1,2,2);
    imshow(abs(Fat),[]);
    title('Fat Image');

    % Decompose Fat/Water Images - Using 3-Point Dixon
    [Water2, Fat2] = ThreePointDixon(S, S2, S3);
    figure();
    subplot(1,2,1);
    imshow(abs(Water2),[]);
    title('Water Image');
    subplot(1,2,2);
    imshow(abs(Fat2),[]);
    title('Fat Image');
        
end

function [SelectPixels, SelectIndices] = SelectTop3Pixels(Images)

    % Initialize
    N = length(Images);
    [Ny, Nx] = size(Images{1});
    SelectPixels  = {zeros(Ny,Nx), zeros(Ny,Nx), zeros(Ny,Nx)};
    SelectIndices = {zeros(Ny,Nx), zeros(Ny,Nx), zeros(Ny,Nx)};
    
    % Select top 3 pixels from Phase Cycled Images
    for r = 1:Ny
        for c = 1:Nx
            
            % Select top 3 pixels
            signal = zeros(1, N);
            for n = 1:length(Images)
                signal(n) = Images{n}(r,c);
            end
            
            [mag, indices] = InsertionSort(abs(signal));

            % Save results
            for n = 1:3
                SelectPixels{n}(r,c) = signal(indices(N+1-n));
                SelectIndices{n}(r,c) = indices(N+1-n);
            end
            
        end
    end
end

function [PhaseState] = CalculatePhaseState(Indices)
    N = length(Indices);
    [Ny, Nx] = size(Indices{1});
    for n = 1:N
        for r = 1:Ny
            for c =1:Nx
                Indices{n}(r,c) = mod(Indices{n}(r,c), 2);
            end
        end
    end
    PhaseState = Indices;
end

function [Array, Indices] = InsertionSort(Array)
 
    Indices = 1:numel(Array);
    for n = (2:numel(Array))
 
        value = Array(n);
        index = Indices(n);
        m = n - 1;
 
        while (m >= 1) && (Array(m) > value)
            Array(m+1) = Array(m);
            Indices(m+1) = Indices(m);
            m = m-1;
        end
 
        Array(m+1) = value;
        Indices(m+1) = index;
    end 
end 

function GenerateSpectraPlots
   
    tissues = GenerateTissueProperties(); 

    % Generate Pulse Parameters
    P = struct('dt', 10e-6, ...                     % dt - Sample Period
               'TR', 10/1000, ...                   % TR - Time between Flips
               'TE', 5/ 1000, ...                   % TE - Time from Flip to Sample
               'Nr', 200, ...                       % Number of Repetitions
               'Alpha', pi/8, ...                   % Tip Angle
               'Phi', 0, ...                        % Phi
               'dPhi', pi);                         % dPhi
    P2 = P;
    P3 = P;
    P4 = P;
    Pm2 = P;
    Pm3 = P;
    P2.TE = P.TE + abs(1/(2*tissues.Fat.Fo));
    P3.TE = P.TE + 2 * abs(1/(2*tissues.Fat.Fo));
    P4.TE = P.TE + 3 * abs(1/(2*tissues.Fat.Fo));
    Pm2.TE = P.TE - abs(1/(2*tissues.Fat.Fo));
    Pm3.TE = P.TE - 2 * abs(1/(2*tissues.Fat.Fo));
    Pm3.dPhi = Pm3.dPhi - pi/2;
    Pm2.dPhi = Pm2.dPhi - pi/4;
    P2.dPhi = P2.dPhi + pi/4;
    P3.dPhi = P3.dPhi + pi/2;
    P4.dPhi = P4.dPhi + 3*pi/4;
    
    % View Off-Resonance Spectra for Tissues
%     PlotOffResonanceSpectrum([tissues.WhiteMatter, tissues.Fat], Pm3, 100);
%     PlotOffResonanceSpectrum([tissues.WhiteMatter, tissues.Fat], Pm2,  100);
    PlotOffResonanceSpectrum([tissues.WhiteMatter, tissues.Fat], P,  100);
%     PlotOffResonanceSpectrum([tissues.WhiteMatter, tissues.Fat], P2, 100);
%     PlotOffResonanceSpectrum([tissues.WhiteMatter, tissues.Fat], P3, 100);
%     PlotOffResonanceSpectrum([tissues.WhiteMatter, tissues.Fat], P4, 100);
    
end

function GeneratePhaseDiffPlots

    % Generate Tissue Types
    tissues = GenerateTissueProperties();

    % Generate Pulse Parameters
    P = struct('dt', 10e-6, ...                     % dt - Sample Period
               'TR', 10/1000, ...                   % TR - Time between Flips
               'TE', 5/ 1000, ...                   % TE - Time from Flip to Sample
               'Nr', 200, ...                       % Number of Repetitions
               'Alpha', pi/8, ...                   % Tip Angle
               'Phi', 0, ...                        % Phi
               'dPhi', pi);                         % dPhi
     
    PlotdPhiSpectrum([tissues.WhiteMatter, tissues.Fat], P, 3/2*pi);
           
    % Generate Spectrum at beta = 0
    % Initalize
    N = 100;
    Sw = zeros(N); Sf = zeros(N); dS = zeros(N);             % Initalize Signals Vectors
    TE = linspace(-abs(1/(2*tissues.Fat.Fo)), abs(1/(2*tissues.Fat.Fo)), N);
    TE = TE + P.TE;
    for n = 1:N
        Pn = P;
        Pn.TE = TE(n);
        Mw = GenerateMxySignal(tissues.WhiteMatter, Pn);     % Get Mxy Signal
        Mf = GenerateMxySignal(tissues.Fat, Pn);             % Get Mxy Signal
        Sw(n) = Mw(1) + 1i * Mw(2);
        Sf(n) = Mf(1) + 1i * Mf(2);
    end
    
    TE = TE - P.TE;
    figure()
    subplot(2,1,1); plot(TE, abs(Sw), TE, abs(Sf));
    subplot(2,1,2); plot(TE, angle(Sw), TE, angle(Sf), TE, angle(Sf)-angle(Sw));
    
    
    % Generate Spectrum at beta = 0
    N = 201;
    Sw = zeros(N); Sf = zeros(N);                            % Initalize Signals Vectors
    dPhi = linspace(-pi, pi, N);
    for n = 1:N
        Pn = P;
        Pn.dPhi = P.dPhi + dPhi(n);
        Mw = GenerateMxySignal(tissues.WhiteMatter, Pn);     % Get Mxy Signal
        Mf = GenerateMxySignal(tissues.Fat, Pn);             % Get Mxy Signal
        Sw(n) = Mw(1) + 1i * Mw(2);
        Sf(n) = Mf(1) + 1i * Mf(2);
    end
    
    figure();
    subplot(2,1,1); plot(dPhi, abs(Sw), dPhi, abs(Sf));
    subplot(2,1,2); plot(dPhi, angle(Sw), dPhi, angle(Sf));

end

function GenerateNewData
    % Set Flags
    UseFoPSF         = 1;       % Use Off-Resonance Point Spread Function (PSF) to create Gaussian Fo Distributions
    UseRGradient     = 1;       % Use Linear Off-Resonance Gradient
    FoSD             = 10;       % Off-Resonance Std. Dev. in Hz
    MaxGradFreqz     = 100;     % Max Frequency for Off-Resonance Linear Gradient
    NoiseSD          = 0.01;    % Signal Noise Std. Dev.
    
    % Generate Tissue Model
    filenames = {'TissueMasks/RectPhantom64Tissues.png', 'TissueMasks/Phantom64Tissues.png', 'TissueMasks/Phantom256Tissues.png'};
    filename = filenames{2};
    tissues = GenerateTissueProperties();
    tissue = GenerateTwoTissueModel(filename, tissues);   

    % Generate Off-Resonance PSF
    if UseFoPSF
       tissue = GenerateResonancePSF(tissue, FoSD); 
    end
    
    % Generate Off-Resonance Gradient
    if UseRGradient
        tissue = GenerateResonanceGradient(tissue, MaxGradFreqz);
    end

    % Generate Pulse Parameters
    P = struct('Nx', length(tissue(1,:)), ...       % Number of Samples in x
               'Ny', length(tissue(:,1)), ...       % Number of Samples in y
               'FOVx', 0.256, 'FOVy', 0.256, ...    % Field of View in meters
               'dt', 10e-6, ...                     % dt - Sample Period
               'tau', 1e-3, ...                     % tau - Time for Phase Encode  
               'B0', 3, ...                         % B0 Field (T)
               'Gamma', 3 * 42.58e6, ...            % Gyromagnetic Ratio (Hz/T) 
               'TR', 10/1000, ...                   % TR - Time between Flips
               'TE', 5/ 1000, ...                   % TE - Time from Flip to Sample
               'Nr', 200, ...                       % Number of Repetitions
               'Alpha', pi/8, ...                   % Tip Angle
               'Phi', 0, ...                        % Phi
               'dPhi', pi);                          % dPhi
    P = GenerateResolutionAndGradientAmplitudes(P);
    P2 = P;
    P3 = P;
    P4 = P;
    Pm2 = P;
    Pm3 = P;
    P2.TE = P.TE + abs(1/(2*tissues.Fat.Fo));
    P3.TE = P.TE + 2 * abs(1/(2*tissues.Fat.Fo));
    P4.TE = P.TE + 3 * abs(1/(2*tissues.Fat.Fo));
    Pm2.TE = P.TE - abs(1/(2*tissues.Fat.Fo)) - (.15 *(1.571-1.405) / (2*pi) * P.TR);
    Pm3.TE = P.TE - 2 * abs(1/(2*tissues.Fat.Fo));
    Pm3.dPhi = Pm3.dPhi - pi/2;
    Pm2.dPhi = Pm2.dPhi - pi/4;
    P2.dPhi = P2.dPhi + pi/4;
    P3.dPhi = P3.dPhi + pi/2;
    P4.dPhi = P4.dPhi + 3*pi/4;
    
%     tissueSingals = zeros(3,1);
%     tissueSingals(1) = GenerateMxySignal(tissues.Empty, P);
%     tissueSingals(2) = GenerateMxySignal(tissues.WhiteMatter, P);
%     tissueSingals(3) = GenerateMxySignal(tissues.Fat, P);
  
    % Generate MRI Singal
%     S  = GenerateSignalFromMxyImage(tissue, P);   % Generate MRI Singal
%     S2 = GenerateSignalFromMxyImage(tissue, P2);  % Generate 2nd MRI Signal  @ pi
%     S3 = GenerateSignalFromMxyImage(tissue, P3);  % Generate 3rd MRI Signal @ 2pi
    S = GenerateMxyImage(tissue, P);                % Generate MRI Singal
    S2 = GenerateMxyImage(tissue, P2);              % Generate 2nd MRI Signal  @ pi
    S3 = GenerateMxyImage(tissue, P3);              % Generate 3rd MRI Signal @ 2pi
    S4 = GenerateMxyImage(tissue, P4);
    Sm2 = GenerateMxyImage(tissue, Pm2);            % Generate 2nd MRI Signal  @ -pi
    Sm3 = GenerateMxyImage(tissue, Pm3);            % Generate 3rd MRI Signal @ -2pi

    % Add Complex Signal Noise
    S = AddComplexSignalNoise(S, 0, NoiseSD);
    S2 = AddComplexSignalNoise(S2, 0, NoiseSD);
    S3 = AddComplexSignalNoise(S3, 0, NoiseSD);

    savefile = '6SSPFsignalWithNoise.mat';
    save(savefile, 'tissue', 'S', 'S2', 'S3', 'S4', 'Sm2', 'Sm3', 'P', 'P2', 'P3', 'P4', 'Pm2', 'Pm3')
    
    % Generate Image From MRI Signal
%     S = ifft2(fftshift(S));
%     S2 = ifft2(fftshift(S2)); 
%     S3 = ifft2(fftshift(S3));
    
end

function GenerateNewDataForDixonSeparation
    % Set Flags
    UseFoPSF         = 0;       % Use Off-Resonance Point Spread Function (PSF) to create Gaussian Fo Distributions
    UseRGradient     = 1;       % Use Linear Off-Resonance Gradient
    FoSD             = 20;       % Off-Resonance Std. Dev. in Hz
    MaxGradFreqz     = 200;      % Max Frequency for Off-Resonance Linear Gradient
    NoiseSD          = 0.01;     % Signal Noise Std. Dev.
    
    % Generate Tissue Model
    filenames = {'RectPhantom64Tissues.png', 'Phantom64Tissues.png', 'Phantom256Tissues.png'};
    filename = filenames{2};
    tissues = GenerateTissueProperties();
    tissue = GenerateTwoTissueModel(filename, tissues);   

    % Generate Off-Resonance PSF
    if UseFoPSF
       tissue = GenerateResonancePSF(tissue, FoSD); 
    end
    
    % Generate Off-Resonance Gradient
    if UseRGradient
        tissue = GenerateResonanceGradient(tissue, MaxGradFreqz);
    end

    % Generate Pulse Parameters
    P = struct('Nx', length(tissue(1,:)), ...       % Number of Samples in x
               'Ny', length(tissue(:,1)), ...       % Number of Samples in y
               'FOVx', 0.256, 'FOVy', 0.256, ...    % Field of View in meters
               'dt', 10e-6, ...                     % dt - Sample Period
               'tau', 1e-3, ...                     % tau - Time for Phase Encode  
               'B0', 3, ...                         % B0 Field (T)
               'Gamma', 3 * 42.58e6, ...            % Gyromagnetic Ratio (Hz/T) 
               'TR', 10/1000, ...                   % TR - Time between Flips
               'TE', 5/ 1000, ...                   % TE - Time from Flip to Sample
               'Nr', 200, ...                       % Number of Repetitions
               'Alpha', pi/8, ...                   % Tip Angle
               'Phi', 0, ...                        % Phi
               'dPhi', pi);                          % dPhi
    P = GenerateResolutionAndGradientAmplitudes(P);
    P2 = P;
    P3 = P;
    P2.TE = P.TE + abs(1/(2*tissues.Fat.Fo));
    P3.TE = P.TE + 2 * abs(1/(2*tissues.Fat.Fo));
  
    % Generate MRI Singal
    S = GenerateMxyImage(tissue, P);                % Generate MRI Singal
    S2 = GenerateMxyImage(tissue, P2);              % Generate 2nd MRI Signal  @ pi
    S3 = GenerateMxyImage(tissue, P3);              % Generate 3rd MRI Signal @ 2pi

    % Add Complex Signal Noise
    S = AddComplexSignalNoise(S, 0, NoiseSD);
    S2 = AddComplexSignalNoise(S2, 0, NoiseSD);
    S3 = AddComplexSignalNoise(S3, 0, NoiseSD);

    savefile = 'DixonWithNoisewithGradients.mat';
    save(savefile, 'tissue', 'S', 'S2', 'S3', 'P', 'P2', 'P3')
end

function [TissueProperties] = GenerateTissueProperties()
    % Tissues (T1, T2, ProtonDensity, Off Resonance Frequency (Hz))
    Empty       = struct('Name', 'Empty', ...
                         'T1', 0, ...            
                         'T2', 0, ...            
                         'ProtonDensity', 0, ... 
                         'Fo', 0);             
    WhiteMatter = struct('Name', 'WhiteMatter', ...
                         'T1', 790/1000, ...         
                         'T2', 92/1000, ...          
                         'ProtonDensity', 1, ... 
                         'Fo', 0);             
    Fat         = struct('Name', 'Fat', ...
                         'T1', 270/1000, ...          
                         'T2', 85/1000, ...            
                         'ProtonDensity', 1, ... 
                         'Fo', -428);               % -428Hz @ 3 T    
                     
    TissueProperties = struct('Empty', Empty, 'WhiteMatter', WhiteMatter, 'Fat', Fat');                 
end

function [Tissue] = GenerateTwoTissueModel(Filename, Tissues)

    % Tissues (T1, T2, ProtonDensity, Off Resonance Frequency (Hz))
    Empty = Tissues.Empty;
    WhiteMatter = Tissues.WhiteMatter;
    Fat = Tissues.Fat;      

    % Tissue Identification Labels
%     black = [0 0 0];        % Empty
    white = [255 255 255];  % Tissue 1
    red = [237 28 36];      % Tissue 2
    
    % Read Tissue Image
    t = imread(Filename);
    
    % Create Tissue Modls
    Rows = size(t, 1);
    Cols = size(t, 2);
    Tissue = Empty;
    for r = 1:Rows
        for c = 1:Cols
            if sum(squeeze(t(r,c,:)) == white') == 3
                Tissue(r,c) = Fat;
            elseif sum(squeeze(t(r,c,:)) == red') == 3
                Tissue(r,c) = WhiteMatter;
            else
                Tissue(r,c) = Empty;
            end
        end
    end           
end

function [Tissue] = GenerateResonancePSF(Tissue, FoSD)
    Rows = size(Tissue,1);
    Cols = size(Tissue,2);
    
    for r = 1:Rows
        for c = 1:Cols
            Tissue(r,c).Fo = Tissue(r,c).Fo + FoSD * randn(1, 1);
        end
    end
end

function [Tissue] = GenerateResonanceGradient(Tissue, MaxFreqz)
    Rows = size(Tissue,1);
    Cols = size(Tissue,2);
    
    phi = linspace(0,MaxFreqz,Rows);
    
    for r = 1:Rows
        for c = 1:Cols
            Tissue(r,c).Fo = Tissue(r,c).Fo + phi(r);
        end
    end
end

function [SequenceParameters] = GenerateResolutionAndGradientAmplitudes(SequenceParameters)
    
    P = SequenceParameters;
        
    % Calculate Resolution and Position
    P.x = linspace(0, P.FOVx, P.Nx);
    P.y = linspace(0, P.FOVy, P.Ny);
    P.dx = P.FOVx / P.Nx;
    P.dy = P.FOVy / P.Ny;
    
    % Calculate kx and ky limits
    % ky = (P.Gamma / (2*pi)) * P.Gyi * tau;
    % kx = (P.Gamma / (2*pi)) * P.GxMax * dt;
    P.dkx = 1 / P.FOVx;
    P.dky = 1 / P.FOVy;
    P.kxMax = 1 / (2 * P.dx);
    P.kyMax = 1 / (2 * P.dy);
    
    % Calculate Gradient Amplitudes
    % FOVx = 1 / ((Gamma/(2*pi))* Gx * dt)
    % FOVy = 1 / ((Gamma/(2*pi))* Gyi * tau)
    P.GxMax = 1/ ( (P.Gamma/(2*pi)) * P.FOVx * P.dt);
    P.Gyi = 1 / ( (P.Gamma/(2*pi)) * P.FOVy * P.tau);
    P.GyMax = P.Gyi * P.Ny / 2;
    
    SequenceParameters = P;
end

function [Mxy] = GenerateMxySignal(Tissue, SequenceParameters)

    P = SequenceParameters;
    TR = P.TR;              % TR
    TE = P.TE;              % TE
    Nr = P.Nr;              % Number of Repetitions
    alpha = P.Alpha;        % Tip Angle
    phi = P.Phi;            % Phi
    dphi = P.dPhi;          % dPhi
    
    % Alpha Tip and Phi Phase Angle
   	Rtip = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;      
            (1-cos(alpha))*cos(phi)*sin(phi)  cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi) ;
            sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];
    
    % Initialize
    M0 = Tissue.ProtonDensity;
    T1 = Tissue.T1;
    T2 = Tissue.T2;
    b = 2 * pi * Tissue.Fo * TR;
    M = [0 0 M0]';
      
    for n = 1:Nr   % Interate through tips
        % Alpha Degree Tip
        M = Rtip * M;  

        % T1, T2 Relaxation
        E = diag([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)]);
        Ez = [0;0;M0*(1-exp(-TR/T1))];
        % Off Resonance Precesion
        P = [cos(b),sin(b),0; -sin(b),cos(b),0; 0,0,1];
        % Single after Relaxation and Off-Resonace Precession
        M = P * E * M + Ez;

        % Increasing Phase each excitation for SSPF
        phi = phi + dphi;
        Rtip = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;      
                (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
                sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];
    end

    % After One Last Tip, and T1,T2 Relaxation and Beta Precession
 	M = Rtip * M;
  	E = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]); Ez = [0;0;M0*(1-exp(-TE/T1))];
 	P = [cos(b * TE/TR),sin(b * TE/TR),0; -sin(b * TE/TR),cos(b * TE/TR),0; 0,0,1];
 	M = P * E * M + Ez;

    % Sample Mxy
    Mxy = [0 0]';
    Mxy(1) = M(1);
    Mxy(2) = M(2);
end

function [Mc] = GenerateMxyImage(Tissue, SequenceParameters)
    % Initialize Signal Vector
    P = SequenceParameters;
    Rows = P.Ny;
    Cols = P.Nx;
    Mxy = zeros(Rows, Cols, 2);
    Mc = zeros(Rows, Cols);
   
    % Generate Mxy for each Sample Location without Gradients
    count = 0
    for r = 1:Rows
        for c = 1:Cols
            count = count + 1
            Mxy(r,c,:) = GenerateMxySignal(Tissue(r,c), P);
            Mc(r,c) = Mxy(r,c,1) + 1i * Mxy(r,c,2);
        end
    end
end

function [S] = GenerateSignalFromMxyImage(Tissue, SequenceParameters)
    
    % Initialize Signal Vector
    P = SequenceParameters;
    Rows = P.Ny;
    Cols = P.Nx;
    Mxy = zeros(Rows, Cols, 2);
    Mc = zeros(Rows, Cols);
    M = zeros(Rows, Cols);
    S = zeros(Rows, Cols);
   
    % Generate Mxy for each Sample Location without Gradients
    count = 0
    for r = 1:Rows
        for c = 1:Cols
            count = count + 1
            Mxy(r,c,:) = GenerateMxySignal(Tissue(r,c), P);
            Mc(r,c) = Mxy(r,c,1) + 1i * Mxy(r,c,2);
        end
    end
    
    % Plot Mxy Singal Image
    if(0)
        PlotMagnitudePhaseImage(Mc)
    end
    
    % Sample Singal 
    GyFactor = linspace(-1, 1, Rows);
    for r = 1:Rows
            % Phase Encode
            dky = (P.Gamma/(2*pi)) * P.Gyi * P.tau;
            ky(r) = (P.Gamma / (2*pi)) * GyFactor(r) * P.GyMax * P.tau; 
%             ky(r) = (r - Rows/2) * dky;
%             ky(r) = r * dky;
    end
    
    Tx = P.dt * P.Nx;
    dkx = (P.Gamma / (2*pi)) * P.GxMax * P.dt;
    kxShift = (P.Gamma / (2*pi)) * P.GxMax * Tx / 2;
    for c = 1:Cols
            % Frequency Encode
            t = c * P.dt;
            kx(c) = (P.Gamma / (2*pi)) * P.GxMax * t - kxShift;
    end
       
%     for r = -Rows/2:1:Rows/2-1
%         ky(r+Rows/2+1) = (r/Rows) * 2 * P.kyMax;
%     end
%     
%     for c = -Cols/2:1:Cols/2-1
%         kx(c+Cols/2+1) = (c/Cols) * 2 * P.kxMax;
%     end
    
    % Get a Sample for each kx and ky
    count2 = 0;
    for kr = 1:Rows
        for kc = 1:Cols
            count2 = count2 + 1
            
            % Get a single sample in k-space
            for r = 1:Rows
                for c = 1:Cols
                    M(r,c) = Mc(r,c) * exp(-2 * pi * 1i * (kx(kc) * P.x(c) + ky(kr) * P.y(r) ));
                end
            end

            % Sample Singal
            S(kr,kc) = sum(sum(M));
        end
    end
    

end

function [Spectrum] = PlotOffResonanceSpectrum(Tissues, SequenceParameters, MaxFreq)
    
    % Initalize
    P = SequenceParameters;
    FoMax = MaxFreq;                    % Set Max Freqz
    Ns = 201;                           % Set Number of Sample Points
    N = length(Tissues);                % Get Number of Tissues
    S = zeros(N,Ns);                    % Initalize Signals Vectors
    phi = linspace(-FoMax,FoMax,Ns);    % Create a Freqz Vector
    figure(); hold on;                  % Initalize Figure
    
    % Generate Spectrum for each Tissue
    for n = 1:N
        for ns = 1:Ns                                % For each Tissue
            tissue = Tissues(n);
            tissue.Fo = tissue.Fo + phi(ns);         % Set Off-Resonant Freqz
            Mxy = GenerateMxySignal(tissue, P);      % Get Mxy Signal
            S(n,ns) = Mxy(1) + 1i * Mxy(2);
        end
    end
    
    subplot(2,1,1);
    plot(phi, abs(S(1,:)), '-', phi, abs(S(2,:)), '--'); 
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Magnitude');
    legend(Tissues(1).Name, Tissues(2).Name);
    subplot(2,1,2);
    plot(phi, angle(S(1,:)), '-', phi, angle(S(2,:)), '--'); 
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Magnitude');
    hold off;
end

function [Spectrum] = PlotdPhiSpectrum(Tissues, SequenceParameters, MaxPhase)
    
    % Initalize
    P = SequenceParameters;
    Ns = 201;                                 % Set Number of Sample Points
    N = length(Tissues);                      % Get Number of Tissues
    S = zeros(N,Ns);                          % Initalize Signals Vectors
    phi = linspace(-MaxPhase,MaxPhase,Ns);    % Create a Phase Vector
    figure(); hold on;                        % Initalize Figure
    
    % Generate Spectrum for each Tissue
    for n = 1:N
        for ns = 1:Ns                                    % For each Tissue
            Pn = P;
            Pn.dPhi = Pn.dPhi + phi(ns);
            Mxy = GenerateMxySignal(Tissues(n), Pn);      % Get Mxy Signal
            S(n,ns) = Mxy(1) + 1i * Mxy(2);
        end
    end
    
    subplot(2,1,1);
    plot(phi, abs(S'));  
    xlabel('Off Spectrum Phase');
    ylabel('Signal Magnitude');
    legend(Tissues(1).Name, Tissues(2).Name);
    subplot(2,1,2);
    plot(phi, angle(S'));
    xlabel('Off Spectrum Phase');
    ylabel('Signal Angle');
    hold off;
end

function [Water, Fat] = TwoPointDixon(Image1, Image2)

    Water = 0.5*(Image1+Image2);
    Fat = 0.5*(Image1-Image2);
end

function [Water, Fat] = ThreePointDixon(Image1, Image2, Image3)

    phi = angle(conj(Image1).*Image3)/2;
    Water = 0.5*(Image1+Image2.*exp(-1i*phi));
    Fat = 0.5*(Image1-Image2.*exp(-1i*phi));
end

function [Water, Fat] = TwoPointDixonWithSort(Image1, Image2, Image3, PhaseState1, PhaseState2, PhaseState3)

    Water = 0.5*(Image1+Image2);
    Fat = 0.5*(Image1-Image2);
end

function [Water, Fat] = ThreePointDixonWithSort(Image1, Image2, Image3, PhaseState1, PhaseState2, PhaseState3)

    [Ny,Nx] = size(Image1);
    Water = zeros(Ny, Nx);
    Fat = zeros(Ny, Nx);
    for r = 1:Ny
        for c = 1:Nx
            
            PhaseCount = PhaseState1(r,c) + PhaseState2(r,c) + PhaseState3(r,c);
            if PhaseCount == 2  % There are 2 In-Phase Images
                
                if PhaseState1(r,c)
                    Pixel1 = Image1(r,c);
                    if PhaseState2(r,c)
                       Pixel2 = Image3(r,c);
                       Pixel3 = Image2(r,c);
                    else
                       Pixel2 = Image2(r,c);
                       Pixel3 = Image3(r,c);
                    end
                    
                else
                    Pixel1 = Image2(r,c);
                    Pixel2 = Image1(r,c);
                    Pixel3 = Image3(r,c);
                end
                
            else              % There are 2 Out-Phase Images
                
                if PhaseState1(r,c)
                    Pixel1 = Image2(r,c);
                    Pixel2 = Image1(r,c);
                    Pixel3 = Image3(r,c);
                else
                    
                    Pixel1 = Image1(r,c);
                    if PhaseState2(r,c)
                       Pixel2 = Image2(r,c);
                       Pixel3 = Image3(r,c);
                    else
                       Pixel2 = Image3(r,c);
                       Pixel3 = Image2(r,c);
                    end
                    
                end

            end
            
            phi = angle(conj(Pixel1)*Pixel3)/2;
            Water(r,c) = 0.5*(Pixel1+Pixel2*exp(-1i*phi));
            Fat(r,c) = 0.5*(Pixel1-Pixel2*exp(-1i*phi));
        end
    end
    

end


function [S] = AddComplexSignalNoise(Signal, Mean, StdDev)
    [Ny, Nx] = size(Signal);
    Noise = Mean + StdDev.*randn(Ny, Nx) + 1i * (Mean + StdDev.*randn(Ny, Nx));
    S = Signal + Noise;
end

function PlotMagnitudePhaseImage(S)
	figure();
  	subplot(1,2,1);
 	imshow(abs(S),[]);         % Signal Magnitude
 	subplot(1,2,2);
 	imshow(angle(S),[]);       % Signal Phase
end

function PlotMagnitiudePhaseFor2Images(S, S2)
    figure();
    subplot(2,2,1);
 	imshow(abs(S),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(2,2,2);
 	imshow(angle(S),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(2,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
    title('Magitude Out-Phase Image');
 	subplot(2,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
end

function PlotMagnitiudePhaseFor3Images(S, S2, S3)
    figure();
    subplot(3,2,1);
 	imshow(abs(S),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(3,2,2);
 	imshow(angle(S),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(3,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
    title('Magitude Out-Phase Image');
 	subplot(3,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(3,2,5);
 	imshow(abs(S3),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(3,2,6);
 	imshow(angle(S3),[]);       % Signal Phase
    title('Phase In-Phase Image');
end

function PlotMagnitiudePhaseFor4Images(S, S2, S3, S4)
    figure();
    subplot(4,2,1);
 	imshow(abs(S),[]);          % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(4,2,2);
 	imshow(angle(S),[]);        % Signal Phase
    title('Phase In-Phase Image');
    subplot(4,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
    title('Magitude Out-Phase Image');
 	subplot(4,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(4,2,5);
 	imshow(abs(S3),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(4,2,6);
 	imshow(angle(S3),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(4,2,7);
    imshow(abs(S4),[]);         % Signal Magnitude
    title('Magitude Out-Phase Image');
 	subplot(4,2,8);
 	imshow(angle(S4),[]);       % Signal Phase
    title('Magitude Out-Phase Image');
end

function PlotMagnitiudePhaseFor5Images(S, S2, S3, S4, S5)
    figure();
    subplot(5,2,1);
 	imshow(abs(S),[]);          % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(5,2,2);
 	imshow(angle(S),[]);        % Signal Phase
    title('Phase In-Phase Image');
    subplot(5,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
    title('Magitude Out-Phase Image');
 	subplot(5,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(5,2,5);
 	imshow(abs(S3),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(5,2,6);
 	imshow(angle(S3),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(5,2,7);
    imshow(abs(S4),[]);         % Signal Magnitude
    title('Magitude Out-Phase Image');
 	subplot(5,2,8);
 	imshow(angle(S4),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(5,2,9);
 	imshow(abs(S5),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(5,2,10);
 	imshow(angle(S5),[]);       % Signal Phase
    title('Phase In-Phase Image');
end

function PlotMagnitiudePhaseFor6Images(S, S2, S3, S4, S5, S6)
    figure();
    subplot(6,2,1);
 	imshow(abs(S),[]);          % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(6,2,2);
 	imshow(angle(S),[]);        % Signal Phase
    title('Phase In-Phase Image');
    subplot(6,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
    title('Magitude Out-Phase Image');
 	subplot(6,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(6,2,5);
 	imshow(abs(S3),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(6,2,6);
 	imshow(angle(S3),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(6,2,7);
    imshow(abs(S4),[]);         % Signal Magnitude
    title('Magitude Out-Phase Image');
 	subplot(6,2,8);
 	imshow(angle(S4),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(6,2,9);
 	imshow(abs(S5),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(6,2,10);
 	imshow(angle(S5),[]);       % Signal Phase
    title('Phase In-Phase Image');
    subplot(6,2,11);
 	imshow(abs(S6),[]);         % Signal Magnitude
    title('Magitude In-Phase Image');
 	subplot(6,2,12);
 	imshow(angle(S6),[]);       % Signal Phase
    title('Phase In-Phase Image');
end

function [SNR] = CalculateSNR(Signal, TissueModel, Tissue)
    Rows = size(TissueModel,1);
    Cols = size(TissueModel,2);
    
    s = [];
    n = [];
    for r = 1:Rows
        for c = 1:Cols
            if strcmp(Tissue.Name, TissueModel(r,c).Name)
                s = [s Signal(r,c)];
            elseif strcmp('Empty', TissueModel(r,c).Name)
                n = [n Signal(r,c)];
            end
        end
    end
    
    mu = mean(s);
    sigmaN = std(n);
    
    SNR = abs(mu)/sigmaN;
end


function [CNR] = CalculateCNR(Signal, TissueModel, Tissue1, Tissue2)
    Rows = size(TissueModel,1);
    Cols = size(TissueModel,2);
    
    s1 = [];
    s2 = [];
    n = [];
    for r = 1:Rows
        for c = 1:Cols
            if strcmp(Tissue1.Name, TissueModel(r,c).Name)
                s1 = [s1 Signal(r,c)];
            elseif strcmp(Tissue2.Name, TissueModel(r,c).Name)
                s2 = [s2 Signal(r,c)];
            elseif strcmp('Empty', TissueModel(r,c).Name)
                n = [n Signal(r,c)];
            end
        end
    end
    
    mean1 = mean(s1);
    mean2 = mean(s2);
    sigmaN = std(n);
    
    CNR = abs(mean1-mean2)/sigmaN;
end
