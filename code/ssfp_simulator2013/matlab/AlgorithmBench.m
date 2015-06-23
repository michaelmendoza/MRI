% Algorithm Bench
clear;

% Inital Parameters
alpha = pi/3; phi = 0; dphi = pi; Nr = 200;
TR = 10/1000; TE = TR/2; f0 = 0;

FoSD             = 0.05;       	% Off-Resonance Std. Dev. in Hz
MaxGradFreq      = 500; %50; 	% Max Frequency for Off-Resonance Linear Gradient
NoiseSD          = 0.004; %0.004;       % Signal Noise Std. Dev.

% Generate Tissue Model
filenames = {'TissueMasks/RectPhantom64Tissues.png', 'TissueMasks/Phantom64Tissues.png', 'TissueMasks/Phantom256Tissues.png'};
filename = filenames{2};
tissues = TissueElements();
tissue = TissueModel(filename, tissues, FoSD, MaxGradFreq);

%% ------------------------------------------------------------------------
if 0
   PlotSpectraComparsion(alpha, phi, dphi, TE, TR); 
end

%% ------------------------------------------------------------------------
if 0
   PlotSpectrumForAbstract(alpha, phi, dphi, TE, TR); 
end

%% ------------------------------------------------------------------------
% Plot Tissue Spectum for Fat and WhiteMatter
if 0
    PlotTissueSpectrum(alpha, phi, dphi, TE, TR);
end

%% ------------------------------------------------------------------------
% Plot Spectrum for fat at different TEs
if 0
    PlotTissueSpectrumVsTE(alpha, dphi, TE, TR);
end

%% ------------------------------------------------------------------------
% Plot TE Spectrum
if 0
    TE_Range = 1 / (2 * 428); f0 = 100;
    PlotTissueSpectrumTE(alpha, dphi, TR, TE, TE_Range, f0);
end

%% ------------------------------------------------------------------------
% Plot Spectrum for different Phi
if 0
    PlotTissueSpectrumVsPhi(alpha, dphi, TE, TR);
end

%% ------------------------------------------------------------------------
% Plot phi Spectrum
if 0
    f0 = 50;
    PlotTissueSpectrumPhi(alpha, TR, TE, f0);
end

%% ------------------------------------------------------------------------
% Plot Spectrum for different dphi
if 0
   PlotTissueSpectrumVs_dphi(alpha, TE, TR)  
end

%% ------------------------------------------------------------------------
% Plot dphi Spectrum
if 0
   f0 = 100;
   PlotTissueSpectrum_dPhi(alpha, TR, TE, f0); 
end

%% ------------------------------------------------------------------------
% Plot Spectrum for different dphi and phi 
if 0
    PlotTissueSpectrumVsPhi_dPhi(alpha, TE, TR); 
end

%% ------------------------------------------------------------------------
% Plot Spectrum for different dphi and dphi adjusted
if 0
    PlotTissueSpectrumVsdphi_dphiAdjusted(alpha, TE, TR); 
end

%% ------------------------------------------------------------------------
% Generate Simple Single Image
if 0
    img = SSFPImage(tissue, alpha, phi, dphi, TR, TE, Nr, NoiseSD);
    PlotMRImage(img);
end

%% ------------------------------------------------------------------------
% Plot Spectrum of 4 Phase Cycled Images
if 0

    phis = -[0 1/4 2/4 3/4] * pi;
    dphis = [0 1/4 2/4 3/4] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
    TEs = [TE - 3*dt/2, TE - dt/2, TE + dt/2, TE + 3*dt/2];

    PlotTissueSpectrum(alpha, phis(1), dphis(1), TEs(1), TR);
    PlotTissueSpectrum(alpha, phis(2), dphis(2), TEs(2), TR);
    PlotTissueSpectrum(alpha, phis(3), dphis(3), TEs(3), TR);
    PlotTissueSpectrum(alpha, phis(4), dphis(4), TEs(4), TR);
end

%% ------------------------------------------------------------------------
% Plot Spectrum of 5 Phase Cycled Images
if 0
    dphis = [0 1/5 2/5 3/5 4/5] * 2* pi; 
    dt = abs(tissues.Fat.F0);
    TEs = [ TE - 1/dt, TE - 1/(2*dt), TE, TE + 1/(2*dt), TE + 1/dt];

    PlotTissueSpectrum(alpha, phis(1), dphis(1), TEs(1), TR);
    PlotTissueSpectrum(alpha, phis(2), dphis(2), TEs(2), TR);
    PlotTissueSpectrum(alpha, phis(3), dphis(3), TEs(3), TR);
    PlotTissueSpectrum(alpha, phis(4), dphis(4), TEs(4), TR);
    PlotTissueSpectrum(alpha, phis(5), dphis(5), TEs(5), TR);
end

%% ------------------------------------------------------------------------
% Generate 4 Phase Cycled Images
if 0
    phis = -[0 1/4 2/4 3/4] * pi;
    dphis = [0 1/4 2/4 3/4] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
    TEs = [TE - 3*dt/2, TE - dt/2, TE + dt/2, TE + 3*dt/2];
    
    img = SSFPImage(tissue, alpha, phis(1), dphis(1), TR, TEs(1), Nr, NoiseSD, f0);
    img2 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(2), Nr, NoiseSD, f0);
    img3 = SSFPImage(tissue, alpha, phis(3), dphis(3), TR, TEs(3), Nr, NoiseSD, f0);
    img4 = SSFPImage(tissue, alpha, phis(4), dphis(4), TR, TEs(4), Nr, NoiseSD, f0);
    PlotMRImage(img, img2, img3, img4);
end

%% ------------------------------------------------------------------------
% Generate 6 Phase Cycled Images
if 1
    phis = -[0 1/4 2/4 3/4] * pi;
    dphis = [0 1/4 2/4 3/4] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
%     TEs = [TE - 3*dt/2, TE - dt/2, TE + dt/2, TE + 3*dt/2, TE - dt/2 - 2*dt, TE + 3*dt/2 - 2*dt];
    TEs = [TE, TE + dt, TE, TE + dt, TE - dt, TE - dt];

    img = SSFPImage(tissue, alpha, phis(1), dphis(1), TR, TEs(1), Nr, NoiseSD, f0);
    img2 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(2), Nr, NoiseSD, f0);
    img3 = SSFPImage(tissue, alpha, phis(3), dphis(3), TR, TEs(3), Nr, NoiseSD, f0);
    img4 = SSFPImage(tissue, alpha, phis(4), dphis(4), TR, TEs(4), Nr, NoiseSD, f0);
    img5 = [];
    img6 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(5), Nr, NoiseSD, f0);
    img7 = [];
    img8 = SSFPImage(tissue, alpha, phis(4), dphis(4), TR, TEs(6), Nr, NoiseSD, f0);
    PlotMRImage(img, img2, img3, img4, img6, img8);
end

%% ------------------------------------------------------------------------
% Generate 5 Phase Cycled Images
if 0   
    phis = -[0 1/5 2/5 3/5 4/5] * pi;
    dphis = [0 1/5 2/5 3/5 4/5] * 2* pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
    TEs = [TE - 2*dt, TE - dt, TE, TE + dt, TE + 2 * dt];
    
    img = SSFPImage(tissue, alpha, phis(1), dphis(1), TR, TEs(1), Nr, NoiseSD, f0);
    img2 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(2), Nr, NoiseSD, f0);
    img3 = SSFPImage(tissue, alpha, phis(3), dphis(3), TR, TEs(3), Nr, NoiseSD, f0);
    img4 = SSFPImage(tissue, alpha, phis(4), dphis(4), TR, TEs(4), Nr, NoiseSD, f0);
    img5 = SSFPImage(tissue, alpha, phis(5), dphis(5), TR, TEs(5), Nr, NoiseSD, f0);
    PlotMRImage(img, img2, img3, img4, img5);
end

%% ------------------------------------------------------------------------
% Plot SSFP Cycled Images
if 0
   PlotMRImage(img, img2, img3, img4); 
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation
if 0
    [ water, fat ] = DixonSeparation(img, img2);
    PlotMRImage(water, fat);
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation
if 0
    [ water, fat ] = DixonSeparation(img, img2, img3);
    PlotMRImage(water, fat);
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation with 4 Images all combinations
if 0
    [ water, fat ] = DixonSeparation(img, img2);
    [ water2, fat2 ] = DixonSeparation(img2, img3);
    [ water3, fat3 ] = DixonSeparation(img3, img4);
    [ water4, fat4 ] = DixonSeparation(img4, img);
    
    PlotMRImage(water, fat);
    PlotMRImage(water2, fat2);
    PlotMRImage(water3, fat3);
    PlotMRImage(water4, fat4);    
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation with 5 Images all combinations
if 0    
    [ water, fat ] = DixonSeparation(img, img2);
    [ water2, fat2 ] = DixonSeparation(img2, img3);
    [ water3, fat3 ] = DixonSeparation(img3, img4);
    [ water4, fat4 ] = DixonSeparation(img4, img5);
    
    PlotMRImage(water, fat);
    PlotMRImage(water2, fat2);
    PlotMRImage(water3, fat3);
    PlotMRImage(water4, fat4);
end

%% ------------------------------------------------------------------------
% Separate water and fat images with SSFP Separation
if 1
    [water, fat ] = SSFPSeparation(img, img2, img3, img4, 1);
    PlotMRImage(water, fat);
    
    figure();
    subplot(1,6,1); imshow(abs(img),[]);
    subplot(1,6,2); imshow(abs(img2),[]);
    subplot(1,6,3); imshow(abs(img3),[]);
    subplot(1,6,4); imshow(abs(img4),[]);
    subplot(1,6,5); imshow(abs(water),[]);
    subplot(1,6,6); imshow(abs(fat),[]);
    
    
%     [waterIwaterEst, waterIfatE, waterIsigmaN] = CalculateTissueSignalAverage(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     [fatIwaterE, fatIfatE, fatIsigmaN] = CalculateTissueSignalAverage(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [waterIwaterEst, waterIfatE, waterIsigmaN fatIwaterE, fatIfatE, fatIsigmaN]
%     
%     [SNRW, sigmaN] = CalculateSNR(water, tissue, tissues.WhiteMatter);
%     SNRF = CalculateSNR(fat, tissue, tissues.Fat);
%     CNRW = CalculateCNR(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRF = CalculateCNR(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [SNRW SNRF CNRW CNRF sigmaN]
end

%% ------------------------------------------------------------------------
% Separate water and fat images with SSFP Separation
if 1
    [water, fat ] = SSFPSeparationWithFieldMap(img, img2, img3, img4, img5, img6, img7, img8, 1);
    PlotMRImage(water, fat);
    
    figure();
    subplot(1,6,1); imshow(abs(img),[]);
    subplot(1,6,2); imshow(abs(img2),[]);
    subplot(1,6,3); imshow(abs(img3),[]);
    subplot(1,6,4); imshow(abs(img4),[]);
    subplot(1,6,5); imshow(abs(water),[]);
    subplot(1,6,6); imshow(abs(fat),[]);
    
    figure();
    subplot(1,2,1); imshow(angle(conj(img6).*img2) / 2, []);
    subplot(1,2,2); imshow(angle(conj(img8).*img4) / 2, []);
    
%     [waterIwaterEst, waterIfatE, waterIsigmaN] = CalculateTissueSignalAverage(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     [fatIwaterE, fatIfatE, fatIsigmaN] = CalculateTissueSignalAverage(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [waterIwaterEst, waterIfatE, waterIsigmaN fatIwaterE, fatIfatE, fatIsigmaN]
%     
%     [SNRW, sigmaN] = CalculateSNR(water, tissue, tissues.WhiteMatter);
%     SNRF = CalculateSNR(fat, tissue, tissues.Fat);
%     CNRW = CalculateCNR(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRF = CalculateCNR(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [SNRW SNRF CNRW CNRF sigmaN]
end

%% ------------------------------------------------------------------------
% Separate water and fat images with SSFP Separation with 5 Images
if 0
	[ water, fat ] = SSFPSeparationWith5Images(img, img2, img3, img4, img5, 0);
    PlotMRImage(water, fat);
    
    [waterIwaterEst, waterIfatE, waterIsigmaN] = CalculateTissueSignalAverage(water, tissue, tissues.WhiteMatter, tissues.Fat);
    [fatIwaterE, fatIfatE, fatIsigmaN] = CalculateTissueSignalAverage(fat, tissue, tissues.WhiteMatter, tissues.Fat);
    [waterIwaterEst, waterIfatE, waterIsigmaN fatIwaterE, fatIfatE, fatIsigmaN]
    
%     [SNRW, sigmaN] = CalculateSNR(water, tissue, tissues.WhiteMatter);
%     SNRF = CalculateSNR(fat, tissue, tissues.Fat);
%     CNRW = CalculateCNR(water, tissue, tissues.WhiteMatter, tissues.Fat);
%     CNRF = CalculateCNR(fat, tissue, tissues.WhiteMatter, tissues.Fat);
%     [SNRW SNRF CNRW CNRF sigmaN]
end

%% ------------------------------------------------------------------------
% Noise Analysis
if 0
    NoiseAnalysis();
end

