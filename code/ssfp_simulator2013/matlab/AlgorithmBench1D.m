% Algorithm Bench
clear; 

% Inital Parameters
alpha = pi/3; phi = 0; dphi = pi; Nr = 200;
TR = 10/1000; TE = TR/2; f0 = 0;

FoSD             = 0.5;           % Off-Resonance Std. Dev. in Hz
MaxGradFreq      = 500;           % Max Frequency for Off-Resonance Linear Gradient
NoiseSD          = 0.004;         % Signal Noise Std. Dev.

% Generate Tissue Model
filenames = {'TissueMasks/RectPhantom64Tissues.png', 'TissueMasks/Phantom64Tissues.png', 'TissueMasks/Phantom256Tissues.png'};
filename = filenames{2};
tissues = TissueElements();
tissue = TissueModel1D(tissues, FoSD, MaxGradFreq);

%% ------------------------------------------------------------------------
% Plot Tissue Spectum for Fat and WhiteMatter
if 0
    PlotTissueSpectrum(alpha, dphi, TE, TR);
end

%% ------------------------------------------------------------------------
% Plot Spectrum of 4 Phase Cycled Images
if 0
    phis = -[0 1/4 2/4 3/4] * pi;
    dphis = [0 1/4 2/4 3/4] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
    TEs = [TE - 3*dt/2, TE - dt/2, TE + dt/2, TE + 3*dt/2];
    TEs = [TE, TE + dt, TE, TE + dt];
    
    PlotTissueSpectrum(alpha, phis(1), dphis(1), TEs(1), TR);
    PlotTissueSpectrum(alpha, phis(2), dphis(2), TEs(2), TR);
    PlotTissueSpectrum(alpha, phis(3), dphis(3), TEs(3), TR);
    PlotTissueSpectrum(alpha, phis(4), dphis(4), TEs(4), TR);
end

%% ------------------------------------------------------------------------
% Generate Simple Single Image
if 0
    img = SSFPImage(tissue, alpha, phi, dphi, TR, TE, Nr, NoiseSD);
    figure(1);
    subplot(2,1,1); plot(abs(img)); title('Magnitude');
    subplot(2,1,2); plot(angle(img)); title('Phase');
end

%% ------------------------------------------------------------------------
% Generate 4 Phase Cycled Images
if 0
    phis = -[0 1/4 2/4 3/4] * pi;
    dphis = [0 1/4 2/4 3/4] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
%     TEs = [TE - 3*dt/2, TE - dt/2, TE + dt/2, TE + 3*dt/2];
    TEs = [TE, TE + dt, TE, TE + dt];


    img = SSFPImage(tissue, alpha, phis(1), dphis(1), TR, TEs(1), Nr, NoiseSD, f0);
    img2 = SSFPImage(tissue, alpha, phis(2), dphis(2), TR, TEs(2), Nr, NoiseSD, f0);
    img3 = SSFPImage(tissue, alpha, phis(3), dphis(3), TR, TEs(3), Nr, NoiseSD, f0);
    img4 = SSFPImage(tissue, alpha, phis(4), dphis(4), TR, TEs(4), Nr, NoiseSD, f0);
    
    figure();
    subplot(4,1,1); plot(abs(img)); title('Magnitude 1');
    subplot(4,1,2); plot(abs(img2)); title('Magnitude 2');
    subplot(4,1,3); plot(abs(img3)); title('Magnitude 3');
    subplot(4,1,4); plot(abs(img4)); title('Magnitude 4');
    figure();
    subplot(4,1,1); plot(angle(img)); title('Phase 1');
    subplot(4,1,2); plot(angle(img2)); title('Phase 2');
    subplot(4,1,3); plot(angle(img3)); title('Phase 3');
    subplot(4,1,4); plot(angle(img4)); title('Phase 4');
    figure();
    subplot(2,1,1); hold on; plot(abs(img),'b'); plot(abs(img2),'g'); plot(abs(img3),'r'); plot(abs(img4),'c'); hold off;
    subplot(2,1,2); hold on; plot(angle(img),'b'); plot(angle(img2),'g'); plot(angle(img3),'r'); plot(angle(img4),'c'); hold off;
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
    
    figure();
    subplot(4,1,1); plot(abs(img)); title('Magnitude 1');
    subplot(4,1,2); plot(abs(img2)); title('Magnitude 2');
    subplot(4,1,3); plot(abs(img3)); title('Magnitude 3');
    subplot(4,1,4); plot(abs(img4)); title('Magnitude 4');
    figure();
    subplot(4,1,1); plot(angle(img)); title('Phase 1');
    subplot(4,1,2); plot(angle(img2)); title('Phase 2');
    subplot(4,1,3); plot(angle(img3)); title('Phase 3');
    subplot(4,1,4); plot(angle(img4)); title('Phase 4');
    figure();
    subplot(2,1,1); hold on; plot(abs(img),'b'); plot(abs(img2),'g'); plot(abs(img3),'r'); plot(abs(img4),'c'); hold off;
    subplot(2,1,2); hold on; plot(angle(img),'b'); plot(angle(img2),'g'); plot(angle(img3),'r'); plot(angle(img4),'c'); hold off;
end

%% ------------------------------------------------------------------------
% Generate 4 Phase Cycled Images Summed and Subtracted
if 0
    N = length(img);
    line = pi * ones(N,1);
    figure();
    subplot(2,1,1); hold on; plot(line,'r--'); plot(0*line,'r--'); plot(-line,'r--'); plot(angle(img)-angle(img2)); title('img1 - img2'); hold off;
    subplot(2,1,2); hold on; plot(line,'r--'); plot(0*line,'r--'); plot(-line,'r--'); plot(angle(img)+angle(img2)); title('img1 + img2'); hold off;
    figure();
    subplot(2,1,1); hold on; plot(line,'r--'); plot(0*line,'r--'); plot(-line,'r--'); plot(angle(img2)-angle(img3)); title('img2 - img3'); hold off;
    subplot(2,1,2); hold on; plot(line,'r--'); plot(0*line,'r--'); plot(-line,'r--'); plot(angle(img2)+angle(img3)); title('img2 + img3'); hold off;
    figure();
    subplot(2,1,1); hold on; plot(line,'r--'); plot(0*line,'r--'); plot(-line,'r--'); plot(angle(img3)-angle(img4)); title('img3 - img4'); hold off;
    subplot(2,1,2); hold on; plot(line,'r--'); plot(0*line,'r--'); plot(-line,'r--'); plot(angle(img3)+angle(img4)); title('img3 + img4'); hold off;
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation
if 0
    [ water, fat ] = DixonSeparation(img, img2);
    
    figure();
    subplot(2,1,1); plot(abs(water)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat)); title('Fat - Magnitude');
end

%% ------------------------------------------------------------------------
% Separate water and fat images with Dixon Separation
if 0
    [ water, fat ] = DixonSeparation(img, img2, img3);
    
    figure();
    subplot(2,1,1); plot(abs(water)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat)); title('Fat - Magnitude');
end

%% ------------------------------------------------------------------------
% Separate water and fat images with 2-Point Dixon Separation all combinations
if 0
    [ water, fat ] = DixonSeparation(img, img2);
    [ water2, fat2 ] = DixonSeparation(img2, img3);
    [ water3, fat3 ] = DixonSeparation(img3, img4);
    [ water4, fat4 ] = DixonSeparation(img4, img);
    
    figure();
    subplot(2,1,1); plot(abs(water)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat)); title('Fat - Magnitude');
    figure();
    subplot(2,1,1); plot(abs(water2)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat2)); title('Fat - Magnitude');
    figure();
    subplot(2,1,1); plot(abs(water3)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat3)); title('Fat - Magnitude');
    figure();
    subplot(2,1,1); plot(abs(water4)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat4)); title('Fat - Magnitude');
end

%% ------------------------------------------------------------------------
% Separate water and fat images with 3-Point Dixon Separation all combinations
if 0
    [img, img2, img3, img4] = AdjustSSFPConstantPhase1D(img, img2, img3, img4, [60 80]);

    [ water, fat ] = DixonSeparationWithFieldMap(img, img2, img, img3);
    [ water2, fat2 ] = DixonSeparationWithFieldMap(img2, img3, img2, img4);
    [ water3, fat3 ] = DixonSeparationWithFieldMap(img3, img4, img, img3);
    [ water4, fat4 ] = DixonSeparationWithFieldMap(img4, img, img2, img4);
    
    figure();
    subplot(2,1,1); plot(abs(water)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat)); title('Fat - Magnitude');
    figure();
    subplot(2,1,1); plot(abs(water2)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat2)); title('Fat - Magnitude');
    figure();
    subplot(2,1,1); plot(abs(water3)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat3)); title('Fat - Magnitude');
    figure();
    subplot(2,1,1); plot(abs(water4)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat4)); title('Fat - Magnitude');
end

%% ------------------------------------------------------------------------
% Separate water and fat images with SSFP Separation
if 1

    [ water2, fat2 ] = SSFPSeparation(img, img2, img3, img4, 1);
    
    figure();
    subplot(2,1,1); plot(abs(water2)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat2)); title('Fat - Magnitude');
end

%% ------------------------------------------------------------------------
% Separate water and fat images with SSFP Separation
if 1

    [ water2, fat2 ] = SSFPSeparationWithFieldMap(img, img2, img3, img4, img5, img6, img7, img8, 1);
    
    figure();
    subplot(2,1,1); plot(abs(water2)); title('Water - Magnitude');
    subplot(2,1,2); plot(abs(fat2)); title('Fat - Magnitude');
end

