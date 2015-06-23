% Noise Analysis
function NoiseAnalysis()

    % Generate Tissue Model
    tissues = TissueElements();

    % Inital Parameters
    alpha = pi/4; phi = 0; Nr = 200; TR = 10/1000; TE = 5/1000; 
    dphis = [0 pi/4 pi/2 3*pi/2]; 
    dt = tissues.Fat.F0;
    TEs = [TE - 1/(2*dt), TE, TE + 1/(2*dt), TE + 1/dt];

    FoSD             = 0;       % Off-Resonance Std. Dev. in Hz
    MaxGradFreq      = 0;     % Max Frequency for Off-Resonance Linear Gradient
    NoiseSD          = 0.005;       % Signal Noise Std. Dev.

    N = 64;
    data = zeros(N,N,8);

    % Generate a Complex Noise Map
    noise = linspace(-3*NoiseSD,3*NoiseSD,N);
    cnoise = zeros(N,N);
    for r = 1:N
        for c = 1:N
            cnoise(r,c) = noise(c) + 1i * noise(r);
        end
    end

    % Generate a Gradient based on f0Noise in x-dir
    f0Noise = linspace(-3*FoSD,3*FoSD,N);
    f0Noise2D = zeros(N,N);
    for r = 1:N
        for c = 1:N
            f0Noise2D(r,c) = f0Noise(c);
        end
    end

    % Generate Off-Resonace Frequency Gradient in y-dir
    f0Gradient = zeros(N,N);
    gradient = linspace(0, MaxGradFreq, N);
    for r = 1:N
        for c = 1:N
            f0Gradient(r,c) = gradient(r);
        end
    end

    for r = 1:N
        for c = 1:N
            M0 = 1; T1 = tissues.WhiteMatter.T1; T2 = tissues.WhiteMatter.T2; beta = 2*pi*TR*(tissues.WhiteMatter.F0 + f0Noise2D(r,c) + f0Gradient(r,c));
            data(r,c,1) = SSFP(beta, M0, alpha, phi, dphis(1), TR, TEs(1), T1, T2, Nr) + cnoise(r,c);
            data(r,c,2) = SSFP(beta, M0, alpha, phi, dphis(2), TR, TEs(2), T1, T2, Nr) + cnoise(r,c);
            data(r,c,3) = SSFP(beta, M0, alpha, phi, dphis(3), TR, TEs(3), T1, T2, Nr) + cnoise(r,c);
            data(r,c,4) = SSFP(beta, M0, alpha, phi, dphis(4), TR, TEs(4), T1, T2, Nr) + cnoise(r,c);
            M0 = 1; T1 = tissues.Fat.T1; T2 = tissues.Fat.T2; beta = 2*pi*TR*(tissues.Fat.F0 + f0Noise2D(r,c) + f0Gradient(r,c));
            data(r,c,5) = SSFP(beta, M0, alpha, phi, dphis(1), TR, TEs(1), T1, T2, Nr) + cnoise(r,c);
            data(r,c,6) = SSFP(beta, M0, alpha, phi, dphis(2), TR, TEs(2), T1, T2, Nr) + cnoise(r,c);
            data(r,c,7) = SSFP(beta, M0, alpha, phi, dphis(3), TR, TEs(3), T1, T2, Nr) + cnoise(r,c);
            data(r,c,8) = SSFP(beta, M0, alpha, phi, dphis(4), TR, TEs(4), T1, T2, Nr) + cnoise(r,c);
        end
    end

    img = [data(:,:,1) data(:,:,5)];
    img2 = [data(:,:,2) data(:,:,6)];
    img3 = [data(:,:,3) data(:,:,7)];
    img4 = [data(:,:,4) data(:,:,8)];

    PlotMRImage(cnoise);
    PlotMRImage(img, img2, img3, img4);
    [ water, fat ] = SSFPSeparation(img, img2, img3, img4, [1,1,64,64]);
    PlotMRImage(water, fat);

end