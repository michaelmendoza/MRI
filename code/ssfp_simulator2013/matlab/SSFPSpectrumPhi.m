function [ Mc, dphi ] = SSFPSpectrumPhi( M0, alpha, dphi, TR, TE, T1, T2, Nr, Ns, f0)
    % SSFPSpectrum() generates an SSFP transverse magnetization spectrum 
    % with default parameters. 
    %
    % SSFPSpectrum( M0, alpha, phi, dphi, TR, TE, T1, T2, Nr, Ns, f0, TE_Max) 
    % generates an SSFP transverse magnetization spectrum with the following 
    % parameters:
    %   M0 = initial magnetization strength
    %   alpha = RF pulse tip in degrees
    %   phi = RF pulse phase in degrees
    %   TR = repetition time
    %   TE = echo time
    %   T1 = longitudal relaxation time constant
    %   T2 = transverse relaxation time constant
    %   Nr = number of excitations
    %   Ns = number of samples for spectrum
    %   f0 = off-resonance precession frequency


    Mc = [];
    count = 1;
    for phi = linspace(0, 2*pi, Ns);                       % Interate through beta values
        beta0 = 2*pi*TR*f0;                                 % Calculate addition phase due to f0
        Mc(count) = SSFP(beta0, M0, alpha, phi, dphi, TR, TE, T1, T2, Nr);
        count = count + 1;
    end
    
    dphi = linspace(0, 2*pi, Ns);
end
