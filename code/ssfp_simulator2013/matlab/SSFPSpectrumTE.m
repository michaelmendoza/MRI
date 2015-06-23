function [ Mc, TE ] = SSFPSpectrumTE( M0, alpha, phi, dphi, TR, TE, T1, T2, Nr, Ns, f0, TE_Range)
    % SSFPSpectrum() generates an SSFP transverse magnetization spectrum 
    % with default parameters. 
    %
    % SSFPSpectrum( M0, alpha, phi, dphi, TR, TE, T1, T2, Nr, Ns, f0, TE_Max) 
    % generates an SSFP transverse magnetization spectrum with the following 
    % parameters:
    %   M0 = initial magnetization strength
    %   alpha = RF pulse tip in degrees
    %   phi = RF pulse phase in degrees
    %   dphi = RF pulse phase incrementation per TR (for phase cycling)
    %   TR = repetition time
    %   TE = echo time
    %   T1 = longitudal relaxation time constant
    %   T2 = transverse relaxation time constant
    %   Nr = number of excitations
    %   Ns = number of samples for spectrum
    %   f0 = off-resonance precession frequency
    %   TE_Range = TE range, the TE used is TE - TE_Range to TE + TE_Range

    if not(exist('M0') && exist('alpha') && exist('phi') && exist('dphi') && exist('TR') && exist('TE') && exist('T1') && exist('T2') && exist('Nr'))
        M0 = 1; alpha = pi/2; phi = 0; dphi = 0;
        TR = 10/1000; TE = TR/2;
        T1 = 790/1000; T2 = 92/1000; 
        Nr = 200;
    end

    if not(exist('Ns'))
       Ns = 200;
    end
    
    if not(exist('f0'))
        f0 = 0;
    end

    Mc = [];
    count = 1;
    for TEn = linspace(-TE_Range, TE_Range, Ns) + TE         % Interate through beta values
        beta0 = 2*pi*TR*f0;                                 % Calculate addition phase due to f0
        Mc(count) = SSFP(beta0, M0, alpha, phi, dphi, TR, TEn, T1, T2, Nr);
        count = count + 1;
    end
    
    TE = linspace(-TE_Range, TE_Range, Ns) + TE;
end
