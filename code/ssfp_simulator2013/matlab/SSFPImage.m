function [img] = SSFPImage(Tissue, alpha, phi, dphi, TR, TE, Nr, NoiseSD, f0)
    % SSFPImage generates an SSPF magnetization image based on supplied
    % tissue and SSFP sequenece paramters.
    %
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
    %   BetaMax = max off-resonance precession phase between excitations
    %   f0 = off-resonance precession frequency
    %   NoiseSD = Std. Dev. of complex noise

    if (~exist('alpha') && ~exist('phi') && ~exist('dphi') && ~exist('TR') && ~exist('T2') && ~exist('Nr'))
        alpha = pi/2; phi = 0; dphi = 0;
        TR = 10/1000; TE = TR/2;
        Nr = 200;
    end
    
    if (~exist('NoiseSD', 'var'))
        NoiseSD = 0.0;
    end
    
    if (~exist('f0', 'var'))
        f0 = 0.0;
    end
    
    Rows = size(Tissue, 1);
    Cols = size(Tissue, 2);
    img = zeros(Rows, Cols);
    for r = 1:Rows
        for c = 1:Cols
            beta0 = 2*pi*TR*f0;
            M0 = Tissue(r,c).ProtonDensity; T1 = Tissue(r,c).T1; T2 = Tissue(r,c).T2; beta = 2*pi*TR*Tissue(r,c).F0;
            img(r,c) = SSFP(beta + beta0, M0, alpha, phi, dphi, TR, TE, T1, T2, Nr);
        end
    end
    
    img = AddComplexSignalNoise(img, 0, NoiseSD);
end

function [S] = AddComplexSignalNoise(Signal, Mean, StdDev)
    [Ny, Nx] = size(Signal);
    Noise = Mean + StdDev.*randn(Ny, Nx) + 1i * (Mean + StdDev.*randn(Ny, Nx));
    S = Signal + Noise;
end
