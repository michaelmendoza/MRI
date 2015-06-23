function PlotTissueSpectrumVsPhi_dPhi(alpha, TE, TR) 
    
    dt = 1/(2*428);
    M0 = 1;
    Nr = 200; Ns = 100; BetaMax = 2*pi;
    tissues = TissueElements();
    
    T1 = [tissues.Fat.T1 tissues.WhiteMatter.T1];
    T2 = [tissues.Fat.T2 tissues.WhiteMatter.T2];
    f0 = [tissues.Fat.F0 tissues.WhiteMatter.F0];
    phi = -[0 1/4 2/4 3/4 4/4] * pi;
    dphi = [0 1/4 2/4 3/4 4/4] * 2 * pi;
    
    for n=1:1
        Mf = SSFPSpectrum( M0, alpha, phi(1), dphi(1), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf2 = SSFPSpectrum( M0, alpha, phi(2), dphi(2), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf3 = SSFPSpectrum( M0, alpha, phi(3), dphi(3), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf4 = SSFPSpectrum( M0, alpha, phi(4), dphi(4), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        [Mf5, Beta] = SSFPSpectrum( M0, alpha, phi(5), dphi(5), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));

        f = Beta / TR / (2*pi);

        figure();
        subplot(2,1,1);
        plot(f, abs(Mf), '-', f, abs(Mf2), '--', f, abs(Mf3), ':', f, abs(Mf4), '-.', f, abs(Mf5)); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)');
        ylabel('Signal Magnitude');
        legend('phi=0', 'phi=pi/2', 'phi=pi', 'phi=3 pi/2','phi=2 pi');
        subplot(2,1,2);
        plot(f, angle(Mf), '-', f, angle(Mf2), '--', f, angle(Mf3), ':', f, angle(Mf4), '-.', f, angle(Mf5)); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)');
        ylabel('Signal Phase');
        hold off;
    end
    
    dt = 1/(2*428);
    M0 = 1;
    Nr = 200; Ns = 100; BetaMax = 2*pi;
    tissues = TissueElements();
    
    T1 = [tissues.Fat.T1 tissues.WhiteMatter.T1];
    T2 = [tissues.Fat.T2 tissues.WhiteMatter.T2];
    f0 = [tissues.Fat.F0 tissues.WhiteMatter.F0];
    phi = -[0 0 0 0 0] * pi;
    dphi = [0 1/4 2/4 3/4 4/4] * 2 * pi;
    phaseAdjust = [0 1/4 2/4 3/4 4/4] * pi;
    
    for n=1:1
        Mf = SSFPSpectrum( M0, alpha, phi(1), dphi(1), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf2 = SSFPSpectrum( M0, alpha, phi(2), dphi(2), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf3 = SSFPSpectrum( M0, alpha, phi(3), dphi(3), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf4 = SSFPSpectrum( M0, alpha, phi(4), dphi(4), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        [Mf5, Beta] = SSFPSpectrum( M0, alpha, phi(5), dphi(5), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));

        K = length(Mf);
        for k = 1:K
            Mf(k) = Mf(k) * exp(-1i*phaseAdjust(1));
            Mf2(k) = Mf2(k) * exp(-1i*phaseAdjust(2));
            Mf3(k) = Mf3(k) * exp(-1i*phaseAdjust(3));
            Mf4(k) = Mf4(k) * exp(-1i*phaseAdjust(4));
            Mf5(k) = Mf5(k) * exp(-1i*phaseAdjust(5));
        end
            
        f = Beta / TR / (2*pi);

        figure();
        subplot(2,1,1);
        plot(f, abs(Mf), '-', f, abs(Mf2), '--', f, abs(Mf3), '--', f, abs(Mf4), '--', f, abs(Mf5), '--'); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)');
        ylabel('Signal Magnitude');
        legend('phi=0', 'phi=pi/2', 'phi=pi', 'phi=3 pi/2','phi=2 pi');
        subplot(2,1,2);
        plot(f, angle(Mf), '-', f, angle(Mf2), '--', f, angle(Mf3), '--', f, angle(Mf4), '--', f, angle(Mf5), '--'); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)');
        ylabel('Signal Phase');
        hold off;
    end
end