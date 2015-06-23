function PlotTissueSpectrumVsdphi_dphiAdjusted(alpha, TE, TR) 
    
    dt = 1/(2*428);
    M0 = 1; phi = 0;
    Nr = 200; Ns = 100; BetaMax = 2*pi;
    tissues = TissueElements();
    
    T1 = [tissues.Fat.T1 tissues.WhiteMatter.T1];
    T2 = [tissues.Fat.T2 tissues.WhiteMatter.T2];
    f0 = [tissues.Fat.F0 tissues.WhiteMatter.F0];
    dphi = [0 1/4 2/4 3/4 4/4] * 2 *pi;
    
    for n=1:1
        Mf = SSFPSpectrum( M0, alpha, phi, dphi(1), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf2 = SSFPSpectrum( M0, alpha, phi, dphi(2), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf3 = SSFPSpectrum( M0, alpha, phi, dphi(3), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf4 = SSFPSpectrum( M0, alpha, phi, dphi(4), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        [Mf5, Beta] = SSFPSpectrum( M0, alpha, phi, dphi(5), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));

        f = Beta / TR / (2*pi);

        figure();
        subplot(2,2,1);
        plot(f, abs(Mf), '-', f, abs(Mf2), '--', f, abs(Mf3), '--', f, abs(Mf4), '--', f, abs(Mf5), '--'); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)');
        ylabel('Signal Magnitude');
        legend('dphi=0', 'dphi=pi/2', 'dphi=pi', 'dphi=3 pi/2','dphi=2 pi');
        subplot(2,2,3);
        plot(f, angle(Mf), '-', f, angle(Mf2), '--', f, angle(Mf3), '--', f, angle(Mf4), '--', f, angle(Mf5), '--'); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)');
        ylabel('Signal Phase');
        hold off;
    end
    
    dt = 1/(2*428);
    M0 = 1; phi = 0;
    Nr = 200; Ns = 100; BetaMax = 2*pi;
    tissues = TissueElements();
    
    T1 = [tissues.Fat.T1 tissues.WhiteMatter.T1];
    T2 = [tissues.Fat.T2 tissues.WhiteMatter.T2];
    f0 = [tissues.Fat.F0 tissues.WhiteMatter.F0];
    dphi = [0 1/4 2/4 3/4 4/4] * 2 *pi;
    phis = -[0 1/4 2/4 3/4 4/4] * pi;
    
    for n=1:1
        Mf = SSFPSpectrum( M0, alpha, phis(1), dphi(1), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf2 = SSFPSpectrum( M0, alpha, phis(2), dphi(2), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf3 = SSFPSpectrum( M0, alpha, phis(3), dphi(3), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf4 = SSFPSpectrum( M0, alpha, phis(4), dphi(4), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        [Mf5, Beta] = SSFPSpectrum( M0, alpha, phis(5), dphi(5), TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));

        f = Beta / TR / (2*pi);

        subplot(2,2,2);
        plot(f, abs(Mf), '-', f, abs(Mf2), '--', f, abs(Mf3), '--', f, abs(Mf4), '--', f, abs(Mf5), '--'); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)');
        ylabel('Signal Magnitude');
        legend('dphi=0', 'dphi=pi/2', 'dphi=pi', 'dphi=3 pi/2','dphi=2 pi');
        subplot(2,2,4);
        plot(f, angle(Mf), '-', f, angle(Mf2), '--', f, angle(Mf3), '--', f, angle(Mf4), '--', f, angle(Mf5), '--'); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)');
        ylabel('Signal Phase');
        hold off;
    end
end