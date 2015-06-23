function PlotTissueSpectrumVsTE(alpha, dphi, TE, TR) 
    
    dt = 1/(2*428);
    M0 = 1; phi = 0;
    Nr = 200; Ns = 100; BetaMax = 2*pi;
    tissues = TissueElements();
    
    T1 = [tissues.Fat.T1 tissues.WhiteMatter.T1];
    T2 = [tissues.Fat.T2 tissues.WhiteMatter.T2];
    f0 = [tissues.Fat.F0 tissues.WhiteMatter.F0];
    
    figure();
    for n=1:2
        Mf = SSFPSpectrum( M0, alpha, phi, dphi, TR, TE, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf2 = SSFPSpectrum( M0, alpha, phi, dphi, TR, TE+dt/2, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf3 = SSFPSpectrum( M0, alpha, phi, dphi, TR, TE+dt, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        Mf4 = SSFPSpectrum( M0, alpha, phi, dphi, TR, TE+3*dt/2, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));
        [Mf5, Beta] = SSFPSpectrum( M0, alpha, phi, dphi, TR, TE+2*dt, T1(n), T2(n), Nr, Ns, BetaMax, f0(n));

        f = Beta / TR / (2*pi);
 
        subplot(2,2,2*n-1);
        plot(f, abs(Mf), '-', f, abs(Mf2), '--', f, abs(Mf3), '--', f, abs(Mf4), '--', f, abs(Mf5), '--'); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)','FontSize', 14);
        ylabel('Signal Magnitude','FontSize', 14);
        legend('TE=5ms', 'TE=5ms+dt/2', 'TE=5ms+dt', 'TE=5ms+3*dt/2','TE=5ms+2*dt','FontSize', 14);
        subplot(2,2,2*n);
        plot(f, angle(Mf), '-', f, angle(Mf2), '--', f, angle(Mf3), '--', f, angle(Mf4), '--', f, angle(Mf5), '--'); 
        set(gca,'FontSize',14);
        xlabel('Off Spectrum Frequency (Hz)','FontSize', 14);
        ylabel('Signal Phase','FontSize', 14);
        hold off;
    end
end