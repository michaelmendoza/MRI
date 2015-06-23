function PlotTissueSpectrum_dPhi(alpha, TR, TE, f0_Offset)
    % PlotSpectra(TE, TR, dphi) generates a SSFP Spectrum plot for Fat and 
    % WhiteMatter
    
    M0 = 1; phi = 0;
    Nr = 200; Ns = 100;
    tissues = TissueElements();
    
    T1 = tissues.Fat.T1; T2 = tissues.Fat.T2; f0 = tissues.Fat.F0 + f0_Offset;
    [ Mf, dphi ] = SSFPSpectrum_dPhi( M0, alpha, phi, TR, TE, T1, T2, Nr, Ns, f0);

    T1 = tissues.WhiteMatter.T1; T2 = tissues.WhiteMatter.T2; f0 = tissues.WhiteMatter.F0 + f0_Offset;
    [ Mw, dphi ] = SSFPSpectrum_dPhi( M0, alpha, phi, TR, TE, T1, T2, Nr, Ns, f0);
    
    figure();
    subplot(2,2,1);
    plot(dphi, abs(Mf), '-', dphi, abs(Mw), '--'); 
    set(gca,'FontSize',14);
    xlabel('dphi (Radians)');
    ylabel('Signal Magnitude');
    xlim([0 2*pi]);
    legend(tissues.Fat.Name, tissues.WhiteMatter.Name);
    subplot(2,2,3);
    plot(dphi, angle(Mf), '-', dphi, angle(Mw), '--'); 
    set(gca,'FontSize',14);
    xlabel('dphi (Radians)');
    ylabel('Signal Phase');
    xlim([0 2*pi]);
    
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
    end
    
    M0 = 1; phi = 0;
    Nr = 200; Ns = 100;
    tissues = TissueElements();
    
    T1 = tissues.Fat.T1; T2 = tissues.Fat.T2; f0 = tissues.Fat.F0 + f0_Offset;
    [ Mf1, dphi ] = SSFPSpectrum_dPhi( M0, pi/8, phi, TR, TE, T1, T2, Nr, Ns, f0);
    [ Mf2, dphi ] = SSFPSpectrum_dPhi( M0, pi/6, phi, TR, TE, T1, T2, Nr, Ns, f0);
    [ Mf3, dphi ] = SSFPSpectrum_dPhi( M0, pi/4, phi, TR, TE, T1, T2, Nr, Ns, f0);
    [ Mf4, dphi ] = SSFPSpectrum_dPhi( M0, pi/3, phi, TR, TE, T1, T2, Nr, Ns, f0);
    [ Mf5, dphi ] = SSFPSpectrum_dPhi( M0, pi/2, phi, TR, TE, T1, T2, Nr, Ns, f0);
    
    figure();
    subplot(2,1,1);
    plot(dphi, abs(Mf1), '-', dphi, abs(Mf2), '--', dphi, abs(Mf3), '--', dphi, abs(Mf4), '--', dphi, abs(Mf5), '--');
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Magnitude');
    legend('alpha =pi/8','alpha =pi/6','alpha =pi/4','alpha =pi/3','alpha =pi/2');
    subplot(2,1,2);
    plot(dphi, angle(Mf1), '-', dphi, angle(Mf2), '--', dphi, angle(Mf3), '--', dphi, angle(Mf4), '--', dphi, angle(Mf5), '--'); 
    set(gca,'FontSize',14);
    xlabel('dphi (Radians)');
    ylabel('Signal Phase');
    hold off;
    
end      