function PlotTissueSpectrum(alpha, phi, dphi, TE, TR)
    % PlotSpectra(TE, TR, dphi) generates a SSFP Spectrum plot for Fat and 
    % WhiteMatter
    
    M0 = 1;
    Nr = 200; Ns = 100; BetaMax = 2*pi;
    tissues = TissueElements();
    
    T1 = tissues.Fat.T1; T2 = tissues.Fat.T2; f0 = tissues.Fat.F0;
    Mf = SSFPSpectrum( M0, alpha, phi, dphi, TR, TE, T1, T2, Nr, Ns, BetaMax, f0);


    T1 = tissues.WhiteMatter.T1; T2 = tissues.WhiteMatter.T2; f0 = tissues.WhiteMatter.F0;
    [Mw, Beta] = SSFPSpectrum( M0, alpha, phi, dphi, TR, TE, T1, T2, Nr, Ns, BetaMax, f0);

    f = Beta / TR / (2*pi);
    
    figure();
    subplot(2,1,1);
    plot(f, abs(Mf), '-', f, abs(Mw), '--'); 
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Magnitude');
    legend(tissues.Fat.Name, tissues.WhiteMatter.Name);
    subplot(2,1,2);
    plot(f, angle(Mf), '-', f, angle(Mw), '--'); 
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Phase');
    hold off;
    
end