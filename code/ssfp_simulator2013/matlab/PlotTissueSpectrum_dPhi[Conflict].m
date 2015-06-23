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
    subplot(2,1,1);
    plot(dphi, abs(Mf), '-', dphi, abs(Mw), '--'); 
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Magnitude');
    legend(tissues.Fat.Name, tissues.WhiteMatter.Name);
    subplot(2,1,2);
    plot(dphi, angle(Mf), '-', dphi, angle(Mw), '--'); 
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Phase');
    hold off;
    
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
    plot(dphi, abs(Mf1), '-', dphi, abs(Mf2), '--', dphi, abs(Mf3), ':', dphi, abs(Mf4), '-.', dphi, abs(Mf5));
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Magnitude');
    legend('alpha =pi/8','alpha =pi/6','alpha =pi/4','alpha =pi/3','alpha =pi/2');
    subplot(2,1,2);
    plot(dphi, angle(Mf1), '-', dphi, angle(Mf2), '--', dphi, angle(Mf3), ':', dphi, angle(Mf4), '-.', dphi, angle(Mf5)); 
    set(gca,'FontSize',14);
    xlabel('Off Spectrum Frequency (Hz)');
    ylabel('Signal Phase');
    hold off;
    
end      