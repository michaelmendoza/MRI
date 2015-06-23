function PlotTissueSpectrumTE(alpha, dphi, TR, TE, TE_Range, f0_Offset)
    % PlotSpectra(TE, TR, dphi) generates a SSFP Spectrum plot for Fat and 
    % WhiteMatter
    
    M0 = 1; phi = 0;
    Nr = 200; Ns = 100;
    tissues = TissueElements();
    
    T1 = tissues.Fat.T1; T2 = tissues.Fat.T2; f0 = tissues.Fat.F0 + f0_Offset;
    [ Mf, TE ] = SSFPSpectrumTE( M0, alpha, phi, dphi, TR, TE, T1, T2, Nr, Ns, f0, TE_Range);

    
    T1 = tissues.WhiteMatter.T1; T2 = tissues.WhiteMatter.T2; f0 = tissues.WhiteMatter.F0 + f0_Offset;
    [ Mw, TE ] = SSFPSpectrumTE( M0, alpha, phi, dphi, TR, TE, T1, T2, Nr, Ns, f0, TE_Range); 
    
    figure();
    plot(TE*1000, angle(Mf), '-', TE*1000, angle(Mw), '--'); 
    set(gca,'FontSize',14);
    xlabel('TE (ms)','FontSize',14);
    ylabel('Signal Phase','FontSize',14);
    legend(tissues.Fat.Name, tissues.WhiteMatter.Name,'FontSize',14);
    
end      