function PlotTissueSpectrumPhi(alpha, TR, TE, f0_Offset)
    % PlotSpectra(TE, TR, dphi) generates a SSFP Spectrum plot for Fat and 
    % WhiteMatter
    
    M0 = 1; dphi = 0;
    Nr = 200; Ns = 100;
    tissues = TissueElements();
    
    T1 = tissues.Fat.T1; T2 = tissues.Fat.T2; f0 = tissues.Fat.F0 + f0_Offset;
    [ Mf, dphi ] = SSFPSpectrumPhi( M0, alpha, dphi, TR, TE, T1, T2, Nr, Ns, f0);

    
    figure();
    plot(dphi, angle(Mf), '-'); 
    set(gca,'FontSize',14);
    xlabel('Phi (Radians)','FontSize',14);
    ylabel('Signal Phase','FontSize',14);
    hold off;    
end      