function PlotSpectrumForAbstract(alpha, phi, dphi, TE, TR)
    % PlotSpectra(TE, TR, dphi) generates a SSFP Spectrum plot
    
    M0 = 1;
    Nr = 200; Ns = 100; BetaMax = 2*pi;
    tissues = TissueElements();
    
%     T1 = tissues.Fat.T1; T2 = tissues.Fat.T2; f0 = tissues.Fat.F0;
    T1 = tissues.WhiteMatter.T1; T2 = tissues.WhiteMatter.T2; f0 = tissues.WhiteMatter.F0;
    [Mf, Beta] = SSFPSpectrum( M0, alpha, phi, dphi, TR, TE, T1, T2, Nr, Ns, BetaMax, f0);

    f = Beta / TR / (2*pi);
    
    figure();
    subplot(2,1,1);
    plot(f, abs(Mf), '-k'); 
   
    hold on; 
    scatter(-50, 0.01, 100, 'red', 'x','LineWidth', 3);
    scatter(-25, 0.12, 100, 'black', 'x','LineWidth', 3);
    scatter(0, 0.1495, 100, 'black', 'x','LineWidth', 3);
    scatter(25, 0.1208, 100, 'black', 'x','LineWidth', 3);
    hold off;
    
    set(gca,'FontSize',10);
    ylabel('Signal');
    subplot(2,1,2);
    plot(f, angle(Mf), '-k'); 
    
    hold on; 
    scatter(-50, 0.01, 100, 'red', 'x','LineWidth', 3);
    scatter(-25, 1.62, 100, 'black', 'x','LineWidth', 3);
    scatter(0, 1.56, 100, 'black', 'x','LineWidth', 3);
    scatter(25, 1.51, 100, 'black', 'x','LineWidth', 3);
    hold off;
    
    
    
    set(gca,'FontSize',10);
    xlabel('Off Spectrum(Hz)');
    ylabel('Phase');
 
   
    
end