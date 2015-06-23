function PlotSpectraComparsion(alpha, phis, dphis, TEs, TR)
    % PlotSpectra(TE, TR, dphi) generates a SSFP Spectrum plot for Fat and 
    % WhiteMatter
    
    tissues = TissueElements();
    TE = TR/2;
    
    phis = -[0 1/4 2/4 3/4] * pi;
    dphis = [0 1/4 2/4 3/4] * 2 * pi; 
    df = abs(tissues.Fat.F0);
    dt = 1 / (2 * df);
%     TEs = [TE - 3*dt/2, TE - dt/2, TE + dt/2, TE + 3*dt/2];
    TEs = [TE, TE + dt, TE, TE + dt];
    
    M0 = 1;
    Nr = 200; Ns = 100; BetaMax = 2*pi;
    
    T1 = tissues.Fat.T1; T2 = tissues.Fat.T2; f0 = tissues.Fat.F0;
    [Mf,Beta] = SSFPSpectrum( M0, alpha, phis(1), dphis(1), TR, TEs(1), T1, T2, Nr, Ns, BetaMax, f0);
    Mf2 = SSFPSpectrum( M0, alpha, phis(2), dphis(2), TR, TEs(2), T1, T2, Nr, Ns, BetaMax, f0);
    Mf3 = SSFPSpectrum( M0, alpha, phis(3), dphis(3), TR, TEs(3), T1, T2, Nr, Ns, BetaMax, f0);
    Mf4 = SSFPSpectrum( M0, alpha, phis(4), dphis(4), TR, TEs(4), T1, T2, Nr, Ns, BetaMax, f0);
    Mf5 = SSFPSpectrum( M0, alpha, phis(1), dphis(1), TR, TEs(2), T1, T2, Nr, Ns, BetaMax, f0);
    Mf6 = SSFPSpectrum( M0, alpha, phis(3), dphis(3), TR, TEs(4), T1, T2, Nr, Ns, BetaMax, f0);
    
    T1 = tissues.WhiteMatter.T1; T2 = tissues.WhiteMatter.T2; f0 = tissues.WhiteMatter.F0;
    Mw = SSFPSpectrum( M0, alpha, phis(1), dphis(1), TR, TEs(1), T1, T2, Nr, Ns, BetaMax, f0);
    Mw2 = SSFPSpectrum( M0, alpha, phis(2), dphis(2), TR, TEs(2), T1, T2, Nr, Ns, BetaMax, f0);
    Mw3 = SSFPSpectrum( M0, alpha, phis(3), dphis(3), TR, TEs(3), T1, T2, Nr, Ns, BetaMax, f0);
    Mw4 = SSFPSpectrum( M0, alpha, phis(4), dphis(4), TR, TEs(4), T1, T2, Nr, Ns, BetaMax, f0);
    Mw5 = SSFPSpectrum( M0, alpha, phis(1), dphis(1), TR, TEs(2), T1, T2, Nr, Ns, BetaMax, f0);
    Mw6 = SSFPSpectrum( M0, alpha, phis(3), dphis(3), TR, TEs(4), T1, T2, Nr, Ns, BetaMax, f0);
    
    f = Beta / TR / (2*pi);
    
    figure();
    subplot(2,1,1);
    plot(f, abs(Mf), '-', f, abs(Mf2), '--', f, abs(Mf3), '--', f, abs(Mf4), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Magnitude','FontSize',16);
    subplot(2,1,2);
    plot(f, angle(Mf), '-', f, angle(Mf2), '--', f, angle(Mf3), '--', f, angle(Mf4), '--');
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Phase','FontSize',16);
    hold off;
    
    figure();
    subplot(2,1,1);
    plot(f, abs(Mw), '-', f, abs(Mw2), '--', f, abs(Mw3), '--',  f, abs(Mw4), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Magnitude','FontSize',16);
    subplot(2,1,2);
    plot(f, angle(Mw), '-', f, angle(Mw2), '--', f, angle(Mw3), '--', f, angle(Mw4), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Phase','FontSize',16);
    hold off;
    
    [waterMw, fatMw] = SSFPSeparation(Mw, Mw2, Mw3, Mw4, 1);
    [waterMf, fatMf] = SSFPSeparation(Mf, Mf2, Mf3, Mf4, 1);
    
    figure();
    subplot(2,1,1);
    plot(f, abs(waterMw), '-', f, abs(fatMw), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Magnitude','FontSize',16);
    subplot(2,1,2);
    plot(f, angle(waterMw), '-', f, angle(fatMw), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Phase','FontSize',16);
    hold off;
    
    figure();
    subplot(2,1,1);
    plot(f, abs(waterMf), '-', f, abs(fatMf), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Magnitude','FontSize',16);
    subplot(2,1,2);
    plot(f, angle(waterMw), '-', f, angle(fatMw), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Phase','FontSize',16);
    hold off;
    
    [waterMw, fatMw] = SSFPSeparation(Mw, Mw2, Mw3, Mw4, 2);
    [waterMf, fatMf] = SSFPSeparation(Mf, Mf2, Mf3, Mf4, 2);
    
    figure();
    subplot(2,1,1);
    plot(f, abs(waterMw), '-', f, abs(fatMw), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Magnitude','FontSize',16);
    subplot(2,1,2);
    plot(f, angle(waterMw), '-', f, angle(fatMw), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Phase','FontSize',16);
    hold off;
    
    figure();
    subplot(2,1,1);
    plot(f, abs(waterMf), '-', f, abs(fatMf), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Magnitude','FontSize',16);
    subplot(2,1,2);
    plot(f, angle(waterMw), '-', f, angle(fatMw), '--'); 
    set(gca,'FontSize',16);
    xlabel('Off Spectrum Frequency (Hz)','FontSize',16);
    ylabel('Signal Phase','FontSize',16);
    hold off;
    
end    

