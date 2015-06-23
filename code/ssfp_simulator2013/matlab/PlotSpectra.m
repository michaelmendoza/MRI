function PlotSpectra( M, TR, BetaMax )
    % PlotSpectra(M, TR, BetaMax) generates a SSFP Spectrum plot
    %   M       - Signal magnization vector at various beta values
    %   TR      - Repetition time
    %   BetaMax - Max beta used in magnetization vector
    
    figure();
    Beta = linspace(0, BetaMax, length(M));
    f = Beta / TR / (2*pi);
    subplot(2,1,1); plot(f, abs(M));
    ylabel('Magnitude');
    xlabel('Frequency (Hertz)');
    subplot(2,1,2); plot(f, angle(M));
    xlabel('Frequency (Hertz)');
    ylabel('Angle');

end
