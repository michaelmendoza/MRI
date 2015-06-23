function SSPFBench()

    % ------------ Tissue (ms) ------------
    % Gray Matter  : T1 - 920, T2 - 100
    % White Matter : T1 - 790, T2 - 92
    % Fat:           T1 - 270, T2 - 85
    % Muscle:        T1 - 870, T2 - 47
    % ------------ Tissue (ms) ------------
    % Fat Water Spin Separation = 228 Hz
        
    alpha = pi/5;
    T1 = 790/1000; T2 = 92/1000; TR = 10/1000; 
    T1F = 270/1000; T2F = 85/1000;
    T1M = 870/1000; T2M = 47/1000;
    
%     SSPF(alpha, dphi, T1, T2, TR, 'b');
% 	 FEMR(alpha, -pi/2, pi/2, T1, T2, TR, 'b');

    M1 = SSPF(alpha, 0, T1, T2, TR, 'b');
    M2 = SSPF(alpha, pi/2, T1, T2, TR, 'r');
    M3 = SSPF(alpha, pi, T1, T2, TR,'g');
    M4 = SSPF(alpha, 3*pi/2, T1, T2, TR,'k');
    
    M1F = SSPF(alpha, 0, T1F, T2F, TR, 'b');   M1F180 = SSPF(alpha, pi, T1F, T2F, TR, 'b');
    M2F = SSPF(alpha, pi/2, T1F, T2F, TR, 'r');
    M3F = SSPF(alpha, pi, T1F, T2F, TR,'g');
    M4F = SSPF(alpha, 3*pi/2, T1F, T2F, TR,'k');

    M1M = SSPF(alpha, 0, T1M, T2M, TR, 'b');
    M2M = SSPF(alpha, pi/2, T1M, T2M, TR, 'r');
    M3M = SSPF(alpha, pi, T1M, T2M, TR,'g');
    M4M = SSPF(alpha, 3*pi/2, T1M, T2M, TR,'k');
    
    % Plot two Spectrum 
    clf; 
    hFig = figure(1);
    set(hFig, 'Position', [500 500 300 300])
    b = linspace(-2*pi, 2*pi, length(M1F));
    f =  b / (2 * pi * TR);
    subplot(2,1,1);  plot(f, abs(M1F180), 'b'); 
    set(gca,'FontSize',14)
    ylabel('Magnitude', 'FontSize', 16);

    subplot(2,1,2);  plot(f, angle(M1F180), 'b');  
    set(gca,'FontSize',14)
    xlabel('Off-Resonance(Hz)', 'FontSize', 14);   ylabel('Phase', 'FontSize', 16);
    
    hFig = figure(2);
    set(hFig, 'Position', [500 500 400 400])
    b = linspace(-pi, pi, length(M1));
    f =  b / (2 * pi * TR);
    subplot(2,2,1); plot(f, abs(M1), 'r'); hold on;
    subplot(2,2,2); plot(f, abs(M2), 'r'); hold on;
    subplot(2,2,3); plot(f, abs(M3), 'r'); hold on;
    subplot(2,2,4); plot(f, abs(M4), 'r'); hold on;
    
    subplot(2,2,1); plot(f, abs(M1F), 'b--'); 
    subplot(2,2,2); plot(f, abs(M2F), 'b--'); 
    subplot(2,2,3); plot(f, abs(M3F), 'b--'); 
    subplot(2,2,4); plot(f, abs(M4F), 'b--'); 
    
    subplot(2,2,1); plot(f, abs(M1M), 'k-.'); hold off;
    set(gca,'FontSize',14)
    ylabel('Magnitude', 'FontSize', 16);
    subplot(2,2,2); plot(f, abs(M2M), 'k-.'); hold off;
    set(gca,'FontSize',14)
    subplot(2,2,3); plot(f, abs(M3M), 'k-.'); hold off;   
    set(gca,'FontSize',14)
    xlabel('Off-Resonance (Hz)', 'FontSize', 14);
    ylabel('Magnitude', 'FontSize', 16);
    subplot(2,2,4); plot(f, abs(M4M), 'k-.'); hold off;
    set(gca,'FontSize',14)
    xlabel('Off-Resonance (Hz)', 'FontSize', 14);
    legend('White Matter', 'Fat', 'Muscle');

    
    M = [M1,M2,M3,M4];
    MF = [M1F,M2F,M3F,M4F];
    MM = [M1M,M2M,M3M,M4M];
    Mc = ComplexSum(M);
    MFc = ComplexSum(MF);
    MMc = ComplexSum(MM);
    
    hFig = figure(3);
    set(hFig, 'Position', [500 500 400 400])    
    f =  b / (2 * pi * TR);
    plot(f, abs(Mc),'r'); hold on;
    plot(f, abs(MFc),'b--'); 
    plot(f, abs(MMc),'k-.'); hold off;
    ylim([0.1 .8]);
    set(gca,'FontSize',14)
    xlabel('Off-Resonance (Hz)', 'FontSize', 14);
    ylabel('Magnitude', 'FontSize', 16);    
    legend('White Matter', 'Fat', 'Muscle');

    
%     M = [M1,M2,M3,M4];
%     figure(2); clf; hold on;
%     plot(MaxIntensity(M),'r');  
%     plot(ComplexSum(M),'g');    
%     plot(MagnitudeSum(M),'b');  
%     plot(SumOfSquares(M),'k');  
%     hold off;
%     legend('Max Intensity', 'Complex Sum', 'Magnitude Sum', 'Sum of Squares');
%     ylim([0 1]);
        
%     X = [ 1 2 0; 2 0 3; 0 3 4] 
%     RemoveNulls(X)
%     M = [M1,M2,M3,M4];
%     M = RemoveNulls(M);
%     figure(3); clf; hold on;
%     plot(MaxIntensity(M),'r');  
%     plot(ComplexSum(M),'g');    
%     plot(MagnitudeSum(M),'b');  
%     plot(SumOfSquares(M),'k');  
%     hold off;
%     legend('Max Intensity', 'Complex Sum', 'Magnitude Sum', 'Sum of Squares');
%     ylim([0 1]);
end

function [M] = SSPF(alpha, dphi, T1, T2, TR, plotstring)

    TE = TR/2;
    
    Nr = 200;  % Number of Repetitions
    Ns = 200; % Number of Samples for Off-Resonance Graph
    BetaMax = 2*pi; % Off Resonance - Free-Precession Angle (rad)
    
    M0 = 1; phi = 0; 
    
    Rx = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;
          (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
          sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];
      
    SampleCount = 0;
    for b = linspace(-BetaMax, BetaMax, Ns)   % Interate through Beta values

        M = [0,0,M0]'; % Starting M Vector
        SampleCount = SampleCount + 1
        for n = 1:Nr   % Interate through tips
            % Alpha Degree Tip
            Mtip = Rx*M;  

            % T1, T2 Relaxation
            E = diag([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)]);
            Ez = [0;0;M0*(1-exp(-TR/T1))];
            % Off Resonance Precesion
            P = [cos(b),sin(b),0; -sin(b),cos(b),0; 0,0,1];
            % Single after Relaxation and Off-Resonace Precession
            M = P * E * Mtip + Ez;

            phi = phi + dphi;
            Rx = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;      
                  (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
                  sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];
        end

        % After One Last Tip, and T1,T2 Relaxation and Beta Precession
        Mtip = Rx*M;
        E = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]); Ez = [0;0;M0*(1-exp(-TE/T1))];
        P = [cos(b * TE/TR),sin(b * TE/TR),0; -sin(b * TE/TR),cos(b * TE/TR),0; 0,0,1];
        M = P * E * Mtip + Ez;

        % Save Sample Point after Steady State is Reached
        Ms(1,SampleCount) = norm(M(1:2));
        Ms(2,SampleCount) = angle(M(1)+j*M(2));

    end

    M = Ms;
    
    Beta = linspace(-BetaMax, BetaMax, Ns);
    f = Beta / TR / (2*pi);
    PlotSpectra(f, M(1,:), M(2,:), plotstring);
    
    % Make a Column Vector
    M = MagAngleToComplex(Ms(1,:),Ms(2,:))';

end

function [Meven, Modd] = FEMR(alpha, dphi1, dphi2, T1, T2, TR, plotstring)

    TE = TR/2;
    
    Nr = 200;  % Number of Repetitions
    Ns = 200; % Number of Samples for Off-Resonance Graph
    BetaMax = 2*pi; % Off Resonance - Free-Precession Angle (rad)
    
    M0 = 1; phi = 0; 
    
    Rx = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;
          (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
          sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];
      
    SampleCount = 0;
    for b = linspace(-BetaMax, BetaMax, Ns)   % Interate through Beta values
        
        % Reset M Vector and Total Phi
        M = [0,0,M0]'; % Starting M Vector
        phi = 0;
        
        SampleCount = SampleCount + 1
        for n = 1:Nr   % Interate through tips
            % Alpha Degree Tip
            Mtip = Rx*M;  

            % T1, T2 Relaxation
            E = diag([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)]);
            Ez = [0;0;M0*(1-exp(-TR/T1))];
            % Off Resonance Precesion
            P = [cos(b),sin(b),0; -sin(b),cos(b),0; 0,0,1];
            % Single after Relaxation and Off-Resonace Precession
            M = P * E * Mtip + Ez;

            if ( mod(n, 2) == 0)
                phi = phi + dphi1;
            else
                phi = phi + dphi2;
            end
            Rx = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;      
                  (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
                  sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];
        
            % Acquire Even and Odd Samples at Steady State
            if (n == Nr)
                % After One Last Tip, and T1,T2 Relaxation and Beta Precession
                Mtip = Rx*M;
                E = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]); Ez = [0;0;M0*(1-exp(-TE/T1))];
                P = [cos(b * TE/TR),sin(b * TE/TR),0; -sin(b * TE/TR),cos(b * TE/TR),0; 0,0,1];
                M1 = P * E * Mtip + Ez;

                % Save Sample Point after Steady State is Reached
                Meven(1,SampleCount) = norm(M1(1:2));
                Meven(2,SampleCount) = angle(M1(1)+j*M1(2));
            elseif (n == Nr - 1)
                % After One Last Tip, and T1,T2 Relaxation and Beta Precession
                Mtip = Rx*M;
                E = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]); Ez = [0;0;M0*(1-exp(-TE/T1))];
                P = [cos(b * TE/TR),sin(b * TE/TR),0; -sin(b * TE/TR),cos(b * TE/TR),0; 0,0,1];
                M2 = P * E * Mtip + Ez;

                % Save Sample Point after Steady State is Reached
                Modd(1,SampleCount) = norm(M2(1:2));
                Modd(2,SampleCount) = angle(M2(1)+j*M2(2));
            end
        
        end
    end

    Beta = linspace(-BetaMax, BetaMax, Ns);
    f = Beta / TR / (2*pi);
    figure(1);
    PlotSpectra(f, Meven(1,:), Meven(2,:), plotstring);
    figure(2);
    PlotSpectra(f, Modd(1,:), Modd(2,:), plotstring);
    
    % Make a Column Vector
    Meven = MagAngleToComplex(Meven(1,:),Meven(2,:))';
    Modd = MagAngleToComplex(Modd(1,:),Modd(2,:))';
end

function PlotSpectra(f, Mag, Angle, plotstring)
    % Plot SSFP Spectrum
    subplot(2,1,1); plot(f, Mag, plotstring); hold on;
    ylabel('Magnitude');
    xlabel('Frequency (Hertz)');
    subplot(2,1,2); plot(f, Angle, plotstring); hold off;
    xlabel('Frequency (Hertz)');
    ylabel('Angle');   
end

function [Z] = MagAngleToComplex(Mag,Ang)
    Z = Mag .* exp(i.*Ang);
end

function [Mag, Ang] = ComplexToMagAng(Z)
    Mag = abs(Z);
    Ang = angle(Z);
end

function [Out] = MaxIntensity( X )
    Out = max(abs(X), [], 2);
end

function [Out] = ComplexSum( X )
    Out = sum(X, 2);
end

function [Out] = MagnitudeSum( X )
    Out = sum( abs(X), 2);
end

function [Out] = SumOfSquares( X ) 
    Out = sqrt( sum(abs(X).^2, 2) );
end 

function [Out] = RemoveNulls(M)
    Ns = length(M);      % Sample Length
    N = length(M(1,:));  % Number of Datasets
    
    Mout = [];
    for ns = 1:Ns
        min = M(ns, 1);
        minIndex = 1;
        for n = 1:N
            if abs(M(ns, n)) < abs(min)
                minIndex = n;
            end
        end
        
        sample = ones(N,1);
        sample(minIndex) = 0;         
        
        count = 1;
        for n = 1:N
           if (sample(n) == 1)
                Mout(ns, count) = M(ns, n);
                count = count + 1;
           end
        end
        
    end
    
    Out = Mout;
end

