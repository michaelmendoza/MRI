function [Meven, Modd] = FEMR(alpha, dphi1, dphi2, T1, T2, TR)

    TE = TR/2;
    Nr = 200;       % Number of Repetitions
    Ns = 200;       % Number of Samples for Off-Resonance Graph
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
        
        SampleCount = SampleCount + 1;
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
                Meven(SampleCount) = M1(1) + 1j * M1(2);
            elseif (n == Nr - 1)
                % After One Last Tip, and T1,T2 Relaxation and Beta Precession
                Mtip = Rx*M;
                E = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]); Ez = [0;0;M0*(1-exp(-TE/T1))];
                P = [cos(b * TE/TR),sin(b * TE/TR),0; -sin(b * TE/TR),cos(b * TE/TR),0; 0,0,1];
                M2 = P * E * Mtip + Ez;

                % Save Sample Point after Steady State is Reached
                Modd(SampleCount) = M2(1) + 1j * M2(2);
            end
        
        end
    end
end

