function [ Mc ] = SSFP(beta, M0, alpha, phi, dphi, TR, TE, T1, T2, Nr)
    % SSFP(beta) generates an SSFP transverse magnetization vector, Mxy,
    % with an off-resonance precession phase in TR, beta.
    %
    % SSFP(beta, M0, alpha, phi, dphi, TR, TE, T1, T2, Nr) generates an
    % SSFP transverse magnetization vector, Mxy, with the following
    % parameters:
    %   beta = additional off-resonance precession phase between excitations
    %   M0 = initial magnetization strength
    %   alpha = RF pulse tip in degrees
    %   phi = RF pulse phase in degrees
    %   dphi = RF pulse phase incrementation per TR - For Phase Cycling
    %   TR = repetition time
    %   TE = echo time
    %   T1 = longitudal relaxation time constant
    %   T2 = transverse relaxation time constant
    %   Nr = number of excitations

% Parameters
if not(exist('beta'))
    beta = 0;
end

if not(exist('M0') && exist('alpha') && exist('phi') && exist('dphi') && exist('TR') && exist('TE') && exist('T1') && exist('T2') && exist('Nr'))
    M0 = 1; 
    alpha = pi/2; phi = 0; dphi = 0;
    TR = 10/1000; TE = TR/2;
    T1 = 790/1000; T2 = 92/1000; 
    Nr = 200;
end

% Set Up Magnetization Vector
M = [0,0,M0]';  % Starting M Vector
Mc = [];

% Rotation Matrices
Rx = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;      % Rx( alpha, phi )
      (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
      sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];

%%
% Interative Method SSFP Simulation
if 1
    M = [0,0,M0]';
    for n = 1:Nr  % Interate through tips
        % Alpha Degree Tip
        Mtip = Rx*M;  

        % T1, T2 Relaxation
        E = diag([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)]);
        Ez = [0;0;M0*(1-exp(-TR/T1))];
        % Off Resonance Precesion

        P = [cos(beta),sin(beta),0; -sin(beta),cos(beta),0; 0,0,1];
        % Single after Relaxation and Off-Resonace Precession
        M = P * E * Mtip + Ez;

        phi = phi + dphi;
        Rx = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;      % Rx( alpha, phi )
              (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
              sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];
    end

    % After One Last Tip, and T1,T2 Relaxation and betaeta Precession
    Mtip = Rx*M;
    E = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]); Ez = [0;0;M0*(1-exp(-TE/T1))];
    P = [cos(beta * TE/TR),sin(beta * TE/TR),0; -sin(beta * TE/TR),cos(beta * TE/TR),0; 0,0,1];
    M = P * E * Mtip + Ez;

    % Save Sample Point after Steady State is Reached
    Mc = M(1) + 1j*M(2);

%%
% Steady State SSFP Simulation
else

    % Inital Magnetization Vector
    M0 = [0,0,1]'; I = eye(3);
    % T1, T2 Relaxation
    E = diag([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)]);
    Ete = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]);
    % Off Resonance Precesion
    P = [cos(beta),sin(beta),0; -sin(beta),cos(beta),0; 0,0,1];
    Pte = [cos(beta*TE/TR), sin(beta*TE/TR), 0; -sin(beta*TE/TR), cos(beta*TE/TR), 0; 0, 0, 1];
    
    Mminus = (I-P*E*Rx)^(-1)*(I-E)*M0;                         % Before Last Tip
    Mplus = Rx * Mminus;                                       % Singal After Tip
    Mte = Pte*Ete*Mplus + (I-Ete)*M0;                          % After Last Relaxation and Precession
    
    % Save Sample Point after Steady State Magnetization
    Mc = M(1) + 1j*M(2);
end


end

