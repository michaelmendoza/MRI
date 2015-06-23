% SSFP Simulation
% Michael Mendoza

clear;

% Parameters
alpha = pi/2; phi = 0; dphi = 0;
T1 = 790/1000; T2 = 92/1000; M0 = 1; 
% T1 = 270/1000; T2 = 85/1000; M0 = 1; 
TR = 10/1000; TE = TR/2;
Nr = 200;  % Number of Repetitions
Ns = 100; % Number of Samples for Off-Resonance Graph
BetaMax = pi; % Off Resonance - Free-Precession Angle (rad)

% Set Up Magnetization Vector
M = [0,0,M0]';  % Starting M Vector
Ms = [];

% Matrices
% Rotation Matrices
% Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];                                    % Rx( alpha ) with phi = 0
Rx = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;      % Rx( alpha, phi )
      (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
      sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];

%%
% Interative Method SSFP Simulation
if 1
count = 0;
for b = linspace(-BetaMax, BetaMax, Ns)   % Interate through Beta values
    
    M = [0,0,M0]';
    count = count + 1;
    for n = 1:Nr  % Interate through tips
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
        Rx = [cos(alpha)*sin(phi)^2+cos(phi)^2 (1-cos(alpha))*cos(phi)*sin(phi) -sin(alpha)*sin(phi) ;      % Rx( alpha, phi )
              (1-cos(alpha))*cos(phi)*sin(phi) cos(alpha)*cos(phi)^2+sin(phi)^2 sin(alpha)*cos(phi)  ;
              sin(alpha)*sin(phi)              -sin(alpha)*cos(phi)             cos(alpha)          ];
    end
    
    % After One Last Tip, and T1,T2 Relaxation and Beta Precession
    Mtip = Rx*M;
    E = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]); Ez = [0;0;M0*(1-exp(-TE/T1))];
    P = [cos(b * TE/TR),sin(b * TE/TR),0; -sin(b * TE/TR),cos(b * TE/TR),0; 0,0,1];
    M = P * E * Mtip + Ez;
    
    % Save Sample Point after Steady State is Reached
    Ms(1,count) = norm(M(1:2));
    Ms(2,count) = angle(M(1)+j*M(2));
    
end

%%
% Steady State SSFP Simulation
else
i = 0;
for b = linspace(0, BetaMax, Ns)
    i = i + 1;
    % Inital Magnetization Vector
    M0 = [0,0,1]'; I = eye(3);
    % T1, T2 Relaxation
    E = diag([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)]);
    Ete = diag([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)]);
    % Off Resonance Precesion
    P = [cos(b),sin(b),0; -sin(b),cos(b),0; 0,0,1];
    Pte = [cos(b*TE/TR), sin(b*TE/TR), 0; -sin(b*TE/TR), cos(b*TE/TR), 0; 0, 0, 1];
    
    Mminus = (I-P*E*Rx)^(-1)*(I-E)*M0;                         % Before Last Tip
    Mplus = Rx * Mminus;                                       % Singal After Tip
    Mte = Pte*Ete*Mplus + (I-Ete)*M0;                          % After Last Relaxation and Precession
    
    % Save Sample Point after Steady State Magnetization
    Ms(1,i) = norm(Mte(1:2));
    Ms(2,i) = angle(Mte(1)+Mte(2)*j);
end
end

%%
% Plot SSFP Spectrum
figure();
Beta = linspace(0, BetaMax, Ns);
f = Beta / TR / (2*pi);
subplot(2,1,1); plot(f, Ms(1,:));
ylabel('Magnitude');
xlabel('Frequency (Hertz)');
subplot(2,1,2); plot(f, Ms(2,:));
xlabel('Frequency (Hertz)');
ylabel('Angle');


