% RF Pulse
x = linspace(-3,3,80);
RF = [3*sin(pi*x)./(pi*x) zeros(1,260)];

% Gz
x = linspace(0,1,20);
x2 = linspace(0,2,40);
Gz = [x ones(1,40) -x2+1 -ones(1,20) x-1 zeros(1,200)];

% Gx
x = linspace(0,1,20);
x2 = linspace(0,2,40);
Gx = [zeros(1,120) -x -ones(1,20) x2-1 ones(1,40) -x2+1 -ones(1,20) x-1 zeros(1,20)];

% Gy
x = linspace(0,1,20);
Gy = [zeros(1,120) x ones(1,20) -x+1 zeros(1,80) -x -ones(1,20) x-1 zeros(1,20)];
Gy2 = .66*Gy;
Gy3 = .33*Gy;
Gy4 = -.33*Gy;
Gy5 = -.66*Gy;
Gy6 = -Gy;

% Sample
x = linspace(0,1,20);
x2 = linspace(0,2,40);
Readout = [zeros(1,200) ones(1,40) zeros(1,100)];


figure(); hold on;
plot(RF+9,'k');
plot(Gz+6,'k');
plot(Gx+3,'k');
plot(Gy,'k-'); plot(Gy2,'k--'); plot(Gy3,'k--'); plot(Gy4,'k--'); plot(Gy5,'k--'); plot(Gy6,'k--'); 
plot(Readout-3,'k');
text(-60, 9, 'RF', 'FontSize', 18);
text(-60, 6, 'Gz', 'FontSize', 18);
text(-60, 3, 'Gx', 'FontSize', 18);
text(-60, 0, 'Gy', 'FontSize', 18);
text(-80, -3, 'Readout', 'FontSize', 18);
text(330, -3.5, 't', 'FontSize', 18);

hold off;
xlim([-90 350]);
ylim([-4 12]);
axis off;