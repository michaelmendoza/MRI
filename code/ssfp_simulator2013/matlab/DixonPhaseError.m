
mw = 1; mf = 0;
theta = linspace(0, 2*pi, 100);
m1 = mw + mf;
m2 = (mw - mf) .* exp(1i * theta);
m3 = (mw + mf) .* exp(2i * theta);

mw_hat = (1/2) * (m1 + m2);
mf_hat = (1/2) * (m1 - m2);
mw_error = abs(mw - mw_hat);
mf_error = abs(mf - mf_hat);

phi_hat = (1/2) * angle(conj(m1) .* m3);
mw_hat3 = (1/2) * (m1 + m2 .* exp(-1i * phi_hat));
mf_hat3 = (1/2) * (m1 - m2 .* exp(-1i * phi_hat));
mw_error3 = abs(mw - mw_hat3);
mf_error3 = abs(mf - mf_hat3);

figure();
subplot(1,2,1);
plot(theta, abs(mw_hat), '-', theta, abs(mf_hat), '--');
set(gca,'FontSize',14)
legend('Water Estimate', 'Fat Estimate','FontSize', 14);
title('Component Estimates for Water Voxel using 2-Point Dixon','FontSize', 14);
xlabel('Phase Error (Radians)','FontSize', 14);
ylabel('Amplitude','FontSize', 14);
xlim([0 2*pi]); ylim([-0.1, 1.1]);

subplot(1,2,2);
plot(theta, abs(mw_hat3), '-', theta, abs(mf_hat3), '--');
set(gca,'FontSize',14)
legend('Water Estimate', 'Fat Estimate','FontSize', 14);
title('Component Estimates for Water Voxel using 3-Point Dixon','FontSize', 14);
xlabel('Phase Error (Radians)','FontSize', 14);
ylabel('Amplitude','FontSize', 14);
xlim([0 2*pi]); ylim([-0.1, 1.1]);

figure();
subplot(1,2,1);
plot(theta, angle(m3), '-', theta, phi_hat, '--');
set(gca,'FontSize',14)
legend('Phase Between 1 and 3 images', 'Phase Estimate Between 1 and 2 image','FontSize', 14);
title('Phase Error Estimation','FontSize', 14);
ylabel('Phase (Radians)','FontSize', 14);
xlabel('Phase (Radians)','FontSize', 14);
xlim([0 2*pi]);

subplot(1,2,2);
plot(theta, mw_error, '-', theta, mw_error3, '--');
set(gca,'FontSize',14)
legend('2-Point Dixon','3-Point Dixon','FontSize', 14);
title('Water Component Estimate Error','FontSize', 14);
xlabel('Phase (Radians)','FontSize', 14);
ylabel('Factional Error','FontSize', 14);
xlim([0 2*pi]); ylim([-0.1, 1.1]);



