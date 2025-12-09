%% SISTEM KONTROL LQR DAN LQG UNTUK STEAM TURBINE
% Author: Control System Design
% Date: 2025

clear all; close all; clc;

%% 1. DEFINISI SISTEM STEAM TURBINE
fprintf('=== SISTEM STEAM TURBINE ===\n\n');

% Parameter fisik
J = 1000;        % Momen inersia rotor [kg.m^2]
tau_p = 0.5;     % Time constant tekanan [s]
tau_T = 2.0;     % Time constant temperatur [s]
K_t = 50;        % Gain turbin
D = 10;          % Damping coefficient

% Matriks sistem state-space
% State: x = [delta_omega; delta_P; delta_T]
% Input: u = valve_opening
A = [-D/J,      K_t/J,        0;
      0,       -1/tau_p,      0.5/tau_p;
      0,        0,          -1/tau_T];

B = [0;
     1/tau_p;
     0.8/tau_T];

C = [1, 0, 0];  % Output: kecepatan rotor
D_sys = 0;

% Noise covariance matrices
G = [0.1; 0.5; 0.3];  % Process noise distribution
Q_noise = 0.1;         % Process noise covariance
R_noise = 0.01;        % Measurement noise covariance

% Buat sistem state-space
sys = ss(A, B, C, D_sys);

fprintf('Matriks A (Dinamika Sistem):\n');
disp(A);
fprintf('Matriks B (Input):\n');
disp(B);
fprintf('Matriks C (Output):\n');
disp(C);

% Cek controllability dan observability
Co = ctrb(A, B);
Ob = obsv(A, C);
rank_Co = rank(Co);
rank_Ob = rank(Ob);

fprintf('Rank Controllability Matrix: %d (sistem %s)\n', rank_Co, ...
    iif(rank_Co == size(A,1), 'CONTROLLABLE', 'NOT CONTROLLABLE'));
fprintf('Rank Observability Matrix: %d (sistem %s)\n\n', rank_Ob, ...
    iif(rank_Ob == size(A,1), 'OBSERVABLE', 'NOT OBSERVABLE'));

%% 2. DESAIN KONTROLER LQR
fprintf('=== DESAIN LQR CONTROLLER ===\n\n');

% Matriks pembobotan LQR
Q_lqr = diag([100, 10, 5]);  % Penalti state (prioritas kecepatan rotor)
R_lqr = 1;                    % Penalti kontrol

% Hitung gain LQR
[K_lqr, S_lqr, E_lqr] = lqr(A, B, Q_lqr, R_lqr);

fprintf('Gain LQR (K):\n');
disp(K_lqr);
fprintf('Closed-loop poles LQR:\n');
disp(E_lqr);

% Sistem closed-loop dengan LQR
A_cl_lqr = A - B*K_lqr;
sys_cl_lqr = ss(A_cl_lqr, B, C, D_sys);

%% 3. DESAIN KONTROLER LQG (Kalman Filter + LQR)
fprintf('\n=== DESAIN LQG CONTROLLER ===\n\n');

% Desain Kalman Filter (Estimator)
Q_kalman = G * Q_noise * G';  % Process noise covariance
R_kalman = R_noise;            % Measurement noise covariance

% Hitung Kalman gain
[K_kalman, P_kalman, E_kalman] = lqe(A, G, C, Q_noise, R_noise);

fprintf('Gain Kalman Filter (L):\n');
disp(K_kalman);
fprintf('Estimator poles:\n');
disp(E_kalman);

% Sistem LQG (Kontroler + Estimator)
% Augmented system untuk LQG
A_lqg = [A-B*K_lqr,           B*K_lqr;
         zeros(size(A)),      A-K_kalman*C];
B_lqg = [B; zeros(size(B))];
C_lqg = [C, zeros(size(C))];
D_lqg = 0;

sys_lqg = ss(A_lqg, B_lqg, C_lqg, D_lqg);

%% 4. ANALISIS KESTABILAN
fprintf('\n=== ANALISIS KESTABILAN ===\n\n');

% Open-loop poles
poles_ol = eig(A);
fprintf('Open-loop poles:\n');
disp(poles_ol);
fprintf('Sistem open-loop: %s\n', iif(all(real(poles_ol) < 0), 'STABIL', 'TIDAK STABIL'));

% LQR closed-loop poles
poles_lqr = eig(A_cl_lqr);
fprintf('\nLQR closed-loop poles:\n');
disp(poles_lqr);
fprintf('Sistem LQR: %s\n', iif(all(real(poles_lqr) < 0), 'STABIL', 'TIDAK STABIL'));

% Stability margins LQR
[Gm_lqr, Pm_lqr, Wcg_lqr, Wcp_lqr] = margin(sys_cl_lqr);
fprintf('Gain Margin LQR: %.2f dB\n', 20*log10(Gm_lqr));
fprintf('Phase Margin LQR: %.2f deg\n', Pm_lqr);

% LQG closed-loop poles (augmented system)
poles_lqg = eig(A_lqg);
fprintf('\nLQG closed-loop poles (Controller + Estimator):\n');
disp(poles_lqg);
fprintf('Sistem LQG: %s\n', iif(all(real(poles_lqg) < 0), 'STABIL', 'TIDAK STABIL'));

% Analisis separation principle verification
fprintf('\n--- Verification of Separation Principle ---\n');
fprintf('Controller poles (LQR):\n');
disp(poles_lqr);
fprintf('Estimator poles (Kalman):\n');
disp(E_kalman);
fprintf('Combined LQG poles should contain both sets above\n');

% Damping ratio dan natural frequency
[wn_lqr, zeta_lqr] = damp(sys_cl_lqr);
fprintf('\n--- LQR Damping Analysis ---\n');
for i = 1:length(wn_lqr)
    fprintf('Pole %d: wn = %.4f rad/s, zeta = %.4f\n', i, wn_lqr(i), zeta_lqr(i));
end

% Untuk LQG, analisis poles controller dan estimator terpisah
sys_estimator = ss(A-K_kalman*C, [B, K_kalman], eye(3), zeros(3,2));
[wn_est, zeta_est] = damp(sys_estimator);
fprintf('\n--- Kalman Filter (Estimator) Damping Analysis ---\n');
for i = 1:length(wn_est)
    fprintf('Pole %d: wn = %.4f rad/s, zeta = %.4f\n', i, wn_est(i), zeta_est(i));
end

% Stability margins untuk sistem LQG lengkap
% Buat transfer function dari reference ke output untuk LQG
sys_lqg_tf = ss(A_lqg, B_lqg, C_lqg, D_lqg);
[Gm_lqg, Pm_lqg, Wcg_lqg, Wcp_lqg] = margin(sys_lqg_tf);
fprintf('\n--- LQG Stability Margins ---\n');
fprintf('Gain Margin LQG: %.2f dB (at %.2f rad/s)\n', 20*log10(Gm_lqg), Wcg_lqg);
fprintf('Phase Margin LQG: %.2f deg (at %.2f rad/s)\n', Pm_lqg, Wcp_lqg);

% Kondisi stabilitas
fprintf('\n--- Kriteria Stabilitas ---\n');
fprintf('LQR:\n');
fprintf('  All poles in LHP: %s\n', iif(all(real(poles_lqr) < 0), 'YES ✓', 'NO ✗'));
fprintf('  Minimum damping ratio: %.4f %s\n', min(zeta_lqr), ...
    iif(min(zeta_lqr) > 0.05, '(Good)', '(Poor)'));
fprintf('  Gain margin > 6 dB: %s\n', iif(20*log10(Gm_lqr) > 6, 'YES ✓', 'NO ✗'));
fprintf('  Phase margin > 45 deg: %s\n', iif(Pm_lqr > 45, 'YES ✓', 'NO ✗'));

fprintf('\nLQG:\n');
fprintf('  All poles in LHP: %s\n', iif(all(real(poles_lqg) < 0), 'YES ✓', 'NO ✗'));
fprintf('  Estimator faster than controller: %s\n', ...
    iif(max(real(E_kalman)) < max(real(poles_lqr)), 'YES ✓', 'NO ✗'));
fprintf('  Gain margin > 6 dB: %s\n', iif(20*log10(Gm_lqg) > 6, 'YES ✓', 'NO ✗'));
fprintf('  Phase margin > 45 deg: %s\n', iif(Pm_lqg > 45, 'YES ✓', 'NO ✗'));

%% 5. SIMULASI RESPONS STEP
fprintf('\n=== SIMULASI RESPONS SISTEM ===\n\n');

t = 0:0.01:20;  % Waktu simulasi [s]
r = 10;         % Step reference [rad/s]

% Open-loop response
[y_ol, t_ol] = step(r*sys, t);

% LQR response (dengan reference tracking)
sys_cl_lqr_ref = ss(A-B*K_lqr, B*K_lqr(1), C, 0);
[y_lqr, t_lqr] = step(r*sys_cl_lqr_ref, t);

% Simulasi dengan gangguan untuk LQG
rng(42);  % Seed untuk reproducibility
u_ref = r * ones(size(t));
w = sqrt(Q_noise) * randn(size(t));  % Process noise
v = sqrt(R_noise) * randn(size(t));  % Measurement noise

% Simulasi LQR dengan noise
x_lqr_noise = zeros(3, length(t));
y_lqr_noise = zeros(1, length(t));
x_lqr_noise(:,1) = [0; 0; 0];

for i = 1:length(t)-1
    u_lqr = -K_lqr * (x_lqr_noise(:,i) - [r; 0; 0]);
    dx = A*x_lqr_noise(:,i) + B*u_lqr + G*w(i);
    x_lqr_noise(:,i+1) = x_lqr_noise(:,i) + dx*0.01;
    y_lqr_noise(i) = C*x_lqr_noise(:,i) + v(i);
end

% Simulasi LQG dengan noise
x_lqg = zeros(3, length(t));
x_hat = zeros(3, length(t));  % Estimated state
y_lqg = zeros(1, length(t));
x_lqg(:,1) = [0; 0; 0];
x_hat(:,1) = [0; 0; 0];

for i = 1:length(t)-1
    % Measurement
    y_meas = C*x_lqg(:,i) + v(i);
    
    % Control law (menggunakan estimated state)
    u_lqg = -K_lqr * (x_hat(:,i) - [r; 0; 0]);
    
    % Plant dynamics
    dx = A*x_lqg(:,i) + B*u_lqg + G*w(i);
    x_lqg(:,i+1) = x_lqg(:,i) + dx*0.01;
    
    % Kalman filter (estimator)
    y_hat = C*x_hat(:,i);
    dx_hat = A*x_hat(:,i) + B*u_lqg + K_kalman*(y_meas - y_hat);
    x_hat(:,i+1) = x_hat(:,i) + dx_hat*0.01;
    
    y_lqg(i) = C*x_lqg(:,i);
end

%% 6. ANALISIS PERFORMANSI
fprintf('=== ANALISIS PERFORMANSI ===\n\n');

% Step response metrics untuk LQR (tanpa noise)
info_lqr = stepinfo(y_lqr, t_lqr, r);
fprintf('LQR Performance:\n');
fprintf('  Rise Time: %.3f s\n', info_lqr.RiseTime);
fprintf('  Settling Time: %.3f s\n', info_lqr.SettlingTime);
fprintf('  Overshoot: %.2f %%\n', info_lqr.Overshoot);
fprintf('  Peak: %.3f rad/s\n', info_lqr.Peak);
fprintf('  Steady-state error: %.3f rad/s\n', abs(r - y_lqr(end)));

% Performance untuk LQG dengan noise
ss_error_lqg = abs(r - mean(y_lqg(end-100:end)));
std_lqg = std(y_lqg(end-100:end));
fprintf('\nLQG Performance (dengan noise):\n');
fprintf('  Steady-state error: %.3f rad/s\n', ss_error_lqg);
fprintf('  Output std deviation: %.3f rad/s\n', std_lqg);

% ISE dan IAE
ise_lqr_noise = sum((u_ref - y_lqr_noise).^2) * 0.01;
iae_lqr_noise = sum(abs(u_ref - y_lqr_noise)) * 0.01;
ise_lqg = sum((u_ref - y_lqg).^2) * 0.01;
iae_lqg = sum(abs(u_ref - y_lqg)) * 0.01;

fprintf('\nIntegral Performance Indices:\n');
fprintf('  ISE LQR (with noise): %.2f\n', ise_lqr_noise);
fprintf('  IAE LQR (with noise): %.2f\n', iae_lqr_noise);
fprintf('  ISE LQG: %.2f\n', ise_lqg);
fprintf('  IAE LQG: %.2f\n', iae_lqg);

% Estimation performance
fprintf('\nKalman Filter Estimation Performance:\n');
for i = 1:3
    est_error = x_lqg(i,1000:end) - x_hat(i,1000:end);
    rmse = sqrt(mean(est_error.^2));
    mae = mean(abs(est_error));
    fprintf('  State %d - RMSE: %.4f, MAE: %.4f\n', i, rmse, mae);
end

% Robustness metrics
fprintf('\nRobustness Metrics:\n');
fprintf('LQR:\n');
fprintf('  Gain Margin: %.2f dB (%.2fx)\n', 20*log10(Gm_lqr), Gm_lqr);
fprintf('  Phase Margin: %.2f deg\n', Pm_lqr);
fprintf('  Robustness: %s\n', iif(Gm_lqr > 2 && Pm_lqr > 45, 'EXCELLENT', ...
    iif(Gm_lqr > 1.5 && Pm_lqr > 30, 'GOOD', 'ADEQUATE')));

fprintf('\nLQG:\n');
fprintf('  Gain Margin: %.2f dB (%.2fx)\n', 20*log10(Gm_lqg), Gm_lqg);
fprintf('  Phase Margin: %.2f deg\n', Pm_lqg);
fprintf('  Robustness: %s\n', iif(Gm_lqg > 2 && Pm_lqg > 45, 'EXCELLENT', ...
    iif(Gm_lqg > 1.5 && Pm_lqg > 30, 'GOOD', 'ADEQUATE')));

% Bandwidth comparison
bw_ol = bandwidth(sys);
bw_lqr = bandwidth(sys_cl_lqr);
bw_lqg = bandwidth(sys_lqg_tf);
fprintf('\nBandwidth Analysis:\n');
fprintf('  Open-loop: %.4f rad/s\n', bw_ol);
fprintf('  LQR: %.4f rad/s (%.1fx faster)\n', bw_lqr, bw_lqr/bw_ol);
fprintf('  LQG: %.4f rad/s (%.1fx faster)\n', bw_lqg, bw_lqg/bw_ol);

%% 7. VISUALISASI
% Figure 1: Step Response Comparison
figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1);
plot(t_ol, y_ol, 'b-', 'LineWidth', 1.5); hold on;
plot(t_lqr, y_lqr, 'r-', 'LineWidth', 1.5);
yline(r, 'k--', 'LineWidth', 1);
grid on;
xlabel('Time [s]');
ylabel('Rotor Speed [rad/s]');
title('Step Response (Tanpa Noise)');
legend('Open-loop', 'LQR', 'Reference', 'Location', 'best');

subplot(1,3,2);
plot(t, y_lqr_noise, 'r-', 'LineWidth', 1); hold on;
plot(t, y_lqg, 'g-', 'LineWidth', 1);
yline(r, 'k--', 'LineWidth', 1);
grid on;
xlabel('Time [s]');
ylabel('Rotor Speed [rad/s]');
title('Response dengan Noise');
legend('LQR', 'LQG', 'Reference', 'Location', 'best');

subplot(1,3,3);
plot(t, x_lqg(1,:), 'g-', 'LineWidth', 1.5); hold on;
plot(t, x_hat(1,:), 'b--', 'LineWidth', 1.5);
yline(r, 'k--', 'LineWidth', 1);
grid on;
xlabel('Time [s]');
ylabel('Rotor Speed [rad/s]');
title('LQG: True State vs Estimated State');
legend('True State', 'Estimated State', 'Reference', 'Location', 'best');

% Figure 2: Pole-Zero Map
figure('Position', [100, 550, 1200, 800]);
subplot(2,3,1);
pzmap(sys);
title('Open-loop Poles');
grid on;

subplot(2,3,2);
pzmap(sys_cl_lqr);
title('LQR Closed-loop Poles');
grid on;

subplot(2,3,3);
pzmap(sys_lqg_tf);
title('LQG System Poles (Controller + Estimator)');
grid on;

subplot(2,3,4);
hold on;
plot(real(poles_ol), imag(poles_ol), 'bx', 'MarkerSize', 12, 'LineWidth', 2);
plot(real(poles_lqr), imag(poles_lqr), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
plot(real(E_kalman), imag(E_kalman), 'g^', 'MarkerSize', 12, 'LineWidth', 2);
xline(0, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Pole Comparison: OL vs LQR vs Estimator');
legend('Open-loop', 'LQR Controller', 'Kalman Estimator', 'Location', 'best');

subplot(2,3,5);
hold on;
plot(real(poles_lqr), imag(poles_lqr), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
plot(real(poles_lqg), imag(poles_lqg), 'md', 'MarkerSize', 10, 'LineWidth', 2);
xline(0, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('LQR vs LQG Poles');
legend('LQR', 'LQG (Combined)', 'Location', 'best');

subplot(2,3,6);
theta = linspace(0, 2*pi, 100);
hold on;
% Damping ratio lines
for zeta = [0.1, 0.3, 0.5, 0.7, 0.9]
    angle = acos(zeta);
    plot([0, -10*cos(angle)], [0, 10*sin(angle)], 'k:', 'LineWidth', 0.5);
    plot([0, -10*cos(angle)], [0, -10*sin(angle)], 'k:', 'LineWidth', 0.5);
end
plot(real(poles_ol), imag(poles_ol), 'bx', 'MarkerSize', 12, 'LineWidth', 2);
plot(real(poles_lqr), imag(poles_lqr), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
plot(real(E_kalman), imag(E_kalman), 'g^', 'MarkerSize', 12, 'LineWidth', 2);
xline(0, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Damping Ratio Visualization');
legend('Open-loop', 'LQR', 'Estimator', 'Location', 'best');
axis([-12 1 -8 8]);

% Figure 3: Bode Diagram
figure('Position', [100, 100, 1400, 900]);

% Bode magnitude dan phase
subplot(3,2,1);
bode(sys, sys_cl_lqr, sys_lqg_tf);
title('Bode Diagram Comparison');
legend('Open-loop', 'LQR', 'LQG', 'Location', 'best');
grid on;

% Nyquist plot
subplot(3,2,2);
nyquist(sys, sys_cl_lqr, sys_lqg_tf);
title('Nyquist Diagram Comparison');
legend('Open-loop', 'LQR', 'LQG', 'Location', 'best');
grid on;

% Sensitivity function
S_lqr = feedback(1, sys_cl_lqr*K_lqr(1));
S_lqg = feedback(1, sys_lqg_tf);
subplot(3,2,3);
bode(S_lqr, S_lqg);
title('Sensitivity Function (Disturbance Rejection)');
legend('LQR', 'LQG', 'Location', 'best');
grid on;

% Complementary sensitivity
T_lqr = feedback(sys_cl_lqr*K_lqr(1), 1);
T_lqg = feedback(sys_lqg_tf, 1);
subplot(3,2,4);
bode(T_lqr, T_lqg);
title('Complementary Sensitivity (Noise Attenuation)');
legend('LQR', 'LQG', 'Location', 'best');
grid on;

% Gain and Phase margins visualization
subplot(3,2,5);
margin(sys_cl_lqr);
title('LQR Gain and Phase Margins');
grid on;

subplot(3,2,6);
margin(sys_lqg_tf);
title('LQG Gain and Phase Margins');
grid on;

% Figure 4: State trajectories
figure('Position', [100, 100, 1400, 800]);
state_names = {'\Delta\omega [rad/s]', '\DeltaP [bar]', '\DeltaT [°C]'};
for i = 1:3
    subplot(3,2,2*i-1);
    plot(t, x_lqr_noise(i,:), 'r-', 'LineWidth', 1.5);
    grid on;
    ylabel(state_names{i});
    if i == 1
        title('LQR State Trajectories');
    end
    if i == 3
        xlabel('Time [s]');
    end
    
    subplot(3,2,2*i);
    plot(t, x_lqg(i,:), 'g-', 'LineWidth', 1.5); hold on;
    plot(t, x_hat(i,:), 'b--', 'LineWidth', 1.5);
    grid on;
    ylabel(state_names{i});
    if i == 1
        title('LQG State Trajectories');
        legend('True', 'Estimated', 'Location', 'best');
    end
    if i == 3
        xlabel('Time [s]');
    end
end

% Figure 5: Estimation Error Analysis
figure('Position', [100, 100, 1400, 600]);
for i = 1:3
    subplot(2,3,i);
    estimation_error = x_lqg(i,:) - x_hat(i,:);
    plot(t, estimation_error, 'b-', 'LineWidth', 1);
    grid on;
    xlabel('Time [s]');
    ylabel(['Error ' state_names{i}]);
    title(['State ' num2str(i) ' Estimation Error']);
    
    subplot(2,3,i+3);
    histogram(estimation_error(1000:end), 30, 'Normalization', 'pdf');
    hold on;
    % Overlay normal distribution
    mu_err = mean(estimation_error(1000:end));
    sigma_err = std(estimation_error(1000:end));
    x_hist = linspace(min(estimation_error), max(estimation_error), 100);
    y_hist = normpdf(x_hist, mu_err, sigma_err);
    plot(x_hist, y_hist, 'r-', 'LineWidth', 2);
    grid on;
    xlabel(['Error ' state_names{i}]);
    ylabel('Probability Density');
    title(['Error Distribution (μ=' num2str(mu_err, '%.3f') ', σ=' num2str(sigma_err, '%.3f') ')']);
    legend('Actual', 'Gaussian Fit', 'Location', 'best');
end

% Figure 6: Robustness Analysis - Root Locus
figure('Position', [100, 100, 1400, 500]);
subplot(1,2,1);
rlocus(sys);
title('Root Locus - Open Loop System');
grid on;

subplot(1,2,2);
% Root locus untuk LQR dengan variasi gain
sys_ol_lqr = ss(A, B*K_lqr(1), C, 0);
rlocus(sys_ol_lqr);
title('Root Locus - LQR Loop Transfer Function');
grid on;

%% 8. ANALISIS KESTABILAN LANJUTAN

fprintf('\n=== ANALISIS KESTABILAN LANJUTAN ===\n\n');

% 8.1 ROUTH-HURWITZ CRITERION
fprintf('--- ROUTH-HURWITZ STABILITY CRITERION ---\n');

% Untuk LQR Closed-loop
fprintf('\nLQR Closed-Loop System:\n');
char_poly_lqr = poly(A_cl_lqr);
fprintf('Characteristic Polynomial:\n');
fprintf('s^3 + %.4f*s^2 + %.4f*s + %.4f\n', -char_poly_lqr(2), char_poly_lqr(3), -char_poly_lqr(4));

% Routh array untuk sistem orde 3
a0 = char_poly_lqr(1);
a1 = char_poly_lqr(2);
a2 = char_poly_lqr(3);
a3 = char_poly_lqr(4);

% Routh array construction
routh = zeros(4, 2);
routh(1, 1) = a0;
routh(1, 2) = a2;
routh(2, 1) = a1;
routh(2, 2) = a3;
routh(3, 1) = (routh(2,1)*routh(1,2) - routh(1,1)*routh(2,2)) / routh(2,1);
routh(3, 2) = 0;
routh(4, 1) = (routh(3,1)*routh(2,2) - routh(2,1)*0) / routh(3,1);
routh(4, 2) = 0;

fprintf('\nRouth Array:\n');
fprintf('s^3 | %.4f  %.4f\n', routh(1,1), routh(1,2));
fprintf('s^2 | %.4f  %.4f\n', routh(2,1), routh(2,2));
fprintf('s^1 | %.4f  %.4f\n', routh(3,1), routh(3,2));
fprintf('s^0 | %.4f  %.4f\n', routh(4,1), routh(4,2));

% Check sign changes
first_column = routh(:, 1);
sign_changes = sum(diff(sign(first_column)) ~= 0);
fprintf('\nSign changes in first column: %d\n', sign_changes);
fprintf('System is %s (No sign changes = Stable)\n', ...
    iif(sign_changes == 0, 'STABLE', 'UNSTABLE'));

% Untuk LQG System
fprintf('\n\nLQG System:\n');
char_poly_lqg = poly(A_lqg);
fprintf('Characteristic Polynomial (Order 6):\n');
fprintf('Coefficients: ');
fprintf('%.4f ', char_poly_lqg);
fprintf('\n');

% For higher order system (order 6)
n = length(char_poly_lqg) - 1;  % n = 6
num_cols = ceil((n+1)/2);  % Number of columns needed

routh_lqg = zeros(n+1, num_cols);

% Fill first two rows
routh_lqg(1, :) = [char_poly_lqg(1:2:end), zeros(1, num_cols - length(char_poly_lqg(1:2:end)))];
routh_lqg(2, :) = [char_poly_lqg(2:2:end), zeros(1, num_cols - length(char_poly_lqg(2:2:end)))];

% Calculate remaining rows
for i = 3:n+1
    for j = 1:num_cols-1
        if routh_lqg(i-1, 1) ~= 0
            routh_lqg(i, j) = (routh_lqg(i-1,1)*routh_lqg(i-2,j+1) - ...
                               routh_lqg(i-2,1)*routh_lqg(i-1,j+1)) / routh_lqg(i-1,1);
        else
            % Handle zero in first column (replace with small epsilon)
            epsilon = 1e-6;
            routh_lqg(i-1, 1) = epsilon;
            routh_lqg(i, j) = (routh_lqg(i-1,1)*routh_lqg(i-2,j+1) - ...
                               routh_lqg(i-2,1)*routh_lqg(i-1,j+1)) / routh_lqg(i-1,1);
        end
    end
end

fprintf('\nRouth Array (first column):\n');
for i = 1:n+1
    fprintf('s^%d | %.4f\n', n+1-i, routh_lqg(i, 1));
end

sign_changes_lqg = sum(diff(sign(routh_lqg(:,1))) ~= 0);
fprintf('\nSign changes in first column: %d\n', sign_changes_lqg);
fprintf('System is %s\n', iif(sign_changes_lqg == 0, 'STABLE', 'UNSTABLE'));

% 8.2 NYQUIST STABILITY CRITERION
fprintf('\n--- NYQUIST STABILITY CRITERION ---\n');

% Untuk LQR
[re_lqr, im_lqr, w_nyq] = nyquist(sys_cl_lqr*K_lqr(1), logspace(-2, 2, 500));
re_lqr = squeeze(re_lqr);
im_lqr = squeeze(im_lqr);

% Count encirclements of -1 point
encirclements_lqr = 0;
for i = 1:length(re_lqr)-1
    % Simple encirclement check
    if (re_lqr(i) < -1 && im_lqr(i) > 0 && im_lqr(i+1) < 0) || ...
       (re_lqr(i) < -1 && im_lqr(i) < 0 && im_lqr(i+1) > 0)
        encirclements_lqr = encirclements_lqg + 1;
    end
end

fprintf('LQR System:\n');
fprintf('  Open-loop unstable poles (P): %d\n', sum(real(poles_ol) > 0));
fprintf('  Encirclements of -1 point (N): %d\n', encirclements_lqr);
fprintf('  Closed-loop unstable poles (Z = N + P): %d\n', encirclements_lqr + sum(real(poles_ol) > 0));
fprintf('  System is %s (Z = 0 for stability)\n', ...
    iif(encirclements_lqr + sum(real(poles_ol) > 0) == 0, 'STABLE', 'UNSTABLE'));

fprintf('\nLQG System:\n');
fprintf('  Open-loop unstable poles (P): %d\n', sum(real(poles_ol) > 0));
fprintf('  Closed-loop unstable poles: %d\n', sum(real(poles_lqg) > 0));
fprintf('  System is %s\n', iif(sum(real(poles_lqg) > 0) == 0, 'STABLE', 'UNSTABLE'));

% 8.3 BODE STABILITY ANALYSIS
fprintf('\n--- BODE STABILITY ANALYSIS ---\n');

fprintf('LQR System:\n');
% Handle infinite margins
if isinf(Gm_lqr)
    fprintf('  Gain Margin: Inf dB (Infinite - Excellent)\n');
else
    fprintf('  Gain Margin: %.2f dB (> 6 dB required)\n', 20*log10(Gm_lqr));
end

if isinf(Pm_lqr)
    fprintf('  Phase Margin: Inf deg (Infinite - Excellent)\n');
else
    fprintf('  Phase Margin: %.2f deg (> 45 deg required)\n', Pm_lqr);
end

if ~isinf(Wcg_lqr) && Wcg_lqr > 0
    fprintf('  Gain Crossover Frequency: %.4f rad/s\n', Wcg_lqr);
end
if ~isinf(Wcp_lqr) && Wcp_lqr > 0
    fprintf('  Phase Crossover Frequency: %.4f rad/s\n', Wcp_lqr);
end

fprintf('  Stability Assessment: %s\n', ...
    iif(isinf(Gm_lqr) || (20*log10(Gm_lqr) > 6 && Pm_lqr > 45), 'ROBUST', ...
    iif(20*log10(Gm_lqr) > 3 && Pm_lqr > 30, 'ADEQUATE', 'MARGINAL')));

fprintf('\nLQG System:\n');
if isinf(Gm_lqg)
    fprintf('  Gain Margin: Inf dB (Infinite - Excellent)\n');
else
    fprintf('  Gain Margin: %.2f dB\n', 20*log10(Gm_lqg));
end

if isinf(Pm_lqg)
    fprintf('  Phase Margin: Inf deg (Infinite - Excellent)\n');
else
    fprintf('  Phase Margin: %.2f deg\n', Pm_lqg);
end

if ~isinf(Wcg_lqg) && Wcg_lqg > 0
    fprintf('  Gain Crossover Frequency: %.4f rad/s\n', Wcg_lqg);
end
if ~isinf(Wcp_lqg) && Wcp_lqg > 0
    fprintf('  Phase Crossover Frequency: %.4f rad/s\n', Wcp_lqg);
end

fprintf('  Stability Assessment: %s\n', ...
    iif(isinf(Gm_lqg) || (20*log10(Gm_lqg) > 6 && Pm_lqg > 45), 'ROBUST', ...
    iif(20*log10(Gm_lqg) > 3 && Pm_lqg > 30, 'ADEQUATE', 'MARGINAL')));

% 8.4 LYAPUNOV STABILITY ANALYSIS
fprintf('\n--- LYAPUNOV STABILITY ANALYSIS ---\n');

% Untuk LQR: Lyapunov equation A'P + PA = -Q
% Jika P > 0 (positive definite), sistem stabil
Q_lyap = eye(3);  % Identity matrix
P_lyap_lqr = lyap(A_cl_lqr', Q_lyap);

fprintf('LQR System:\n');
fprintf('Lyapunov Matrix P (from A''P + PA + Q = 0):\n');
disp(P_lyap_lqr);

% Check positive definiteness
eig_P_lqr = eig(P_lyap_lqr);
fprintf('Eigenvalues of P: ');
fprintf('%.4f ', eig_P_lqr);
fprintf('\n');
fprintf('P is positive definite: %s\n', iif(all(eig_P_lqr > 0), 'YES', 'NO'));
fprintf('Lyapunov Stability: %s\n', iif(all(eig_P_lqr > 0), 'ASYMPTOTICALLY STABLE', 'NOT STABLE'));

% Lyapunov function: V(x) = x'Px
fprintf('\nLyapunov Function: V(x) = x''Px\n');
fprintf('  V(x) > 0 for all x ≠ 0: %s\n', iif(all(eig_P_lqr > 0), 'YES ✓', 'NO ✗'));
fprintf('  V̇(x) = x''(A''P + PA)x = -x''Qx < 0: YES ✓\n');
fprintf('  Therefore, system is asymptotically stable at origin\n');

% Untuk LQG System
fprintf('\n\nLQG System:\n');
Q_lyap_lqg = eye(6);
P_lyap_lqg = lyap(A_lqg', Q_lyap_lqg);

eig_P_lqg = eig(P_lyap_lqg);
fprintf('Eigenvalues of Lyapunov matrix P: ');
fprintf('%.4f ', eig_P_lqg);
fprintf('\n');
fprintf('P is positive definite: %s\n', iif(all(eig_P_lqg > 0), 'YES', 'NO'));
fprintf('Lyapunov Stability: %s\n', iif(all(eig_P_lqg > 0), 'ASYMPTOTICALLY STABLE', 'NOT STABLE'));

%% 9. ANALISIS PERFORMANSI DETAIL

fprintf('\n=== ANALISIS PERFORMANSI DETAIL ===\n\n');

% 9.1 TIME RESPONSE ANALYSIS
fprintf('--- TIME RESPONSE METRICS ---\n\n');

% Open-loop performance
info_ol = stepinfo(y_ol, t_ol, r);
fprintf('OPEN-LOOP System:\n');
fprintf('  Rise Time (10%%-90%%): %.3f s\n', info_ol.RiseTime);
fprintf('  Settling Time (2%% band): %.3f s\n', info_ol.SettlingTime);
fprintf('  Peak Time: %.3f s\n', info_ol.PeakTime);
fprintf('  Overshoot: %.2f %%\n', info_ol.Overshoot);
fprintf('  Peak Value: %.3f rad/s\n', info_ol.Peak);
fprintf('  Steady-State Value: %.3f rad/s\n', y_ol(end));
fprintf('  Steady-State Error: %.3f rad/s (%.2f %%)\n', ...
    abs(r - y_ol(end)), abs(r - y_ol(end))/r*100);

% LQR performance (without noise)
fprintf('\n\nLQR CLOSED-LOOP System (Ideal):\n');
fprintf('  Rise Time (10%%-90%%): %.3f s\n', info_lqr.RiseTime);
fprintf('  Settling Time (2%% band): %.3f s\n', info_lqr.SettlingTime);
fprintf('  Peak Time: %.3f s\n', info_lqr.PeakTime);
fprintf('  Overshoot: %.2f %%\n', info_lqr.Overshoot);
fprintf('  Peak Value: %.3f rad/s\n', info_lqr.Peak);
fprintf('  Steady-State Value: %.3f rad/s\n', y_lqr(end));
fprintf('  Steady-State Error: %.3f rad/s (%.2f %%)\n', ...
    abs(r - y_lqr(end)), abs(r - y_lqr(end))/r*100);

% Performance improvement
fprintf('\n\nPerformance Improvement (LQR vs Open-Loop):\n');
fprintf('  Rise Time: %.1f%% faster\n', ...
    (info_ol.RiseTime - info_lqr.RiseTime)/info_ol.RiseTime*100);
fprintf('  Settling Time: %.1f%% faster\n', ...
    (info_ol.SettlingTime - info_lqr.SettlingTime)/info_ol.SettlingTime*100);
fprintf('  Overshoot Reduction: %.2f%% → %.2f%%\n', ...
    info_ol.Overshoot, info_lqr.Overshoot);
fprintf('  Steady-State Error Reduction: %.2f%% → %.2f%%\n', ...
    abs(r - y_ol(end))/r*100, abs(r - y_lqr(end))/r*100);

% LQR with noise
fprintf('\n\nLQR System (with Noise):\n');
ss_val_lqr = mean(y_lqr_noise(end-100:end));
ss_error_lqr_noise = abs(r - ss_val_lqr);
fprintf('  Steady-State Value: %.3f ± %.3f rad/s\n', ss_val_lqr, std(y_lqr_noise(end-100:end)));
fprintf('  Steady-State Error: %.3f rad/s (%.2f %%)\n', ...
    ss_error_lqr_noise, ss_error_lqr_noise/r*100);
fprintf('  Output Variance: %.4f\n', var(y_lqr_noise(end-100:end)));

% LQG with noise
fprintf('\n\nLQG System (with Noise):\n');
ss_val_lqg = mean(y_lqg(end-100:end));
fprintf('  Steady-State Value: %.3f ± %.3f rad/s\n', ss_val_lqg, std_lqg);
fprintf('  Steady-State Error: %.3f rad/s (%.2f %%)\n', ...
    ss_error_lqg, ss_error_lqg/r*100);
fprintf('  Output Variance: %.4f\n', var(y_lqg(end-100:end)));

% Noise rejection comparison
fprintf('\n\nNoise Rejection Performance (LQR vs LQG):\n');
fprintf('  LQR Output Std Dev: %.4f rad/s\n', std(y_lqr_noise(end-100:end)));
fprintf('  LQG Output Std Dev: %.4f rad/s\n', std_lqg);
fprintf('  Noise Reduction by LQG: %.1f%%\n', ...
    (std(y_lqr_noise(end-100:end)) - std_lqg)/std(y_lqr_noise(end-100:end))*100);

% 9.2 TRANSIENT RESPONSE CHARACTERISTICS
fprintf('\n--- TRANSIENT RESPONSE CHARACTERISTICS ---\n\n');

% Calculate time constants from dominant poles
dominant_pole_ol = max(real(poles_ol));
dominant_pole_lqr = max(real(poles_lqr));

fprintf('Time Constants (τ = -1/Re(λ_dominant)):\n');
fprintf('  Open-Loop: τ = %.3f s\n', -1/dominant_pole_ol);
fprintf('  LQR: τ = %.3f s\n', -1/dominant_pole_lqr);
fprintf('  Speed-up factor: %.1fx\n', (-1/dominant_pole_ol)/(-1/dominant_pole_lqr));

% Damping characteristics
fprintf('\nDamping Characteristics:\n');
fprintf('  LQR Minimum Damping Ratio: %.4f\n', min(zeta_lqr));
fprintf('  LQR Maximum Damping Ratio: %.4f\n', max(zeta_lqr));
fprintf('  LQR Average Damping Ratio: %.4f\n', mean(zeta_lqr));
fprintf('  Assessment: %s\n', ...
    iif(min(zeta_lqr) > 0.7, 'OVERDAMPED (Slow, No Overshoot)', ...
    iif(min(zeta_lqr) > 0.4, 'WELL-DAMPED (Good Balance)', ...
    iif(min(zeta_lqr) > 0.1, 'UNDERDAMPED (Fast, Some Overshoot)', 'POORLY DAMPED'))));

% Natural frequencies
fprintf('\nNatural Frequencies:\n');
fprintf('  LQR Lowest: %.4f rad/s\n', min(wn_lqr));
fprintf('  LQR Highest: %.4f rad/s\n', max(wn_lqr));

% 9.3 PERFORMANCE SPECIFICATIONS CHECK
fprintf('\n--- PERFORMANCE SPECIFICATIONS COMPLIANCE ---\n\n');

spec_settling = 5.0;  % Required settling time < 5s
spec_overshoot = 15.0;  % Required overshoot < 15%
spec_ss_error = 0.5;  % Required steady-state error < 0.5%

fprintf('Specification Requirements:\n');
fprintf('  Settling Time < %.1f s\n', spec_settling);
fprintf('  Overshoot < %.1f %%\n', spec_overshoot);
fprintf('  Steady-State Error < %.1f %%\n', spec_ss_error);

fprintf('\nLQR Compliance:\n');
fprintf('  Settling Time: %.3f s - %s\n', info_lqr.SettlingTime, ...
    iif(info_lqr.SettlingTime < spec_settling, 'PASS ✓', 'FAIL ✗'));
fprintf('  Overshoot: %.2f %% - %s\n', info_lqr.Overshoot, ...
    iif(info_lqr.Overshoot < spec_overshoot, 'PASS ✓', 'FAIL ✗'));
fprintf('  Steady-State Error: %.2f %% - %s\n', abs(r - y_lqr(end))/r*100, ...
    iif(abs(r - y_lqr(end))/r*100 < spec_ss_error, 'PASS ✓', 'FAIL ✗'));

fprintf('\nLQG Compliance (with noise):\n');
fprintf('  Steady-State Error: %.2f %% - %s\n', ss_error_lqg/r*100, ...
    iif(ss_error_lqg/r*100 < spec_ss_error, 'PASS ✓', 'FAIL ✗'));
fprintf('  Noise Level: %.4f rad/s - %s\n', std_lqg, ...
    iif(std_lqg < 0.5, 'ACCEPTABLE ✓', 'HIGH ⚠'));

fprintf('\n=== SIMULASI SELESAI ===\n');

% Helper function
function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end