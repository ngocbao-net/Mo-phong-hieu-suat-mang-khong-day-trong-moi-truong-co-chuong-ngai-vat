%% ================== PRE-RUN ==================
clear; close all;
warning('off','all')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ================== OUTPUTS ==================
out_dir = fullfile(pwd, 'outputs_clean');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ================== THAM SỐ THAM CHIẾU ==================
lambda_ref = 0.00075;
Lmax_ref   = 20;
R          = 150;
N_data     = 80;
K          = 3;
N_mesh     = 12;

%% ================== THAM SỐ VẬT LÝ (3GPP TR 38.901) ==================
c          = 3e8;
k_B        = 1.38e-23;
T_K        = 290;
B          = 10e6;
NF_dB      = 9;
P_tx_dBm   = 30;
G_ant_dBi  = 5;
freq_GHz   = 2.0;
freq       = freq_GHz * 1e9;

N0_W      = k_B * T_K * B;
Noise_dBm = 10*log10(N0_W) + 30 + NF_dB;

h_BS = 25; h_UT = 1.5; h_E = 1.0;
h_BS_eff = h_BS - h_E;
h_UT_eff = h_UT - h_E;
d_BP     = 4 * h_BS_eff * h_UT_eff * freq / c;

fprintf('=== THAM SO 3GPP TR 38.901 ===\n');
fprintf('Noise floor: %.2f dBm\n', Noise_dBm);
fprintf('d_BP = %.1f m\n', d_BP);

%% ================== HÀM PATH LOSS ==================
PL_LOS_func = @(d) (d < d_BP) .* (28.0 + 22*log10(max(d,1e-3)) + 20*log10(freq_GHz)) + ...
                   (d >= d_BP).* (28.0 + 40*log10(max(d,1e-3)) + 20*log10(freq_GHz) ...
                                  - 9*log10(d_BP^2 + (h_BS - h_UT)^2));
PL_NLOS_func = @(d) max(PL_LOS_func(d), ...
                        13.54 + 39.08*log10(max(d,1e-3)) + 20*log10(freq_GHz) - 0.6*(h_UT - 1.5));

%% ================== MỞ FILE LOG ==================
log_file  = fopen('log_results.txt', 'a');
csv_file  = fopen('bang_ket_qua.csv', 'a');
ngay_chay = datestr(now, 'dd/mm/yyyy HH:MM:SS');

fseek(csv_file, 0, 'eof');
if ftell(csv_file) == 0
    fprintf(csv_file, ['Ngay_chay,Phan,N_data,K,lambda_scene,Lmax_scene,R,freq_GHz,B_MHz,P_tx_dBm,' ...
        'NF_dB,Noise_dBm,d_BP_m,N_run,' ...
        'Acc_kNN_pct,Std_kNN_pct,Acc_kNMAP_pct,Std_kNMAP_pct,' ...
        'DiemSai_kNN,DiemSai_kNMAP,TongSoDiem,' ...
        'PL_LOS_mean,PL_LOS_std,PL_LOS_min,PL_LOS_max,' ...
        'PL_NLOS_mean,PL_NLOS_std,PL_NLOS_min,PL_NLOS_max,' ...
        'SNR_LOS_dB,SNR_NLOS_dB,' ...
        'Spd_LOS_mean,Spd_LOS_std,Spd_LOS_min,Spd_LOS_max,' ...
        'Spd_NLOS_mean,Spd_NLOS_std,Spd_NLOS_min,Spd_NLOS_max,' ...
        'SpeedDiff_Mbps,SpeedDiff_pct,' ...
        'MAE_PL_kNN,MAE_PL_kNMAP,MAE_Spd_kNN,MAE_Spd_kNMAP,' ...
        'AccGain_kNMAP_minus_kNN_pct,PL_MAE_Gain_kNN_minus_kNMAP,Spd_MAE_Gain_kNN_minus_kNMAP,Ghi_chu\n']);
end

fprintf(log_file, '==========================================\n');
fprintf(log_file, 'Ngay chay         : %s\n', ngay_chay);
fprintf(log_file, 'Script            : main_simulation.m (patched)\n');
fprintf(log_file, '\n--- Tham so tham chieu ---\n');
fprintf(log_file, 'lambda_ref/Lmax_ref/R : %.5f / %d m / %d m\n', lambda_ref, Lmax_ref, R);
fprintf(log_file, '\n--- Tham so thuat toan ---\n');
fprintf(log_file, 'N_data / K / N_mesh   : %d / %d / %d\n', N_data, K, N_mesh);
fprintf(log_file, '\n--- Tham so 3GPP TR 38.901 ---\n');
fprintf(log_file, 'freq / B / P_tx   : %.1f GHz / %.0f MHz / %d dBm\n', freq_GHz, B/1e6, P_tx_dBm);
fprintf(log_file, 'G_ant / NF        : %d dBi / %d dB\n', G_ant_dBi, NF_dB);
fprintf(log_file, 'h_BS/h_UT/h_E     : %d/%.1f/%.1f m\n', h_BS, h_UT, h_E);
fprintf(log_file, 'Noise floor       : %.2f dBm\n', Noise_dBm);
fprintf(log_file, 'd_BP              : %.1f m\n', d_BP);
fprintf(log_file, 'Ghi chu           : Downstream duoc tinh tu ca ground truth va nhan du doan.\n');

%% ================== SINGLE RUN ==================
[lambda_scene, Lmax_scene, buildings] = chicago_buildings(R);
L            = 2*R;
total_shadow = generate_shadows(R, buildings);
data         = generate_data(N_data, total_shadow, R);

figure('Color','w');
plot_estimation(L, R, total_shadow, data, [], [], 'Data and shadows', buildings);
export_current_figure(out_dir, '01_data_and_shadows');

[x, y]  = meshgrid(linspace(-R,R,N_mesh), linspace(-R,R,N_mesh));
in_disk = x.^2 + y.^2 <= R^2;
mesh    = [x(in_disk), y(in_disk)];
total_pts = size(mesh, 1);

estimated_kNN   = kNN_estimator(data, mesh, K);
estimated_kNMAP = kNMAP_estimator(data, lambda_scene, Lmax_scene, mesh, K);
true_labels     = compute_true_labels_from_shadow(total_shadow, mesh);

precision_kNN   = estimation_precision(mesh, estimated_kNN,   total_shadow);
precision_kNMAP = estimation_precision(mesh, estimated_kNMAP, total_shadow);
err_kNN   = sum(estimated_kNN(:)   ~= true_labels(:));
err_kNMAP = sum(estimated_kNMAP(:) ~= true_labels(:));

% Ground-truth downstream
[PL_true, ~, SNR_true_lin, ~, Speed_true, dist_sr] = compute_downstream_from_labels(true_labels, mesh, ...
    PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);

% Predicted downstream
[PL_kNN, ~, ~, ~, Speed_kNN_pred] = compute_downstream_from_labels(estimated_kNN, mesh, ...
    PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);
[PL_kNMAP, ~, ~, ~, Speed_kNMAP_pred] = compute_downstream_from_labels(estimated_kNMAP, mesh, ...
    PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);

% Truth-based summaries
sr_PL_LOS   = PL_true(true_labels==1);
sr_PL_NLOS  = PL_true(true_labels==0);
sr_Spd_LOS  = Speed_true(true_labels==1);
sr_Spd_NLOS = Speed_true(true_labels==0);
sr_SNR_LOS_dB  = 10*log10(mean(SNR_true_lin(true_labels==1)));
sr_SNR_NLOS_dB = 10*log10(mean(SNR_true_lin(true_labels==0)));
sr_SpeedDiff     = mean(sr_Spd_LOS) - mean(sr_Spd_NLOS);
sr_SpeedDiff_pct = sr_SpeedDiff / mean(sr_Spd_LOS) * 100;

% Predicted downstream errors against ground truth
sr_mae_PL_kNN     = mean(abs(PL_kNN - PL_true));
sr_mae_PL_kNMAP   = mean(abs(PL_kNMAP - PL_true));
sr_mae_Spd_kNN    = mean(abs(Speed_kNN_pred - Speed_true));
sr_mae_Spd_kNMAP  = mean(abs(Speed_kNMAP_pred - Speed_true));

fprintf('\n--- SINGLE RUN ---\n');
fprintf('Scene lambda/Lmax: %.5f / %.2f m\n', lambda_scene, Lmax_scene);
fprintf('Tong diem: %d | kNN: %.2f%% (%d sai) | kNMAP: %.2f%% (%d sai)\n', ...
    total_pts, precision_kNN*100, err_kNN, precision_kNMAP*100, err_kNMAP);
fprintf('PL  LOS : %.2f+/-%.2f [%.2f,%.2f] dB\n', mean(sr_PL_LOS), std(sr_PL_LOS), min(sr_PL_LOS), max(sr_PL_LOS));
fprintf('PL  NLOS: %.2f+/-%.2f [%.2f,%.2f] dB\n', mean(sr_PL_NLOS),std(sr_PL_NLOS),min(sr_PL_NLOS),max(sr_PL_NLOS));
fprintf('SNR LOS : %.2f dB | SNR NLOS: %.2f dB\n', sr_SNR_LOS_dB, sr_SNR_NLOS_dB);
fprintf('Spd LOS : %.2f+/-%.2f [%.2f,%.2f] Mbps\n', mean(sr_Spd_LOS), std(sr_Spd_LOS), min(sr_Spd_LOS), max(sr_Spd_LOS));
fprintf('Spd NLOS: %.2f+/-%.2f [%.2f,%.2f] Mbps\n', mean(sr_Spd_NLOS),std(sr_Spd_NLOS),min(sr_Spd_NLOS),max(sr_Spd_NLOS));
fprintf('Spd diff: %.2f Mbps (%.1f%%)\n', sr_SpeedDiff, sr_SpeedDiff_pct);
fprintf('Downstream MAE | PL: kNN %.2f dB, kNMAP %.2f dB | Speed: kNN %.2f Mbps, kNMAP %.2f Mbps\n', ...
    sr_mae_PL_kNN, sr_mae_PL_kNMAP, sr_mae_Spd_kNN, sr_mae_Spd_kNMAP);

fprintf(log_file, '\n--- Single run ---\n');
fprintf(log_file, 'Scene lambda/Lmax       : %.5f / %.2f m\n', lambda_scene, Lmax_scene);
fprintf(log_file, 'Tong diem luoi          : %d\n', total_pts);
fprintf(log_file, 'kNN   acc | diem sai    : %.2f%% | %d/%d\n', precision_kNN*100,   err_kNN,   total_pts);
fprintf(log_file, 'kNMAP acc | diem sai    : %.2f%% | %d/%d\n', precision_kNMAP*100, err_kNMAP, total_pts);
fprintf(log_file, 'PL LOS  mean/std/min/max: %.2f/%.2f/%.2f/%.2f dB\n', mean(sr_PL_LOS), std(sr_PL_LOS), min(sr_PL_LOS), max(sr_PL_LOS));
fprintf(log_file, 'PL NLOS mean/std/min/max: %.2f/%.2f/%.2f/%.2f dB\n', mean(sr_PL_NLOS),std(sr_PL_NLOS),min(sr_PL_NLOS),max(sr_PL_NLOS));
fprintf(log_file, 'SNR LOS  mean           : %.2f dB\n', sr_SNR_LOS_dB);
fprintf(log_file, 'SNR NLOS mean           : %.2f dB\n', sr_SNR_NLOS_dB);
fprintf(log_file, 'Speed LOS  mean/std/min/max : %.2f/%.2f/%.2f/%.2f Mbps\n', mean(sr_Spd_LOS), std(sr_Spd_LOS), min(sr_Spd_LOS), max(sr_Spd_LOS));
fprintf(log_file, 'Speed NLOS mean/std/min/max : %.2f/%.2f/%.2f/%.2f Mbps\n', mean(sr_Spd_NLOS),std(sr_Spd_NLOS),min(sr_Spd_NLOS),max(sr_Spd_NLOS));
fprintf(log_file, 'Speed diff LOS-NLOS     : %.2f Mbps (%.1f%%)\n', sr_SpeedDiff, sr_SpeedDiff_pct);
fprintf(log_file, 'Downstream MAE PL       : kNN %.2f dB | kNMAP %.2f dB\n', sr_mae_PL_kNN, sr_mae_PL_kNMAP);
fprintf(log_file, 'Downstream MAE Speed    : kNN %.2f Mbps | kNMAP %.2f Mbps\n', sr_mae_Spd_kNN, sr_mae_Spd_kNMAP);

fprintf(csv_file, '%s,Single_run,%d,%d,%.5f,%.2f,%d,%.1f,%.0f,%d,%d,%.2f,%.1f,1,', ...
    ngay_chay,N_data,K,lambda_scene,Lmax_scene,R,freq_GHz,B/1e6,P_tx_dBm,NF_dB,Noise_dBm,d_BP);
fprintf(csv_file, '%.2f,0,%.2f,0,%d,%d,%d,', precision_kNN*100,precision_kNMAP*100,err_kNN,err_kNMAP,total_pts);
fprintf(csv_file, '%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,', ...
    mean(sr_PL_LOS),std(sr_PL_LOS),min(sr_PL_LOS),max(sr_PL_LOS), ...
    mean(sr_PL_NLOS),std(sr_PL_NLOS),min(sr_PL_NLOS),max(sr_PL_NLOS));
fprintf(csv_file, '%.2f,%.2f,', sr_SNR_LOS_dB, sr_SNR_NLOS_dB);
fprintf(csv_file, '%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,', ...
    mean(sr_Spd_LOS),std(sr_Spd_LOS),min(sr_Spd_LOS),max(sr_Spd_LOS), ...
    mean(sr_Spd_NLOS),std(sr_Spd_NLOS),min(sr_Spd_NLOS),max(sr_Spd_NLOS));
fprintf(csv_file, '%.2f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,Single run patched\n', ...
    sr_SpeedDiff, sr_SpeedDiff_pct, ...
    sr_mae_PL_kNN, sr_mae_PL_kNMAP, sr_mae_Spd_kNN, sr_mae_Spd_kNMAP, ...
    (precision_kNMAP-precision_kNN)*100, sr_mae_PL_kNN-sr_mae_PL_kNMAP, sr_mae_Spd_kNN-sr_mae_Spd_kNMAP);

figure('Color','w');
subplot(1,2,1);
plot_estimation(L, R, total_shadow, data, mesh, estimated_kNN,   'kNN',   buildings);
subplot(1,2,2);
plot_estimation(L, R, total_shadow, data, mesh, estimated_kNMAP, 'kNMAP', buildings);
export_current_figure(out_dir, '02_prediction_knn_vs_knmap');

figure('Color','w');
plot_estimation(L, R, total_shadow, data, mesh, estimated_kNMAP, 'Error kNMAP', buildings);
hold on;
err = estimated_kNMAP(:) ~= true_labels(:);
scatter(mesh(err,1), mesh(err,2), 40, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Errors');
hold off;
export_current_figure(out_dir, '03_error_knmap');

figure('Color','w');
plot_estimation(L, R, total_shadow, data, mesh, estimated_kNN, 'Error kNN', buildings);
hold on;
err2 = estimated_kNN(:) ~= true_labels(:);
scatter(mesh(err2,1), mesh(err2,2), 40, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Errors');
hold off;
export_current_figure(out_dir, '04_error_knn');

figure('Color','w'); hold on;
scatter(mesh(:,1), mesh(:,2), 55, Speed_true, 'filled', 'MarkerEdgeColor','k');
plot_estimation(L, R, total_shadow, data, [], [], 'True speed map', buildings);
title('Heatmap toc do (ground truth LOS/NLOS) - 3GPP TR 38.901');
axis equal; xlim([-R R]); ylim([-R R]);
cb = colorbar; cb.Label.String = 'Toc do (Mbps)'; colormap(jet); hold off;
export_current_figure(out_dir, '05_speed_heatmap_true');

figure('Color','w');
subplot(1,2,1); hold on;
scatter(mesh(:,1), mesh(:,2), 55, Speed_kNN_pred, 'filled', 'MarkerEdgeColor','k');
plot_estimation(L, R, total_shadow, data, [], [], 'Predicted speed map: kNN', buildings);
axis equal; xlim([-R R]); ylim([-R R]); cb1 = colorbar; cb1.Label.String = 'Toc do (Mbps)'; colormap(jet); hold off;
subplot(1,2,2); hold on;
scatter(mesh(:,1), mesh(:,2), 55, Speed_kNMAP_pred, 'filled', 'MarkerEdgeColor','k');
plot_estimation(L, R, total_shadow, data, [], [], 'Predicted speed map: kNMAP', buildings);
axis equal; xlim([-R R]); ylim([-R R]); cb2 = colorbar; cb2.Label.String = 'Toc do (Mbps)'; colormap(jet); hold off;
export_current_figure(out_dir, '06_speed_heatmap_predicted_knn_vs_knmap');

figure('Color','w');
b = bar([sr_mae_PL_kNN, sr_mae_PL_kNMAP; sr_mae_Spd_kNN, sr_mae_Spd_kNMAP]);
set(gca, 'XTickLabel', {'MAE Path Loss (dB)','MAE Speed (Mbps)'}, 'FontSize', 11);
legend({'kNN','kNMAP'}, 'Location', 'best');
grid on; box on;
title('Single-run downstream error from predicted labels', 'Interpreter', 'none');
export_current_figure(out_dir, '07_downstream_mae_single_run');

%% ================== MONTE CARLO ==================
N_run = 100;
fprintf(log_file, '\n--- Monte Carlo (N_run = %d) ---\n', N_run);

results_kNN = zeros(N_run,1); results_kNMAP = zeros(N_run,1);
Speed_LOS_all=[]; Speed_NLOS_all=[]; PL_LOS_all=[]; PL_NLOS_all=[];
d_LOS_all=[]; d_NLOS_all=[];
mean_LOS_per_run=zeros(N_run,1); mean_NLOS_per_run=zeros(N_run,1);
lambda_scene_all = zeros(N_run,1); Lmax_scene_all = zeros(N_run,1);
mae_PL_kNN_runs = zeros(N_run,1); mae_PL_kNMAP_runs = zeros(N_run,1);
mae_Spd_kNN_runs = zeros(N_run,1); mae_Spd_kNMAP_runs = zeros(N_run,1);

for i = 1:N_run
    [lambda_scene_i, Lmax_scene_i, buildings_i] = chicago_buildings(R);
    lambda_scene_all(i) = lambda_scene_i;
    Lmax_scene_all(i)   = Lmax_scene_i;

    total_shadow_i = generate_shadows(R, buildings_i);
    data_i         = generate_data(N_data, total_shadow_i, R);
    true_labels_i  = compute_true_labels_from_shadow(total_shadow_i, mesh);

    est_kNN_i   = kNN_estimator(data_i, mesh, K);
    est_kNMAP_i = kNMAP_estimator(data_i, lambda_scene_i, Lmax_scene_i, mesh, K);
    results_kNN(i)   = estimation_precision(mesh, est_kNN_i,   total_shadow_i);
    results_kNMAP(i) = estimation_precision(mesh, est_kNMAP_i, total_shadow_i);

    [PL_true_i, ~, ~, ~, Speed_true_i, dist_i] = compute_downstream_from_labels(true_labels_i, mesh, ...
        PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);
    [PL_kNN_i, ~, ~, ~, Speed_kNN_i] = compute_downstream_from_labels(est_kNN_i, mesh, ...
        PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);
    [PL_kNMAP_i, ~, ~, ~, Speed_kNMAP_i] = compute_downstream_from_labels(est_kNMAP_i, mesh, ...
        PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);

    mae_PL_kNN_runs(i)   = mean(abs(PL_kNN_i - PL_true_i));
    mae_PL_kNMAP_runs(i) = mean(abs(PL_kNMAP_i - PL_true_i));
    mae_Spd_kNN_runs(i)   = mean(abs(Speed_kNN_i - Speed_true_i));
    mae_Spd_kNMAP_runs(i) = mean(abs(Speed_kNMAP_i - Speed_true_i));

    Speed_LOS_all  = [Speed_LOS_all;  Speed_true_i(true_labels_i==1)];
    Speed_NLOS_all = [Speed_NLOS_all; Speed_true_i(true_labels_i==0)];
    mean_LOS_per_run(i)  = mean(Speed_true_i(true_labels_i==1));
    mean_NLOS_per_run(i) = mean(Speed_true_i(true_labels_i==0));
    PL_LOS_all  = [PL_LOS_all;  PL_true_i(true_labels_i==1)];
    PL_NLOS_all = [PL_NLOS_all; PL_true_i(true_labels_i==0)];
    d_LOS_all   = [d_LOS_all;   dist_i(true_labels_i==1)];
    d_NLOS_all  = [d_NLOS_all;  dist_i(true_labels_i==0)];
end

mean_kNN=mean(results_kNN); std_kNN=std(results_kNN);
mean_kNMAP=mean(results_kNMAP); std_kNMAP=std(results_kNMAP);
mean_LOS_total=mean(mean_LOS_per_run); std_LOS_total=std(mean_LOS_per_run);
mean_NLOS_total=mean(mean_NLOS_per_run); std_NLOS_total=std(mean_NLOS_per_run);
mc_PL_LOS_mean=mean(PL_LOS_all); mc_PL_LOS_std=std(PL_LOS_all);
mc_PL_LOS_min=min(PL_LOS_all);   mc_PL_LOS_max=max(PL_LOS_all);
mc_PL_NLOS_mean=mean(PL_NLOS_all); mc_PL_NLOS_std=std(PL_NLOS_all);
mc_PL_NLOS_min=min(PL_NLOS_all);   mc_PL_NLOS_max=max(PL_NLOS_all);
mc_Spd_LOS_min=min(Speed_LOS_all); mc_Spd_LOS_max=max(Speed_LOS_all);
mc_Spd_NLOS_min=min(Speed_NLOS_all); mc_Spd_NLOS_max=max(Speed_NLOS_all);
mc_SpeedDiff = mean_LOS_total - mean_NLOS_total;
mc_SpeedDiff_pct = mc_SpeedDiff / mean_LOS_total * 100;
mc_SNR_LOS_dB  = (P_tx_dBm+G_ant_dBi) - mc_PL_LOS_mean  - Noise_dBm;
mc_SNR_NLOS_dB = (P_tx_dBm+G_ant_dBi) - mc_PL_NLOS_mean - Noise_dBm;
mc_lambda_mean = mean(lambda_scene_all); mc_lambda_std = std(lambda_scene_all);
mc_Lmax_mean   = mean(Lmax_scene_all);   mc_Lmax_std   = std(Lmax_scene_all);
mc_mae_PL_kNN = mean(mae_PL_kNN_runs); mc_mae_PL_kNMAP = mean(mae_PL_kNMAP_runs);
mc_mae_Spd_kNN = mean(mae_Spd_kNN_runs); mc_mae_Spd_kNMAP = mean(mae_Spd_kNMAP_runs);

fprintf('\n--- KET QUA MONTE CARLO ---\n');
fprintf('Scene lambda mean/std : %.5f +/- %.5f\n', mc_lambda_mean, mc_lambda_std);
fprintf('Scene Lmax   mean/std : %.2f +/- %.2f m\n', mc_Lmax_mean, mc_Lmax_std);
fprintf('kNN  : %.2f%%+/-%.2f%% | kNMAP: %.2f%%+/-%.2f%% | kNMAP hon: %.2f%%\n', ...
    mean_kNN*100,std_kNN*100,mean_kNMAP*100,std_kNMAP*100,(mean_kNMAP-mean_kNN)*100);
fprintf('PL LOS  : %.2f+/-%.2f [%.2f,%.2f] dB\n', mc_PL_LOS_mean,mc_PL_LOS_std,mc_PL_LOS_min,mc_PL_LOS_max);
fprintf('PL NLOS : %.2f+/-%.2f [%.2f,%.2f] dB\n', mc_PL_NLOS_mean,mc_PL_NLOS_std,mc_PL_NLOS_min,mc_PL_NLOS_max);
fprintf('SNR LOS : %.2f dB | SNR NLOS: %.2f dB\n', mc_SNR_LOS_dB, mc_SNR_NLOS_dB);
fprintf('Spd LOS : %.2f+/-%.2f [%.2f,%.2f] Mbps\n', mean_LOS_total,std_LOS_total,mc_Spd_LOS_min,mc_Spd_LOS_max);
fprintf('Spd NLOS: %.2f+/-%.2f [%.2f,%.2f] Mbps\n', mean_NLOS_total,std_NLOS_total,mc_Spd_NLOS_min,mc_Spd_NLOS_max);
fprintf('Spd diff: %.2f Mbps (%.1f%%)\n', mc_SpeedDiff, mc_SpeedDiff_pct);
fprintf('Downstream MAE | PL: kNN %.2f dB, kNMAP %.2f dB | Speed: kNN %.2f Mbps, kNMAP %.2f Mbps\n', ...
    mc_mae_PL_kNN, mc_mae_PL_kNMAP, mc_mae_Spd_kNN, mc_mae_Spd_kNMAP);

fprintf(log_file, 'Scene lambda mean/std   : %.5f +/- %.5f\n', mc_lambda_mean, mc_lambda_std);
fprintf(log_file, 'Scene Lmax   mean/std   : %.2f +/- %.2f m\n', mc_Lmax_mean, mc_Lmax_std);
fprintf(log_file, 'kNN   acc mean/std      : %.2f%% +/- %.2f%%\n', mean_kNN*100,   std_kNN*100);
fprintf(log_file, 'kNMAP acc mean/std      : %.2f%% +/- %.2f%%\n', mean_kNMAP*100, std_kNMAP*100);
fprintf(log_file, 'kNMAP cao hon kNN       : %.2f%%\n', (mean_kNMAP-mean_kNN)*100);
fprintf(log_file, 'PL LOS  mean/std/min/max: %.2f/%.2f/%.2f/%.2f dB\n', mc_PL_LOS_mean,mc_PL_LOS_std,mc_PL_LOS_min,mc_PL_LOS_max);
fprintf(log_file, 'PL NLOS mean/std/min/max: %.2f/%.2f/%.2f/%.2f dB\n', mc_PL_NLOS_mean,mc_PL_NLOS_std,mc_PL_NLOS_min,mc_PL_NLOS_max);
fprintf(log_file, 'PL NLOS cao hon LOS     : %.2f dB\n', mc_PL_NLOS_mean-mc_PL_LOS_mean);
fprintf(log_file, 'SNR LOS  mean           : %.2f dB\n', mc_SNR_LOS_dB);
fprintf(log_file, 'SNR NLOS mean           : %.2f dB\n', mc_SNR_NLOS_dB);
fprintf(log_file, 'SNR LOS cao hon NLOS    : %.2f dB\n', mc_SNR_LOS_dB-mc_SNR_NLOS_dB);
fprintf(log_file, 'Speed LOS  mean/std/min/max : %.2f/%.2f/%.2f/%.2f Mbps\n', mean_LOS_total,std_LOS_total,mc_Spd_LOS_min,mc_Spd_LOS_max);
fprintf(log_file, 'Speed NLOS mean/std/min/max : %.2f/%.2f/%.2f/%.2f Mbps\n', mean_NLOS_total,std_NLOS_total,mc_Spd_NLOS_min,mc_Spd_NLOS_max);
fprintf(log_file, 'Speed diff LOS-NLOS     : %.2f Mbps (%.1f%%)\n', mc_SpeedDiff, mc_SpeedDiff_pct);
fprintf(log_file, 'Downstream MAE PL       : kNN %.2f dB | kNMAP %.2f dB\n', mc_mae_PL_kNN, mc_mae_PL_kNMAP);
fprintf(log_file, 'Downstream MAE Speed    : kNN %.2f Mbps | kNMAP %.2f Mbps\n', mc_mae_Spd_kNN, mc_mae_Spd_kNMAP);

fprintf(csv_file, '%s,Monte_Carlo,%d,%d,%.5f,%.2f,%d,%.1f,%.0f,%d,%d,%.2f,%.1f,%d,', ...
    ngay_chay,N_data,K,mc_lambda_mean,mc_Lmax_mean,R,freq_GHz,B/1e6,P_tx_dBm,NF_dB,Noise_dBm,d_BP,N_run);
fprintf(csv_file, '%.2f,%.2f,%.2f,%.2f,0,0,%d,', mean_kNN*100,std_kNN*100,mean_kNMAP*100,std_kNMAP*100,total_pts);
fprintf(csv_file, '%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,', ...
    mc_PL_LOS_mean,mc_PL_LOS_std,mc_PL_LOS_min,mc_PL_LOS_max, ...
    mc_PL_NLOS_mean,mc_PL_NLOS_std,mc_PL_NLOS_min,mc_PL_NLOS_max);
fprintf(csv_file, '%.2f,%.2f,', mc_SNR_LOS_dB, mc_SNR_NLOS_dB);
fprintf(csv_file, '%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,', ...
    mean_LOS_total,std_LOS_total,mc_Spd_LOS_min,mc_Spd_LOS_max, ...
    mean_NLOS_total,std_NLOS_total,mc_Spd_NLOS_min,mc_Spd_NLOS_max);
fprintf(csv_file, '%.2f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,Monte Carlo baseline patched\n', ...
    mc_SpeedDiff,mc_SpeedDiff_pct, ...
    mc_mae_PL_kNN, mc_mae_PL_kNMAP, mc_mae_Spd_kNN, mc_mae_Spd_kNMAP, ...
    (mean_kNMAP-mean_kNN)*100, mc_mae_PL_kNN-mc_mae_PL_kNMAP, mc_mae_Spd_kNN-mc_mae_Spd_kNMAP);

%% Biểu đồ Monte Carlo
nbins_s=15; all_speeds=[Speed_LOS_all;Speed_NLOS_all];
edges_S=linspace(min(all_speeds),max(all_speeds),nbins_s+1);
[cnt_S_LOS,~]=histcounts(Speed_LOS_all,edges_S);
[cnt_S_NLOS,~]=histcounts(Speed_NLOS_all,edges_S);
bin_S=(edges_S(1:end-1)+edges_S(2:end))/2;
figure('Color','w');
b=bar(bin_S,[cnt_S_LOS(:),cnt_S_NLOS(:)],'grouped','BarWidth',1.0);
b(1).FaceColor=[0 0.45 0.74]; b(2).FaceColor=[0.85 0.33 0.1];
b(1).FaceAlpha=0.85; b(2).FaceAlpha=0.85;
xlabel('Toc do (Mbps)','FontSize',12); ylabel('So diem','FontSize',12);
title('Phan bo toc do LOS vs NLOS (Monte Carlo)','FontSize',14,'FontWeight','bold');
legend({'LOS','NLOS'},'Location','best','FontSize',11); grid on; box on;
export_current_figure(out_dir, '08_speed_histogram_mc');

figure('Color','w'); hold on;
plot(1:N_run,mean_LOS_per_run,'b-o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor',[0 0.45 0.74],'DisplayName','LOS');
plot(1:N_run,mean_NLOS_per_run,'r-s','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor',[0.85 0.33 0.1],'DisplayName','NLOS');
yline(mean_LOS_total,'b--','LineWidth',1.5,'DisplayName',sprintf('TB LOS = %.1f Mbps',mean_LOS_total));
yline(mean_NLOS_total,'r--','LineWidth',1.5,'DisplayName',sprintf('TB NLOS = %.1f Mbps',mean_NLOS_total));
xlabel('Lan chay','FontSize',12); ylabel('Toc do TB (Mbps)','FontSize',12);
title('Toc do TB LOS vs NLOS theo tung lan Monte Carlo','FontSize',14,'FontWeight','bold');
legend('Location','best','FontSize',10); xlim([1 N_run]); grid on; box on; hold off;
export_current_figure(out_dir, '09_speed_mean_per_run_mc');

figure('Color','w');
b3=bar([mean_LOS_total,mean_NLOS_total],'BarWidth',0.5,'FaceColor','flat');
b3.CData(1,:)=[0 0.45 0.74]; b3.CData(2,:)=[0.85 0.33 0.1]; b3.FaceAlpha=0.85;
set(gca,'XTickLabel',{'LOS','NLOS'},'FontSize',12);
ylabel('Toc do TB (Mbps)','FontSize',12);
title('Toc do TB LOS vs NLOS (Tong hop Monte Carlo)','FontSize',14,'FontWeight','bold');
grid on; box on; ylim([0,max(mean_LOS_total,mean_NLOS_total)*1.25]);
text(1,mean_LOS_total+0.5,sprintf('%.2f Mbps',mean_LOS_total),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
text(2,mean_NLOS_total+0.5,sprintf('%.2f Mbps',mean_NLOS_total),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
export_current_figure(out_dir, '10_speed_mean_bar_mc');

figure('Color','w');
b2=bar([mean_kNN*100 mean_kNMAP*100],'BarWidth',0.6,'FaceColor','flat');
b2.CData(1,:)=[0.3 0.75 0.9]; b2.CData(2,:)=[0.85 0.33 0.1];
set(gca,'XTickLabel',{'kNN','kNMAP'},'FontSize',12);
ylabel('Precision (%)','FontSize',12); title('Monte Carlo Average Precision','FontSize',14,'FontWeight','bold');
grid on; box on;
text(1,mean_kNN*100,sprintf('%.2f%%',mean_kNN*100),'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
text(2,mean_kNMAP*100,sprintf('%.2f%%',mean_kNMAP*100),'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
export_current_figure(out_dir, '11_precision_bar_mc');

figure('Color','w');
b4 = bar([mc_mae_PL_kNN, mc_mae_PL_kNMAP; mc_mae_Spd_kNN, mc_mae_Spd_kNMAP]);
set(gca, 'XTickLabel', {'MAE Path Loss (dB)','MAE Speed (Mbps)'}, 'FontSize', 11);
legend({'kNN','kNMAP'}, 'Location', 'best');
grid on; box on;
title('Monte Carlo downstream error from predicted labels', 'Interpreter', 'none');
export_current_figure(out_dir, '12_downstream_mae_mc');

d_range=linspace(1,R,500);
PL_LOS_theory=PL_LOS_func(d_range); PL_NLOS_theory=PL_NLOS_func(d_range);
figure('Color','w'); hold on;
scatter(d_LOS_all,PL_LOS_all,8,[0 0.45 0.74],'filled','MarkerFaceAlpha',0.12,'DisplayName','LOS (data)');
scatter(d_NLOS_all,PL_NLOS_all,8,[0.85 0.33 0.1],'filled','MarkerFaceAlpha',0.12,'DisplayName','NLOS (data)');
plot(d_range,PL_LOS_theory,'b-','LineWidth',2.5,'DisplayName','LOS (3GPP)');
plot(d_range,PL_NLOS_theory,'r-','LineWidth',2.5,'DisplayName','NLOS (3GPP)');
xline(d_BP,'k--','LineWidth',1.5,'DisplayName',sprintf('d_{BP}=%.0f m',d_BP));
xlabel('Khoang cach (m)','FontSize',12); ylabel('Path Loss (dB)','FontSize',12);
title('Suy hao LOS vs NLOS - 3GPP TR 38.901','FontSize',14,'FontWeight','bold');
legend('Location','northwest','FontSize',10); grid on; box on; xlim([0 R]); hold off;
export_current_figure(out_dir, '13_pathloss_vs_distance_mc');

nbins_pl=15;
[cnt_LOS,edges_LOS]=histcounts(PL_LOS_all,nbins_pl);
[cnt_NLOS,edges_NLOS]=histcounts(PL_NLOS_all,nbins_pl);
bin_LOS=(edges_LOS(1:end-1)+edges_LOS(2:end))/2;
bin_NLOS=(edges_NLOS(1:end-1)+edges_NLOS(2:end))/2;
figure('Color','w'); hold on;
barh(bin_LOS,cnt_LOS,1.0,'FaceAlpha',0.7,'FaceColor',[0 0.45 0.74],'DisplayName','LOS');
barh(bin_NLOS,-cnt_NLOS,1.0,'FaceAlpha',0.7,'FaceColor',[0.85 0.33 0.1],'DisplayName','NLOS');
ylabel('Path Loss (dB)','FontSize',12); xlabel('So diem','FontSize',12);
title('Phan bo Path Loss LOS vs NLOS - 3GPP TR 38.901','FontSize',14,'FontWeight','bold');
legend('Location','best','FontSize',10);
xt2=get(gca,'XTick'); set(gca,'XTickLabel',abs(xt2));
grid on; box on; hold off;
export_current_figure(out_dir, '14_pathloss_histogram_mc');

%% ================== SWEEP N_DATA ==================
fprintf('\n--- SWEEP N_DATA ---\n');
fprintf(log_file, '\n--- Sweep N_data (N_run_sweep = 100) ---\n');

Ndata_list=[10,20,40,80,120]; N_run_sweep=100;
acc_kNN_sweep=zeros(length(Ndata_list),1); acc_kNMAP_sweep=zeros(length(Ndata_list),1);
std_kNN_sweep=zeros(length(Ndata_list),1); std_kNMAP_sweep=zeros(length(Ndata_list),1);
mae_pl_knn_sweep=zeros(length(Ndata_list),1); mae_pl_knmap_sweep=zeros(length(Ndata_list),1);
mae_spd_knn_sweep=zeros(length(Ndata_list),1); mae_spd_knmap_sweep=zeros(length(Ndata_list),1);
lambda_sweep_mean=zeros(length(Ndata_list),1); Lmax_sweep_mean=zeros(length(Ndata_list),1);

for idx = 1:length(Ndata_list)
    N_data_sw=Ndata_list(idx);
    acc_kNN_buf=zeros(N_run_sweep,1); acc_kNMAP_buf=zeros(N_run_sweep,1);
    mae_pl_knn_buf=zeros(N_run_sweep,1); mae_pl_knmap_buf=zeros(N_run_sweep,1);
    mae_spd_knn_buf=zeros(N_run_sweep,1); mae_spd_knmap_buf=zeros(N_run_sweep,1);
    lambda_buf=zeros(N_run_sweep,1); Lmax_buf=zeros(N_run_sweep,1);
    for j = 1:N_run_sweep
        [lambda_sw,Lmax_sw,buildings_sw]=chicago_buildings(R);
        lambda_buf(j)=lambda_sw; Lmax_buf(j)=Lmax_sw;
        total_shadow_sw=generate_shadows(R,buildings_sw);
        data_sw=generate_data(N_data_sw,total_shadow_sw,R);
        true_labels_sw = compute_true_labels_from_shadow(total_shadow_sw, mesh);
        est_kNN_sw  =kNN_estimator(data_sw,mesh,K);
        est_kNMAP_sw=kNMAP_estimator(data_sw,lambda_sw,Lmax_sw,mesh,K);
        acc_kNN_buf(j)  =estimation_precision(mesh,est_kNN_sw,  total_shadow_sw);
        acc_kNMAP_buf(j)=estimation_precision(mesh,est_kNMAP_sw,total_shadow_sw);
        [PL_true_sw, ~, ~, ~, Speed_true_sw] = compute_downstream_from_labels(true_labels_sw, mesh, ...
            PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);
        [PL_kNN_sw, ~, ~, ~, Speed_kNN_sw] = compute_downstream_from_labels(est_kNN_sw, mesh, ...
            PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);
        [PL_kNMAP_sw, ~, ~, ~, Speed_kNMAP_sw] = compute_downstream_from_labels(est_kNMAP_sw, mesh, ...
            PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B);
        mae_pl_knn_buf(j)   = mean(abs(PL_kNN_sw - PL_true_sw));
        mae_pl_knmap_buf(j) = mean(abs(PL_kNMAP_sw - PL_true_sw));
        mae_spd_knn_buf(j)   = mean(abs(Speed_kNN_sw - Speed_true_sw));
        mae_spd_knmap_buf(j) = mean(abs(Speed_kNMAP_sw - Speed_true_sw));
    end
    acc_kNN_sweep(idx)  =mean(acc_kNN_buf)  *100; std_kNN_sweep(idx)  =std(acc_kNN_buf)  *100;
    acc_kNMAP_sweep(idx)=mean(acc_kNMAP_buf)*100; std_kNMAP_sweep(idx)=std(acc_kNMAP_buf)*100;
    mae_pl_knn_sweep(idx) = mean(mae_pl_knn_buf); mae_pl_knmap_sweep(idx) = mean(mae_pl_knmap_buf);
    mae_spd_knn_sweep(idx) = mean(mae_spd_knn_buf); mae_spd_knmap_sweep(idx) = mean(mae_spd_knmap_buf);
    lambda_sweep_mean(idx)=mean(lambda_buf); Lmax_sweep_mean(idx)=mean(Lmax_buf);

    fprintf('Ndata=%3d | kNN: %.2f%%(+/-%.2f) | kNMAP: %.2f%%(+/-%.2f) | MAE PL: %.2f/%.2f dB | MAE Spd: %.2f/%.2f Mbps\n', ...
        N_data_sw,acc_kNN_sweep(idx),std_kNN_sweep(idx),acc_kNMAP_sweep(idx),std_kNMAP_sweep(idx), ...
        mae_pl_knn_sweep(idx), mae_pl_knmap_sweep(idx), mae_spd_knn_sweep(idx), mae_spd_knmap_sweep(idx));
    fprintf(log_file,'Ndata=%3d | kNN: %.2f%%(+/-%.2f) | kNMAP: %.2f%%(+/-%.2f) | MAE PL: %.2f/%.2f dB | MAE Spd: %.2f/%.2f Mbps\n', ...
        N_data_sw,acc_kNN_sweep(idx),std_kNN_sweep(idx),acc_kNMAP_sweep(idx),std_kNMAP_sweep(idx), ...
        mae_pl_knn_sweep(idx), mae_pl_knmap_sweep(idx), mae_spd_knn_sweep(idx), mae_spd_knmap_sweep(idx));

    fprintf(csv_file,'%s,Sweep_Ndata=%d,%d,%d,%.5f,%.2f,%d,%.1f,%.0f,%d,%d,%.2f,%.1f,%d,', ...
        ngay_chay,N_data_sw,N_data_sw,K,lambda_sweep_mean(idx),Lmax_sweep_mean(idx),R,freq_GHz,B/1e6,P_tx_dBm,NF_dB,Noise_dBm,d_BP,N_run_sweep);
    fprintf(csv_file,'%.2f,%.2f,%.2f,%.2f,0,0,%d,', ...
        acc_kNN_sweep(idx),std_kNN_sweep(idx),acc_kNMAP_sweep(idx),std_kNMAP_sweep(idx),total_pts);
    fprintf(csv_file,'0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,Sweep N_data patched\n', ...
        mae_pl_knn_sweep(idx), mae_pl_knmap_sweep(idx), mae_spd_knn_sweep(idx), mae_spd_knmap_sweep(idx), ...
        acc_kNMAP_sweep(idx)-acc_kNN_sweep(idx), mae_pl_knn_sweep(idx)-mae_pl_knmap_sweep(idx), mae_spd_knn_sweep(idx)-mae_spd_knmap_sweep(idx));
end

figure('Color','w'); hold on;
plot(Ndata_list,acc_kNN_sweep,'b-o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0 0.45 0.74],'DisplayName','kNN');
plot(Ndata_list,acc_kNMAP_sweep,'r-s','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0.85 0.33 0.1],'DisplayName','kNMAP');
xlabel('So diem do N\_data','FontSize',12); ylabel('Accuracy (%)','FontSize',12);
title('Accuracy theo so diem do (Sweep N\_data)','FontSize',14,'FontWeight','bold');
legend('Location','southeast','FontSize',11); grid on; box on; hold off;
export_current_figure(out_dir, '15_accuracy_vs_ndata');

figure('Color','w'); hold on;
plot(Ndata_list,mae_pl_knn_sweep,'b-o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0 0.45 0.74],'DisplayName','kNN');
plot(Ndata_list,mae_pl_knmap_sweep,'r-s','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0.85 0.33 0.1],'DisplayName','kNMAP');
xlabel('So diem do N\_data','FontSize',12); ylabel('MAE Path Loss (dB)','FontSize',12);
title('MAE path loss theo so diem do (Sweep N\_data)','FontSize',14,'FontWeight','bold');
legend('Location','northeast','FontSize',11); grid on; box on; hold off;
export_current_figure(out_dir, '16_mae_pathloss_vs_ndata');

figure('Color','w'); hold on;
plot(Ndata_list,mae_spd_knn_sweep,'b-o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0 0.45 0.74],'DisplayName','kNN');
plot(Ndata_list,mae_spd_knmap_sweep,'r-s','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[0.85 0.33 0.1],'DisplayName','kNMAP');
xlabel('So diem do N\_data','FontSize',12); ylabel('MAE Speed (Mbps)','FontSize',12);
title('MAE toc do theo so diem do (Sweep N\_data)','FontSize',14,'FontWeight','bold');
legend('Location','northeast','FontSize',11); grid on; box on; hold off;
export_current_figure(out_dir, '17_mae_speed_vs_ndata');

%% ================== TỔNG KẾT & ĐÓNG FILE ==================
fprintf(log_file, '\n--- Tong ket ---\n');
fprintf(log_file, 'kNMAP cao hon kNN         : %.2f%%\n',  (mean_kNMAP-mean_kNN)*100);
fprintf(log_file, 'Speed LOS/NLOS ratio      : %.2f\n',    mean_LOS_total/mean_NLOS_total);
fprintf(log_file, 'Sut giam speed NLOS       : %.1f%% so voi LOS\n', mc_SpeedDiff_pct);
fprintf(log_file, 'PL NLOS cao hon LOS       : %.2f dB\n', mc_PL_NLOS_mean-mc_PL_LOS_mean);
fprintf(log_file, 'SNR LOS cao hon NLOS      : %.2f dB\n', mc_SNR_LOS_dB-mc_SNR_NLOS_dB);
fprintf(log_file, 'Downstream MAE gain (PL)  : %.2f dB, nghiêng ve kNMAP neu duong\n', mc_mae_PL_kNN-mc_mae_PL_kNMAP);
fprintf(log_file, 'Downstream MAE gain (Spd) : %.2f Mbps, nghiêng ve kNMAP neu duong\n', mc_mae_Spd_kNN-mc_mae_Spd_kNMAP);
fprintf(log_file, 'Nguong bao hoa N_data     : ~80 diem (xem do thi)\n');
fprintf(log_file, 'Thu muc hinh sach         : outputs_clean/\n');
fprintf(log_file, '==========================================\n\n');

fclose(log_file);
fclose(csv_file);

fprintf('\n===================================\n');
fprintf('Log day du : log_results.txt\n');
fprintf('Bang CSV   : bang_ket_qua.csv\n');
fprintf('Hinh sach  : outputs_clean/\n');
fprintf('===================================\n');
