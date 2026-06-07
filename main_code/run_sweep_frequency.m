
clear; close all;
warning('off','all')

R = 150;
N_data = 80;
K = 3;
N_mesh = 12;
Lmax = 20;
N_run = 50;


lambda = 0.00088;

% Sweep frequency
freq_list_GHz = [2, 28, 39, 60];

% Constants
c = 3e8; k_B = 1.38e-23; T_K = 290;
B = 10e6; NF_dB = 9; P_tx_dBm = 30; G_ant_dBi = 5;

h_BS = 25; h_UT = 1.5; h_E = 1.0;

[x, y]  = meshgrid(linspace(-R,R,N_mesh), linspace(-R,R,N_mesh));
in_disk = x.^2 + y.^2 <= R^2;
mesh    = [x(in_disk), y(in_disk)];

% Result buffers
acc_kNN = zeros(size(freq_list_GHz));
acc_kNMAP = zeros(size(freq_list_GHz));
mae_pl_kNN = zeros(size(freq_list_GHz));
mae_pl_kNMAP = zeros(size(freq_list_GHz));
mae_spd_kNN = zeros(size(freq_list_GHz));
mae_spd_kNMAP = zeros(size(freq_list_GHz));

for idx = 1:length(freq_list_GHz)

    freq_GHz = freq_list_GHz(idx);
    freq = freq_GHz * 1e9;

    % Recompute noise
    Noise_dBm = 10*log10(k_B*T_K*B) + 30 + NF_dB;

    % Breakpoint distance
    d_BP = 4 * (h_BS-h_E) * (h_UT-h_E) * freq / c;

    % Path loss functions
    PL_LOS_func = @(d) (d < d_BP) .* ...
        (28.0 + 22*log10(max(d,1e-3)) + 20*log10(freq_GHz)) + ...
        (d >= d_BP).* ...
        (28.0 + 40*log10(max(d,1e-3)) + 20*log10(freq_GHz) ...
        - 9*log10(d_BP^2 + (h_BS - h_UT)^2));

    PL_NLOS_func = @(d) max(PL_LOS_func(d), ...
        13.54 + 39.08*log10(max(d,1e-3)) + 20*log10(freq_GHz) ...
        - 0.6*(h_UT - 1.5));

    % Buffers
    acc_buf_kNN = zeros(N_run,1); 
    acc_buf_kNMAP = zeros(N_run,1);

    mae_pl_buf_kNN = zeros(N_run,1); 
    mae_pl_buf_kNMAP = zeros(N_run,1);

    mae_spd_buf_kNN = zeros(N_run,1); 
    mae_spd_buf_kNMAP = zeros(N_run,1);

    for i = 1:N_run

        buildings = PPP_buildings(lambda, Lmax, R);
        total_shadow = generate_shadows(R, buildings);
        data = generate_data(N_data, total_shadow, R);

        true_labels = compute_true_labels_from_shadow(total_shadow, mesh);

        est_kNN = kNN_estimator(data, mesh, K);
        est_kNMAP = kNMAP_estimator(data, lambda, Lmax, mesh, K);

        acc_buf_kNN(i) = estimation_precision(mesh, est_kNN, total_shadow);
        acc_buf_kNMAP(i) = estimation_precision(mesh, est_kNMAP, total_shadow);

        [PL_true, ~, ~, ~, Speed_true] = ...
            compute_downstream_from_labels(true_labels, mesh, ...
            PL_LOS_func, PL_NLOS_func, ...
            P_tx_dBm, G_ant_dBi, Noise_dBm, B);

        [PL_knn, ~, ~, ~, Speed_knn] = ...
            compute_downstream_from_labels(est_kNN, mesh, ...
            PL_LOS_func, PL_NLOS_func, ...
            P_tx_dBm, G_ant_dBi, Noise_dBm, B);

        [PL_knmap, ~, ~, ~, Speed_knmap] = ...
            compute_downstream_from_labels(est_kNMAP, mesh, ...
            PL_LOS_func, PL_NLOS_func, ...
            P_tx_dBm, G_ant_dBi, Noise_dBm, B);

        mae_pl_buf_kNN(i) = mean(abs(PL_knn - PL_true));
        mae_pl_buf_kNMAP(i) = mean(abs(PL_knmap - PL_true));

        mae_spd_buf_kNN(i) = mean(abs(Speed_knn - Speed_true));
        mae_spd_buf_kNMAP(i) = mean(abs(Speed_knmap - Speed_true));
    end

    acc_kNN(idx) = mean(acc_buf_kNN) * 100;
    acc_kNMAP(idx) = mean(acc_buf_kNMAP) * 100;

    mae_pl_kNN(idx) = mean(mae_pl_buf_kNN);
    mae_pl_kNMAP(idx) = mean(mae_pl_buf_kNMAP);

    mae_spd_kNN(idx) = mean(mae_spd_buf_kNN);
    mae_spd_kNMAP(idx) = mean(mae_spd_buf_kNMAP);

    fprintf(['f = %d GHz | Acc kNN = %.2f | Acc kNMAP = %.2f ' ...
        '| MAE PL = %.2f/%.2f | MAE Spd = %.2f/%.2f\n'], ...
        freq_GHz, acc_kNN(idx), acc_kNMAP(idx), ...
        mae_pl_kNN(idx), mae_pl_kNMAP(idx), ...
        mae_spd_kNN(idx), mae_spd_kNMAP(idx));
end

%% Tạo thư mục lưu ảnh
output_folder = 'Ket_qua_hinh_anh';

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% Plot 1 - Accuracy
fig1 = figure('Color','w'); 
hold on;

plot(freq_list_GHz, acc_kNN, 'b-o', 'LineWidth', 2);
plot(freq_list_GHz, acc_kNMAP, 'r-s', 'LineWidth', 2);

xlabel('Tần số (GHz)');
ylabel('Độ chính xác (%)');

title('Độ chính xác theo tần số');

legend('kNN','kNMAP');
grid on;

saveas(fig1, fullfile(output_folder, 'Do_chinh_xac_theo_tan_so.png'));

%% Plot 2 - Path Loss MAE
fig2 = figure('Color','w'); 
hold on;

plot(freq_list_GHz, mae_pl_kNN, 'b-o', 'LineWidth', 2);
plot(freq_list_GHz, mae_pl_kNMAP, 'r-s', 'LineWidth', 2);

xlabel('Tần số (GHz)');
ylabel('Sai số MAE suy hao đường truyền (dB)');

title('Sai số suy hao đường truyền theo tần số');

legend('kNN','kNMAP');
grid on;

saveas(fig2, fullfile(output_folder, 'MAE_suy_hao_duong_truyen.png'));

%% Plot 3 - Speed MAE
fig3 = figure('Color','w'); 
hold on;

plot(freq_list_GHz, mae_spd_kNN, 'b-o', 'LineWidth', 2);
plot(freq_list_GHz, mae_spd_kNMAP, 'r-s', 'LineWidth', 2);

xlabel('Tần số (GHz)');
ylabel('Sai số MAE tốc độ (Mbps)');

title('Sai số tốc độ theo tần số');

legend('kNN','kNMAP');
grid on;

saveas(fig3, fullfile(output_folder, 'MAE_toc_do.png'));

disp('Da luu tat ca hinh anh vao thu muc: Ket_qua_hinh_anh');