%% ================== PRE-RUN ==================
clear; close all;
warning('off','all')
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

%% ================== THAM SỐ ==================
lambda = 0.00075;
Lmax   = 20;
R      = 150;
N_data = 80;
K      = 3;

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

% Nhiễu nền
N0_W      = k_B * T_K * B;
Noise_dBm = 10*log10(N0_W) + 30 + NF_dB;

% Tham số ăng-ten (3GPP TR 38.901 UMa)
h_BS = 25;
h_UT = 1.5;
h_E  = 1.0;

h_BS_eff = h_BS - h_E;
h_UT_eff = h_UT - h_E;
d_BP     = 4 * h_BS_eff * h_UT_eff * freq / c;

fprintf('=== THAM SO 3GPP TR 38.901 ===\n');
fprintf('Noise floor: %.2f dBm\n', Noise_dBm);
fprintf('Breakpoint distance d_BP = %.1f m\n', d_BP);

%% ================== HÀM TÍNH PATH LOSS 3GPP ==================
PL_LOS_func = @(d) (d < d_BP) .* (28.0 + 22*log10(max(d,1e-3)) + 20*log10(freq_GHz)) + ...
                   (d >= d_BP).* (28.0 + 40*log10(max(d,1e-3)) + 20*log10(freq_GHz) ...
                                  - 9*log10(d_BP^2 + (h_BS - h_UT)^2));

PL_NLOS_func = @(d) max(PL_LOS_func(d), ...
                        13.54 + 39.08*log10(max(d,1e-3)) + 20*log10(freq_GHz) - 0.6*(h_UT - 1.5));

%% ================== MÔ PHỎNG 1 LẦN ==================
[lambda, Lmax, buildings] = chicago_buildings(R);

L            = 2*R;
total_shadow = generate_shadows(R, buildings);
data         = generate_data(N_data, total_shadow, R);  

figure(); 
plot_estimation(L, R, total_shadow, data, [], [], 'Data and shadows', buildings);

%% Mesh
N_mesh = 12;
[x, y] = meshgrid(linspace(-R,R,N_mesh), linspace(-R,R,N_mesh));
in_disk = x.^2 + y.^2 <= R^2;
mesh    = [x(in_disk), y(in_disk)];

%% Estimation
estimated_kNN   = kNN_estimator(data, mesh, K);  
estimated_kNMAP = kNMAP_estimator(data, lambda, Lmax, mesh, K);  

%% Ground truth
in_shadow = false(size(mesh,1), 1); 
if ismethod(total_shadow, 'regions')
    shadow_regions = regions(total_shadow);
    for k = 1:length(shadow_regions)
        v         = shadow_regions(k).Vertices;
        in_shadow = in_shadow | inpolygon(mesh(:,1), mesh(:,2), v(:,1), v(:,2));
    end
else
    v         = total_shadow.Vertices;
    in_shadow = inpolygon(mesh(:,1), mesh(:,2), v(:,1), v(:,2));
end
true_labels = ~in_shadow;

%% Precision
precision_kNN   = estimation_precision(mesh, estimated_kNN,   total_shadow);
precision_kNMAP = estimation_precision(mesh, estimated_kNMAP, total_shadow);

fprintf('--- SINGLE RUN ---\n');
fprintf('kNN:   %.2f %%\n', precision_kNN*100);
fprintf('kNMAP: %.2f %%\n', precision_kNMAP*100);

%% Plot estimation
figure;
subplot(1,2,1);
plot_estimation(L, R, total_shadow, data, mesh, estimated_kNN,   'kNN',   buildings);
subplot(1,2,2);
plot_estimation(L, R, total_shadow, data, mesh, estimated_kNMAP, 'kNMAP', buildings);

%% Error map
figure;
plot_estimation(L, R, total_shadow, data, mesh, estimated_kNMAP, 'Error kNMAP', buildings);
hold on;
err = estimated_kNMAP(:) ~= true_labels(:);
plot(mesh(err,1), mesh(err,2), 'r.', 'MarkerSize', 15);
legend('show'); hold off;

figure;
plot_estimation(L, R, total_shadow, data, mesh, estimated_kNN, 'Error kNN', buildings);
hold on;
err2 = estimated_kNN(:) ~= true_labels(:);
plot(mesh(err2,1), mesh(err2,2), 'r.', 'MarkerSize', 15);
legend('show'); hold off;

%% ================== HEATMAP TỐC ĐỘ (SINGLE RUN - 3GPP) ==================
dist_sr = sqrt(mesh(:,1).^2 + mesh(:,2).^2);
dist_sr(dist_sr == 0) = 1e-3;

PL_sr = zeros(size(dist_sr));
PL_sr(true_labels == 1) = PL_LOS_func( dist_sr(true_labels == 1));
PL_sr(true_labels == 0) = PL_NLOS_func(dist_sr(true_labels == 0));

P_rx_sr  = (P_tx_dBm + G_ant_dBi) - PL_sr;
SNR_sr   = 10.^((P_rx_sr - Noise_dBm) / 10);
Speed_sr = (B * log2(1 + SNR_sr)) / 1e6;

figure('Color','w');
hold on;
scatter(mesh(:,1), mesh(:,2), 50, Speed_sr, 'filled', 'MarkerEdgeColor', 'k');
plot(0, 0, 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'y');
plot_estimation(L, R, total_shadow, data, [], [], '', buildings);
title('Heatmap toc do (LOS + NLOS) - 3GPP TR 38.901');
axis equal; xlim([-R R]); ylim([-R R]);
grid on;
cb = colorbar; cb.Label.String = 'Toc do (Mbps)';
colormap(jet);
hold off;

%% ================== MONTE CARLO ==================
N_run = 10;

results_kNN   = zeros(N_run, 1);
results_kNMAP = zeros(N_run, 1);

Speed_LOS_all  = [];
Speed_NLOS_all = [];
PL_LOS_all     = [];
PL_NLOS_all    = [];
d_LOS_all      = [];
d_NLOS_all     = [];

% Lưu tốc độ trung bình từng lần để vẽ đường theo run
mean_LOS_per_run  = zeros(N_run, 1);
mean_NLOS_per_run = zeros(N_run, 1);

for i = 1:N_run

    [lambda, Lmax, buildings] = chicago_buildings(R);
    total_shadow = generate_shadows(R, buildings);
    data         = generate_data(N_data, total_shadow, R);  

    % Ground truth
in_shadow = false(size(mesh,1), 1);
if ismethod(total_shadow, 'regions')
    shadow_regions = regions(total_shadow);
    for kk = 1:length(shadow_regions)
        v         = shadow_regions(kk).Vertices;
        in_shadow = in_shadow | inpolygon(mesh(:,1), mesh(:,2), v(:,1), v(:,2));
    end
else
    v         = total_shadow.Vertices;
    in_shadow = inpolygon(mesh(:,1), mesh(:,2), v(:,1), v(:,2));
end
true_labels = ~in_shadow;

    % Estimation
    est_kNN   = kNN_estimator(data, mesh, K);  
    est_kNMAP = kNMAP_estimator(data, lambda, Lmax, mesh, K);  

    results_kNN(i)   = estimation_precision(mesh, est_kNN,   total_shadow);
    results_kNMAP(i) = estimation_precision(mesh, est_kNMAP, total_shadow);

    % ===== TÍNH SPEED + PATH LOSS (3GPP TR 38.901) =====
    dist = sqrt(mesh(:,1).^2 + mesh(:,2).^2);
    dist(dist == 0) = 1e-3;

    PL_LOS_vec  = PL_LOS_func(dist);
    PL_NLOS_vec = PL_NLOS_func(dist);

    PL_mc = zeros(size(dist));
    PL_mc(true_labels == 1) = PL_LOS_vec( true_labels == 1);
    PL_mc(true_labels == 0) = PL_NLOS_vec(true_labels == 0);

    P_rx  = (P_tx_dBm + G_ant_dBi) - PL_mc;
    SNR   = 10.^((P_rx - Noise_dBm) / 10);
    Speed = (B * log2(1 + SNR)) / 1e6;

    % Lưu tốc độ tích lũy
    Speed_LOS_all  = [Speed_LOS_all;  Speed(true_labels == 1)];
    Speed_NLOS_all = [Speed_NLOS_all; Speed(true_labels == 0)];

    % Lưu tốc độ trung bình từng lần
    mean_LOS_per_run(i)  = mean(Speed(true_labels == 1));
    mean_NLOS_per_run(i) = mean(Speed(true_labels == 0));

    % Lưu path loss và khoảng cách
    PL_LOS_all  = [PL_LOS_all;  PL_LOS_vec( true_labels == 1)];
    PL_NLOS_all = [PL_NLOS_all; PL_NLOS_vec(true_labels == 0)];
    d_LOS_all   = [d_LOS_all;   dist(true_labels == 1)];
    d_NLOS_all  = [d_NLOS_all;  dist(true_labels == 0)];

end

%% ================== BIỂU ĐỒ PHÂN BỐ TỐC ĐỘ (grouped bar) ==================
nbins_s    = 15;
all_speeds = [Speed_LOS_all; Speed_NLOS_all];
edges_S    = linspace(min(all_speeds), max(all_speeds), nbins_s + 1);

[cnt_S_LOS,  ~] = histcounts(Speed_LOS_all,  edges_S);
[cnt_S_NLOS, ~] = histcounts(Speed_NLOS_all, edges_S);
bin_S = (edges_S(1:end-1) + edges_S(2:end)) / 2;

figure('Color','w');
b = bar(bin_S, [cnt_S_LOS(:), cnt_S_NLOS(:)], 'grouped', 'BarWidth', 1.2);
b(1).FaceColor = [0 0.45 0.74];
b(2).FaceColor = [0.85 0.33 0.1];
b(1).FaceAlpha = 0.85;
b(2).FaceAlpha = 0.85;
xlabel('Toc do (Mbps)', 'FontSize', 12);
ylabel('So diem',       'FontSize', 12);
title('Phan bo toc do LOS vs NLOS (Monte Carlo)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'LOS','NLOS'}, 'Location', 'best', 'FontSize', 11);
grid on; box on;

%% ================== BIỂU ĐỒ TỐC ĐỘ TRUNG BÌNH THEO TỪNG LẦN MONTE CARLO ==================
figure('Color','w');
hold on;

% Đường tốc độ trung bình từng lần
plot(1:N_run, mean_LOS_per_run,  'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0 0.45 0.74],   'DisplayName', 'LOS');
plot(1:N_run, mean_NLOS_per_run, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0.85 0.33 0.1], 'DisplayName', 'NLOS');

% Đường trung bình tổng (ngang)
yline(mean(mean_LOS_per_run),  'b--', 'LineWidth', 1.5, ...
      'DisplayName', sprintf('TB LOS = %.1f Mbps',  mean(mean_LOS_per_run)));
yline(mean(mean_NLOS_per_run), 'r--', 'LineWidth', 1.5, ...
      'DisplayName', sprintf('TB NLOS = %.1f Mbps', mean(mean_NLOS_per_run)));

xlabel('Lan chay (run)',        'FontSize', 12);
ylabel('Toc do trung binh (Mbps)', 'FontSize', 12);
title('Toc do trung binh LOS vs NLOS theo tung lan Monte Carlo', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
xlim([1 N_run]);
grid on; box on;
hold off;

%% ================== BIỂU ĐỒ BAR TỐC ĐỘ TRUNG BÌNH TỔNG ==================
figure('Color','w');
mean_LOS_total  = mean(mean_LOS_per_run);
mean_NLOS_total = mean(mean_NLOS_per_run);

b3 = bar([mean_LOS_total, mean_NLOS_total], 'BarWidth', 0.5, 'FaceColor', 'flat');
b3.CData(1,:) = [0 0.45 0.74];
b3.CData(2,:) = [0.85 0.33 0.1];
b3.FaceAlpha  = 0.85;

set(gca, 'XTickLabel', {'LOS','NLOS'}, 'FontSize', 12);
ylabel('Toc do trung binh (Mbps)', 'FontSize', 12);
title('Toc do trung binh LOS vs NLOS (Tong hop Monte Carlo)', 'FontSize', 14, 'FontWeight', 'bold');
grid on; box on;
ylim([0, max(mean_LOS_total, mean_NLOS_total) * 1.25]);
text(1, mean_LOS_total  + 0.5, sprintf('%.2f Mbps', mean_LOS_total),  ...
     'HorizontalAlignment','center', 'FontSize', 12, 'FontWeight', 'bold');
text(2, mean_NLOS_total + 0.5, sprintf('%.2f Mbps', mean_NLOS_total), ...
     'HorizontalAlignment','center', 'FontSize', 12, 'FontWeight', 'bold');

%% ================== MONTE CARLO - PRECISION ==================
mean_kNN   = mean(results_kNN);
mean_kNMAP = mean(results_kNMAP);

figure('Color','w');
b2 = bar([mean_kNN*100 mean_kNMAP*100], 'BarWidth', 0.6, 'FaceColor', 'flat');
b2.CData(1,:) = [0.3 0.75 0.9];
b2.CData(2,:) = [0.85 0.33 0.1];
set(gca, 'XTickLabel', {'kNN','kNMAP'}, 'FontSize', 12);
ylabel('Precision (%)', 'FontSize', 12);
title('Monte Carlo Average Precision', 'FontSize', 14, 'FontWeight', 'bold');
grid on; box on;
text(1, mean_kNN*100,   sprintf('%.2f%%', mean_kNN*100),   ...
     'HorizontalAlignment','center', 'FontSize', 11, 'FontWeight', 'bold');
text(2, mean_kNMAP*100, sprintf('%.2f%%', mean_kNMAP*100), ...
     'HorizontalAlignment','center', 'FontSize', 11, 'FontWeight', 'bold');

%% ================== BIỂU ĐỒ 1: PL vs KHOẢNG CÁCH (3GPP) ==================
d_range        = linspace(1, R, 500);
PL_LOS_theory  = PL_LOS_func(d_range);
PL_NLOS_theory = PL_NLOS_func(d_range);

figure('Color','w'); hold on;

scatter(d_LOS_all,  PL_LOS_all,  8, [0 0.45 0.74],   'filled', ...
        'MarkerFaceAlpha', 0.15, 'DisplayName', 'LOS (data)');
scatter(d_NLOS_all, PL_NLOS_all, 8, [0.85 0.33 0.1], 'filled', ...
        'MarkerFaceAlpha', 0.15, 'DisplayName', 'NLOS (data)');
plot(d_range, PL_LOS_theory,  'b-', 'LineWidth', 2.5, 'DisplayName', 'LOS (3GPP ly thuyet)');
plot(d_range, PL_NLOS_theory, 'r-', 'LineWidth', 2.5, 'DisplayName', 'NLOS (3GPP ly thuyet)');
xline(d_BP, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('d_{BP} = %.0f m', d_BP));

xlabel('Khoang cach (m)', 'FontSize', 12);
ylabel('Path Loss (dB)',   'FontSize', 12);
title('Suy hao duong truyen LOS vs NLOS - 3GPP TR 38.901 (Monte Carlo)', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
grid on; box on;
xlim([0 R]);
hold off;

%% ================== BIỂU ĐỒ 2: HISTOGRAM PHÂN BỐ PL ==================
nbins_pl = 15;

edges = linspace(min([PL_LOS_all; PL_NLOS_all]), ...
                 max([PL_LOS_all; PL_NLOS_all]), nbins_pl+1);

[cnt_LOS,  ~] = histcounts(PL_LOS_all,  edges);
[cnt_NLOS, ~] = histcounts(PL_NLOS_all, edges);

bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

figure('Color','w'); hold on;

%  CỘT NGANG (cùng phía, không đối xứng)
bar(bin_centers, cnt_LOS,  0.4, 'FaceColor', [0 0.45 0.74], ...
     'DisplayName', 'LOS');

bar(bin_centers, cnt_NLOS, 0.4, 'FaceColor', [0.85 0.33 0.1], ...
     'DisplayName', 'NLOS');

xlabel('So diem', 'FontSize', 12);
ylabel('Path Loss (dB)', 'FontSize', 12);

title('Phan bo Path Loss LOS vs NLOS - 3GPP TR 38.901', ...
      'FontSize', 14, 'FontWeight', 'bold');

legend('Location','best');
grid on; box on;
hold off;