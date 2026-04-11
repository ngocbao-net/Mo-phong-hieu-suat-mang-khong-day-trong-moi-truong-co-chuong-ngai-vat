%% ============================================================
% kNN vs kNMAP - Lambda Comparison 
% ============================================================

clear; close all;
warning('off','all');

% Random mỗi lần chạy khác nhau
rng('shuffle');

%% PARAMETERS
R = 150;           % Cell radius
Lmax = 20;         % Building size
N_data = 80;       % Số điểm training
N_mesh = 12;       % Grid
k = 3;             % kNN

num_lambda = 10;    
num_runs = 500;    % Monte Carlo (tăng nếu muốn mượt hơn)

% Lambda tăng dần
lambda_list = linspace(0.0005, 0.0035, num_lambda);

%% CREATE MESH
[x, y] = meshgrid(linspace(-R,R,N_mesh), linspace(-R,R,N_mesh));
mask = x.^2 + y.^2 <= R^2;
mesh = [x(mask), y(mask)];

%% STORAGE
precision_kNN   = zeros(num_lambda,1);
precision_kNMAP = zeros(num_lambda,1);

%% MAIN SIMULATION
for i = 1:num_lambda
    
    lambda = lambda_list(i);
    
    acc_kNN = 0;
    acc_kNMAP = 0;
    
    for run = 1:num_runs
        
        % --- RANDOM BUILDINGS ---
        buildings = PPP_buildings(lambda, Lmax, R);
        
        % --- SHADOW ---
        shadow = generate_shadows(R, buildings);
        
        % --- DATA ---
        data = generate_data(N_data, shadow, R);
        
        % --- ESTIMATION ---
        est_kNN   = kNN_estimator(data, mesh, k);
        est_kNMAP = kNMAP_estimator(data, lambda, Lmax, mesh, k);
        
        % --- PRECISION ---
        acc_kNN   = acc_kNN   + estimation_precision(mesh, est_kNN, shadow);
        acc_kNMAP = acc_kNMAP + estimation_precision(mesh, est_kNMAP, shadow);
        
    end
    
    % --- AVERAGE ---
    precision_kNN(i)   = (acc_kNN / num_runs) * 100;
    precision_kNMAP(i) = (acc_kNMAP / num_runs) * 100;
    
end

%% ============================================================
% PLOT 
% ============================================================

figure;
plot(lambda_list, precision_kNN, '-o','LineWidth',2); hold on;
plot(lambda_list, precision_kNMAP, '-s','LineWidth',2);

xlabel('\lambda');
ylabel('Precision (%)');
title('kNN vs kNMAP (Monte Carlo Simulation)');

legend('kNN','kNMAP','Location','best');
grid on;

% TRỤC Y
ylim([0 100]);

% TRỤC X
xticks(lambda_list);
xticklabels(string(num2str(lambda_list','%.4f')));

% HIỂN THỊ GIÁ TRỊ
for i = 1:length(lambda_list)
    
    text(lambda_list(i), precision_kNN(i)+1, ...
        sprintf('%.1f', precision_kNN(i)), ...
        'HorizontalAlignment','center');
    
    text(lambda_list(i), precision_kNMAP(i)+1, ...
        sprintf('%.1f', precision_kNMAP(i)), ...
        'HorizontalAlignment','center');
end