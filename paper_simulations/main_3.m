%% ============================================================
%  MO PHONG MANG KHONG DAY NGOAI TRO - PPP BUILDINGS (CHUAN BAO)
%  Ket hop:
%   - Mo hinh PPP (bai bao IEEE)
%   - Node di chuyen (bai 1)
%   - Hieu suat mang (RSSI, Throughput, Loss, Delay)
%% ============================================================

clc; clear; close all;
rng(42);

%% ===== THAM SO =====
R        = 100;        % Ban kinh khu vuc
freq     = 2.4e9;
Pt_dBm   = 20;
lambda_b = 0.0008;     % mat do vat can PPP
Lmax     = 20;         % kich thuoc vat can
nNodes   = 5;
nSteps   = 150;
sensitivity = -75;

AP = [0,0];

%% ===== SINH VAT CAN (HINH CHU NHAT) =====
N_build = poissrnd(lambda_b * pi * R^2);

buildings = zeros(N_build,4); % [x y w h]

for i=1:N_build
    x = (rand*2-1)*R;
    y = (rand*2-1)*R;
    w = rand*Lmax;
    h = rand*Lmax;
    buildings(i,:) = [x y w h];
end

%% ===== TAO SHADOW (DA GIAC) =====
shadows = cell(N_build,1);

for i=1:N_build
    x = buildings(i,1);
    y = buildings(i,2);
    w = buildings(i,3);
    h = buildings(i,4);

    rect = [
        x-w/2 y-h/2;
        x+w/2 y-h/2;
        x+w/2 y+h/2;
        x-w/2 y+h/2
    ];

    scale = 3;
    shadow = rect * scale;
    shadows{i} = shadow;
end

%% ===== KHOI TAO NODE =====
nodeX = (rand(nNodes,1)*2-1)*R;
nodeY = (rand(nNodes,1)*2-1)*R;

vx = (rand(nNodes,1)-0.5)*2;
vy = (rand(nNodes,1)-0.5)*2;

colors = lines(nNodes);

%% ===== LUU LICH SU =====
rssi_hist = zeros(nNodes,nSteps);
tput_hist = zeros(nNodes,nSteps);
loss_hist = zeros(nNodes,nSteps);
delay_hist= zeros(nNodes,nSteps);

%% ===== FIGURE =====
figure('Position',[50 50 1400 800],'Color','white');

%% ============================================================
%  VONG LAP MO PHONG
%% ============================================================
for step=1:nSteps

    % ===== DI CHUYEN =====
    nodeX = nodeX + vx;
    nodeY = nodeY + vy;

    % Gioi han trong vung tron
    for i=1:nNodes
        if norm([nodeX(i),nodeY(i)]) > R
            vx(i) = -vx(i);
            vy(i) = -vy(i);
        end
    end

    %% ===== TINH TOAN =====
    for n=1:nNodes

        px = nodeX(n);
        py = nodeY(n);

        d = max(sqrt(px^2+py^2),1);

        % ===== KIEM TRA NLOS =====
        nlos = false;
        for s=1:length(shadows)
            poly = shadows{s};
            if inpolygon(px,py,poly(:,1),poly(:,2))
                nlos = true;
                break;
            end
        end

        % ===== PATH LOSS =====
        lambda_w = 3e8/freq;
        PL0 = 20*log10(4*pi/lambda_w);

        if nlos
            n_path = 3.5;
        else
            n_path = 2.0;
        end

        PL = PL0 + 10*n_path*log10(d);
        rssi = Pt_dBm - PL;

        rssi_hist(n,step) = rssi;

        %% ===== HIEU SUAT =====
        if rssi >= -60
            tput = 60;
            loss = 1;
            delay = 5;
        elseif rssi >= -67
            tput = 40;
            loss = 5;
            delay = 10;
        elseif rssi >= -72
            tput = 20;
            loss = 20;
            delay = 30;
        elseif rssi >= sensitivity
            tput = 5;
            loss = 50;
            delay = 80;
        else
            tput = 0;
            loss = 100;
            delay = 200;
        end

        tput_hist(n,step) = tput;
        loss_hist(n,step) = loss;
        delay_hist(n,step)= delay;
    end

    clf;

    %% =====================================================
    %  O 1: BAN DO + VAT CAN + DUONG TRUYEN
    %% =====================================================
    subplot(2,3,[1,2]);
    hold on;

    % Ve vung tron
    th = linspace(0,2*pi,200);
    plot(R*cos(th),R*sin(th),'k--');

    % Ve vat can
    for i=1:N_build
        x = buildings(i,1);
        y = buildings(i,2);
        w = buildings(i,3);
        h = buildings(i,4);

        rectangle('Position',[x-w/2,y-h/2,w,h],...
            'FaceColor',[0.3 0.3 0.3],'EdgeColor','k');
    end

    % Ve AP
    plot(0,0,'r^','MarkerSize',12,'MarkerFaceColor','r');

    % Ve node + duong truyen
    for n=1:nNodes

        rssi = rssi_hist(n,step);

        if rssi >= -67
            col = 'g'; ls='-';
        elseif rssi >= sensitivity
            col = [1 0.5 0]; ls='--';
        else
            col = 'r'; ls=':';
        end

        plot([0 nodeX(n)],[0 nodeY(n)],...
            'Color',col,'LineStyle',ls,'LineWidth',2);

        scatter(nodeX(n),nodeY(n),120,colors(n,:),'filled');
    end

    title(sprintf('Mo hinh PPP ngoai tro (Step %d)',step));
    axis equal; xlim([-R R]); ylim([-R R]);

    %% =====================================================
    %  O 2: THROUGHPUT
    %% =====================================================
    subplot(2,3,3); hold on;
    for n=1:nNodes
        plot(1:step,tput_hist(n,1:step),'LineWidth',2);
    end
    title('Throughput (Mbps)');
    ylim([0 70]); grid on;

    %% =====================================================
    %  O 3: RSSI
    %% =====================================================
    subplot(2,3,4); hold on;
    for n=1:nNodes
        plot(1:step,rssi_hist(n,1:step),'LineWidth',2);
    end
    yline(sensitivity,'r--');
    title('RSSI (dBm)');
    ylim([-120 -20]); grid on;

    %% =====================================================
    %  O 4: PACKET LOSS
    %% =====================================================
    subplot(2,3,5); hold on;
    for n=1:nNodes
        plot(1:step,loss_hist(n,1:step),'LineWidth',2);
    end
    title('Packet Loss (%)');
    ylim([0 100]); grid on;

    %% =====================================================
    %  O 5: DELAY
    %% =====================================================
    subplot(2,3,6); hold on;
    for n=1:nNodes
        plot(1:step,delay_hist(n,1:step),'LineWidth',2);
    end
    title('Delay (ms)');
    ylim([0 200]); grid on;

    drawnow;
end