function [PL_dB, P_rx_dBm, SNR_lin, SNR_dB, Speed_Mbps, dist_sr] = compute_downstream_from_labels(labels, mesh, PL_LOS_func, PL_NLOS_func, P_tx_dBm, G_ant_dBi, Noise_dBm, B)
    dist_sr = sqrt(mesh(:,1).^2 + mesh(:,2).^2);
    dist_sr(dist_sr == 0) = 1e-3;

    labels = logical(labels(:));
    PL_LOS_vec = PL_LOS_func(dist_sr);
    PL_NLOS_vec = PL_NLOS_func(dist_sr);

    PL_dB = zeros(size(dist_sr));
    PL_dB(labels==1) = PL_LOS_vec(labels==1);
    PL_dB(labels==0) = PL_NLOS_vec(labels==0);

    P_rx_dBm = (P_tx_dBm + G_ant_dBi) - PL_dB;
    SNR_lin  = 10.^((P_rx_dBm - Noise_dBm) / 10);
    SNR_dB   = 10 * log10(max(SNR_lin, 1e-12));
    Speed_Mbps = (B * log2(1 + SNR_lin)) / 1e6;
end
