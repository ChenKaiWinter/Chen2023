function [H_ric_d, H_ric_u, H_ric_u_d, H_RSI] = gen_h_rci(G_dB, beta, a_tx, a_rx, N_t, N_r, K, M, ...
    DUE_agl, UUE_agl, DUE_dis, UUE_dis, D_U_dis)

% doesn't matter with the mento carlo method
% Rician factor
G = 10^(G_dB/10);
tx_ric_vec = zeros(N_t,K);
rx_ric_vec = zeros(N_r,M);
for k = 1:K
    tx_ric_vec(:,k) = a_tx(:,90+DUE_agl(k));
end
for m = 1:M
    rx_ric_vec(:,m) = a_rx(:,90+UUE_agl(m));
end

% ray part
H_ray_d = zeros(N_t,K);
for k = 1:K
    H_ray_d(:,k) = (randn(N_t,1)+1i*randn(N_t,1))/sqrt(2);
end
H_ray_u = zeros(N_r,M);
for m = 1:M
    H_ray_u(:,m) = (randn(N_r,1)+1i*randn(N_r,1))/sqrt(2);
end
H_ray_u_d = zeros(M,K);
for k = 1:K
    H_ray_u_d(:,k) = (randn(M,1)+1i*randn(M,1))/sqrt(2);
end

% ric part
H_ric_d = zeros(N_t,K);
for k = 1 : K
    H_ric_d(:,k) = sqrt(DUE_dis(k)^(-beta)) * (sqrt(G/(G+1)) * tx_ric_vec(:,k) + sqrt(1/(G+1)) * H_ray_d(:,k));
end
H_ric_u = zeros(N_r,M);
for m = 1 : M
    H_ric_u(:,m) = sqrt(UUE_dis(m)^(-beta)) * (sqrt(G/(G+1)) * rx_ric_vec(:,m) + sqrt(1/(G+1)) * H_ray_u(:,m));
end
H_ric_u_d = zeros(M,K);
for m = 1 : M
    for k = 1 : K
        H_ric_u_d(m,k) = sqrt(D_U_dis(m,k)^(-beta)) * (H_ray_u_d(m,k));
    end
end

% H_RSI
eta_RIS_dB=3; % ...dBm
eta_RIS=10^(eta_RIS_dB/10);
hat_H_RSI=(randn(N_r,N_t)+1i*randn(N_r,N_t))/sqrt(2);
H_RSI=sqrt(eta_RIS)*hat_H_RSI;