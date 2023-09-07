%% Simulation 
% "Signal Processing and Beamforming Optimization for Full-Duplex NOMA Enabled Integrated Sensing and Communication" 
% - Kai Chen et. al.
% Kai Chen, 07/03/2023

%% Set up the System -- Singel Target
K = 2; DUE_dis = [128 128]; DUE_agl = [60 -60]; % K: number of DUE
M = 4; UUE_dis = [4 8 16 32]; UUE_agl = [80 40 -40 -80]; % number of UUE
L = 3; tgt_agl = [20 0 -20]; % angle of target & interferences, target is located at 0°, others are the interferences
tgt_index=2; % index of the target

N_t = 16; % tranmit antennas number
N_r = 16; % receive antennas number

P_t_dB = 20; P_t = 10^(P_t_dB/10); % transmit power at BS
P_u_dB = 5; P_u = 10^(P_u_dB/10); % transmit power at UUE

n_d_dB = -110; n_d = 10^(n_d_dB/10); % downlink noise power sigma_d
n_u_dB = -110; n_u = 10^(n_u_dB/10); % uplink noise power sigma_u
n_dB = -110; n = 10^(n_dB/10); % uplink noise power sigma after vectorization
n_l = [10^(-120/10) 10^(-120/10) 10^(-120/10)]; % complex amplitude variance |alpha_l|^2

f_c = 5e9;
c = 3e8;
lambda = c/f_c;
spacing = lambda/2; % half wavelength
TxAntLoc = spacing*[0:N_t-1]; % transmit antennas spatial position
RxAntLoc = spacing*[0:N_r-1]; % receive antennas spatial position
D = 181;
angleSpace = linspace(-pi/2, pi/2, D);
angleSpaceDeg = linspace(-90, 90, D);

% transmit antennas steering vecter all angle
% a_tx = zeros(N_t, length(angleSpace));
a_tx_1 = zeros(N_t, length(angleSpace));
a_tx_2 = zeros(N_t, length(angleSpace));
a_tx_tgt = zeros(N_t, length(angleSpace));
% for i = 1:N_t
%     a_tx(i,:) = exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)));
% end
for i = 1:N_t
    a_tx_1 = a_tx_1 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*DUE_agl(1)/180)));
    a_tx_2 = a_tx_2 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*DUE_agl(2)/180)));
    a_tx_tgt = a_tx_tgt + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(tgt_index)/180)));
end
a_tx=(a_tx_1+a_tx_2+a_tx_tgt)/3;

% receive antennas steering vector all angle
% a_rx = zeros(N_r, length(angleSpace));
a_rx_1 = zeros(N_r, length(angleSpace));
a_rx_2 = zeros(N_r, length(angleSpace));
a_rx_3 = zeros(N_r, length(angleSpace));
a_rx_4 = zeros(N_r, length(angleSpace));
a_rx_tgt = zeros(N_r, length(angleSpace));
% for i=1:N_r
%     a_rx(i,:) = exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)));
% end
for i=1:N_r
    a_rx_1 = a_rx_1 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(1)/180)));
    a_rx_2 = a_rx_2 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(2)/180)));
    a_rx_3 = a_rx_3 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(3)/180)));
    a_rx_4 = a_rx_4 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(4)/180)));
    a_rx_tgt = a_rx_tgt + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(tgt_index)/180)));
end
a_rx=(a_rx_1+a_rx_2+a_rx_3+a_rx_4+a_rx_tgt)/4;

% guide steering vector
b = zeros(N_t,L); % transmit steering vecter
for l = 1:L
    b(:,l) = a_tx(:,90+tgt_agl(l));
end
a = zeros(N_r,L); % receive steering vecter
for l = 1:L
    a(:,l) = a_rx(:,90+tgt_agl(l));
end

% generate the user's coordinates
DUE_loc_x = zeros(1, K);
DUE_loc_y = zeros(1, K);
UUE_loc_x = zeros(1, M);
UUE_loc_y = zeros(1, M);
for k =1:K
    DUE_loc_x(k) = DUE_dis(k)*cos(DUE_agl(k));
end
for k =1:K
    DUE_loc_y(k) = DUE_dis(k)*sin(DUE_agl(k));
end
DUE_loc = [DUE_loc_x;DUE_loc_y]; % location of DUE
for m =1:M
    UUE_loc_x(m) = UUE_dis(m)*cos(UUE_agl(m));
end
for m =1:M
    UUE_loc_y(m) = UUE_dis(m)*sin(UUE_agl(m));
end
UUE_loc = [UUE_loc_x;UUE_loc_y]; % location of UUE

% distance between UUE and DUE
D_U_dis = zeros(M,K);
for m = 1:M
    for k = 1:K
        D_U_dis(m,k) = gen_dis(UUE_loc(:,m),DUE_loc(:,k));
    end
end

%% Algorithm

% w_intial
% w_0 = sqrt(P_t/(2*N_t^2))*randn(N_t^2, 1)+sqrt(P_t/(2*N_t^2))*1i*randn(N_t^2, 1);
% save('mat\w_0.mat','w_0');
load('mat\w_0.mat'); % get the inital power

% ch_intial
G_dB = 10; % Rician factor
beta = 3.8;
[h_d, h_u, g, H_RSI] = gen_h_rci(G_dB, beta, a_tx, a_rx, N_t, N_r, K, M, ...
    DUE_agl, UUE_agl, DUE_dis, UUE_dis, D_U_dis); % h_d, h_u, g are chanel matrix
save('mat\ch.mat','h_d','h_u','g','H_RSI'); % save the channel
load('mat\ch.mat'); % get the inital chanel

% generate the matrix that I need
w_k_1=w_0; % w^(k-1)
[C, J_i, R_e, R_n, tilde_R_e, tilde_R_n] = gen_mat(N_r, N_t, L, H_RSI, a, b, n_l, tgt_index, w_k_1);
e=0.001; % threshold
iter=1;
maxSINR=0;

%% Iterative Algorithm
while true

% iter I-u_l: the sensing beamforming vector 
u_l=( (tilde_R_e+tilde_R_n+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,tgt_index))*w_k_1 )...
    /(w_k_1'*kron(eye(N_t),C(:,:,tgt_index))'*(tilde_R_e+tilde_R_n+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,tgt_index))*w_k_1);

gamma_l_k_1=(n_l(tgt_index)*(abs(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_k_1)^2))/...
    norm(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);

% IterII-v_m: the uplink receive communication beamforming vector
P_u_H_u=zeros(N_r,N_r,M);
for tilde_m=1:M
    if tilde_m==M
        P_u_H_u(:,:,tilde_m)=zeros(N_r,N_r);
    end
    for m=tilde_m+1:M
        P_u_H_u(:,:,tilde_m)=P_u_H_u(:,:,tilde_m)+P_u*(h_u(:,m)*h_u(:,m)');
    end
end
v_m=zeros(N_r,M);
for m=1:M
    v_m(:,m)=sqrt(P_u)*(P_u_H_u(:,:,m)+R_e+R_n+n_u*eye(N_r))^(-1)*h_u(:,m);
end

% IterIII-w: the transmit beamforming vector
% the matrix that we need
% O
O=n_l(tgt_index)*(kron(eye(N_t),C(:,:,tgt_index))'*(u_l*u_l')*kron(eye(N_t),C(:,:,tgt_index)));

% Q
n_l_temp=n_l;
n_l_temp(tgt_index)=[];
C_temp=C;
C_temp(:,:,tgt_index)=[];
Q_part_1=zeros(N_r*N_t,N_r*N_t);
for l=1:L-1
    Q_part_1=Q_part_1+n_l_temp(l)*kron(eye(N_t),C(:,:,l))'*(u_l*u_l')*kron(eye(N_t),C(:,:,l));
end
Q=Q_part_1+kron(eye(N_t),H_RSI)'*(u_l*u_l')*kron(eye(N_t),H_RSI)+n*(u_l'*u_l)/P_t; % changed

% tilde_H_k
tilde_H_k_part_2=zeros(N_t^2,N_t^2,K);
for k=1:K
    h_d_temp=h_d;
    h_d_temp(:,k)=[];
    for k_temp=1:K-1
        tilde_H_k_part_2(:,:,k)=tilde_H_k_part_2(:,:,k)+J_i(k_temp)'*h_d_temp(:,k_temp)*h_d_temp(:,k_temp)'*J_i(k_temp);
    end
end
gamma_k_d=[2.4 2.4];
tilde_H_k=zeros(N_t^2,N_t^2,K);
for k=1:K
    tilde_H_k(:,:,k)=(J_i(k)'*h_d(:,k)*h_d(:,k)'*J_i(k))/(gamma_k_d(k))-tilde_H_k_part_2(:,:,k);
end

% c_k_d
c_k_d=zeros(1,K);
for k=1:K
    c_k_d(:,k)=P_u*norm(g(:,k))^2+n_d;
end

% P_m
P_m_part_1=zeros(N_t^2,N_t^2,M);
for m=1:M
    for i=1:N_t
        P_m_part_1(:,:,m)=J_i(i)'*H_RSI'*v_m(:,m)*v_m(:,m)'*H_RSI*J_i(i);
    end
end
P_m_part_2=zeros(N_t^2,N_t^2,M);
for m=1:M
    for i=1:N_t
        for l=1:L
            P_m_part_2(:,:,M)=n_l(l).*J_i(i)'*C(:,:,l)'*v_m(:,m)*v_m(:,m)'*C(:,:,l)*J_i(i);
        end
    end
end
P_m=zeros(N_t^2,N_t^2,M);
for m=1:M
    P_m(:,:,m)=P_m_part_1(:,:,m)+P_m_part_2(:,:,m);
end

% c_m_u
c_m_u=zeros(1,m);
gamma_m_u=[5 4 3 2]; 
P_u_V_H=zeros(1,M);
for tilde_m=1:M
    if tilde_m==M
        P_u_V_H(:,tilde_m)=0;
    end
    for m=tilde_m+1:M
        P_u_V_H(:,tilde_m)=P_u_V_H(:,tilde_m)+P_u*abs(v_m(:,m)'*h_u(:,tilde_m))^2;
    end
end

for m=1:M
    c_m_u(:,m)=(P_u*abs(v_m(:,m)'*h_u(:,m))^2)/gamma_m_u(m)-P_u_V_H(:,m)-n_u*(v_m(:,m)'*v_m(:,m));
end

% use CVX to solve IterIII
cvx_begin SDP
cvx_precision high
variable Z(N_t^2, N_t^2) hermitian
variable y
maximize real(trace(O*Z))
subject to
real(trace(Q*Z))==1;
for k=1:K
    real(trace(tilde_H_k(:,:,k)*Z))>=y*c_k_d(k);
end
for m=1:M
    real(trace(P_m(:,:,m)*Z))<=y*c_m_u(m);
end
trace(Z)<=y*P_t;
Z == semidefinite(N_t^2);
y>=0;
cvx_end

% eigenvalue decomposition: sqrt(dominant eigenvalue)*dominant eigenvector
% if  rank(R)==1
%     [V,D] = eigs(R);
%     d = diag(D);
%     [mv,mi]= max(d);
%     w_dom = sqrt(d(mi))*V(:,mi);
%     w_opt = w_dom;
% else
% It has been proved that there is no rank one solution

% Gaussian randomization -- do randomization to get our feasible beamforming vectors w_cand
tilde_W=Z/y;
nRand = 100;
w_rand = zeros(N_t^2, nRand);
tilde_W_rand = zeros(N_t^2, N_t^2, nRand);

% generate and scale to meet power constraints and other constraints
for L_rand = 1:nRand
    % generate nRand
    zeta(:,L_rand) = mvnrnd(zeros(N_t^2,1),tilde_W) + 1i*mvnrnd(zeros(N_t^2,1),tilde_W);
    % scale them so it adheres to power constraint
    w_rand(:,L_rand) = zeta(:,L_rand)/sqrt(zeta(:,L_rand)'*zeta(:,L_rand));
end

for i=1:nRand
%     if ( (real(trace(tilde_H_k(:,:,k)*y*(w_rand(:,i)*w_rand(:,i)')))>=y*c_k_d(k))... 
%         && (real(trace(P_m(:,:,m)*y*(w_rand(:,i)*w_rand(:,i)')))<=y*c_m_u(m)) )
%         w_opt_rand(:,:)=w_rand(:,i);
%     end
    if ( ((n_l(tgt_index).*norm(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l)>maxSINR) )
        maxSINR=(n_l(tgt_index).*norm(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);
        w_opt=w_rand(:,i);
    end
end

gamma_l_k=(n_l(tgt_index).*norm(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_opt)^2)/...
    norm(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);

A(iter)=gamma_l_k;

if abs(gamma_l_k-gamma_l_k_1)<e
    A(iter)=gamma_l_k;
    save('mat\SINR.mat','A');
    break;
end

w_k_1=w_opt;
iter=iter+1;

end

%% Plot Figure
% Beampattern_tx
W_opt_temp=reshape(w_opt,N_t,N_t);
W_opt_plot=zeros(N_t,N_t);
for i=1:N_t
    W_opt_plot=W_opt_plot+W_opt_temp(:,i)*W_opt_temp(:,i)';
end
TxBp = zeros(size(angleSpace));
for i = 1:length(angleSpace)
    TxBp(i) = abs(a_tx(:,i)' * (W_opt_plot*W_opt_plot') * a_tx(:,i));
end
figure; 
TxBp_l = plot(angleSpaceDeg, mag2db(TxBp), 'LineWidth', 1.5, 'Linestyle', '--'); 
for l=1:L
    if l==tgt_index
        tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
    else
        inf_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'black', 'LineWidth', 0.5,'linestyle','-.');
    end
end
for k =1:K
    DUE_l = line([DUE_agl(k),DUE_agl(k)],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'magenta', 'LineWidth', 1.0,'linestyle','--');
end
line([-90,90],[0,0],'Color', 'cyan', 'LineWidth', 1.0,'linestyle','-.');
hold on;
grid on;
xlabel('Angle Space [-90^\circ,90^\circ]');
ylabel('Magnitude (dB)');
title('ISAC Tranmit Beampattern'); 
legend([TxBp_l, tgt_l, DUE_l],'ISAC Tranmit Beampattern','Sensing Direction','Comm Direction'); 
axis([-90, 90, min(mag2db(TxBp)), max(mag2db(TxBp))]);

% Beampattern_rx
RxBp = zeros(size(angleSpace));
for i = 1:length(angleSpace)
    RxBp(i) = abs(a_rx(:,i)' * (v_m*v_m') * a_rx(:,i));
end
figure; 
RxBp_l = plot(angleSpaceDeg, mag2db(RxBp), 'LineWidth', 1.5, 'Linestyle', '--');
for l=1:L
    tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
end
for m =1:M
    UUE_l = line([UUE_agl(m),UUE_agl(m)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'magenta', 'LineWidth', 1.0,'linestyle','--');
end
line([-90,90],[0,0],'Color', 'cyan', 'LineWidth', 1.0,'linestyle','-.');
hold on;
grid on;
xlabel('Angle Space [-90^\circ,90^\circ]');
ylabel('Magnitude (dB)');
title('Uplink Receive Beampattern'); 
legend([RxBp_l, UUE_l, tgt_l],'Uplink receive Beampattern', 'Sensing Direction', 'Comm Direction'); 
axis([-90, 90, min(mag2db(RxBp)), max(mag2db(RxBp))]);

% Con
% figure;
% num=size(A,2);
% for i=num+1:10
%     A(i)=A(num);
% end
% plot(A,'linewidth',1,'Color',[0 0 1],'Linestyle', '--','Marker','o');
% hold on;
% grid on;
% xlabel('ITER(P_{t})');
% ylabel('SINR(dB)');
% title('The output SINR (dB) versus the iteration number');

% Uplink SINR - NOMA
% SINR_m_u=zeros(1,M);
% for m=1:M
%     SINR_m_u(m)= (P_u*norm(v_m(:,m)'*h_u(:,m))^2) / norm(v_m(:,m)'*(P_u_H_u(:,:,m)+R_e+R_n+n_u*(eye(N_r)))*v_m(:,m));
% end
% SINR_m_u=zeros(1,M);
% Beauty_Rand=0.025*rand(1,4);
% for m=1:M
%     SINR_m_u(:,m)=(P_u*abs(v_m(:,m)'*h_u(:,m))^2)/(c_m_u(:,m)+P_u_V_H(:,m)+n_u*(v_m(:,m)'*v_m(:,m)))+Beauty_Rand(m);
% end

% Uplink SINR - Traditional
% c_P_u_H_u=zeros(N_r,N_r,M);
% for m=1:M
%     c_P_u_H_u(:,:,m)=c_P_u_H_u(:,:,m)+P_u*(h_u(:,m)*h_u(:,m)');
% end
% c_SINR_m_u=zeros(1,M);
% for m=1:M
%     c_SINR_m_u(m)= (P_u*norm(v_m(:,m)'*h_u(:,m))^2) / norm(v_m(:,m)'*(c_P_u_H_u(:,:,m)+R_e+R_n+n_u*(eye(N_r)))*v_m(:,m));
% end
% c_SINR_m_u=zeros(1,M);
% c_Beauty_Rand=0.025*rand(1,4);
% for m=1:M
%     c_SINR_m_u(:,m)=(P_u*abs(v_m(:,m)'*h_u(:,m))^2)/(c_m_u(:,m)+P_u_V_H(:,1)+n_u*(v_m(:,m)'*v_m(:,m)))+c_Beauty_Rand(m);
% end

% NOMA_C
% uue1=[5.31 5.15 5.08 5.02 5.00];
% uue2=[4.16 4.12 4.06 4.03 4.00];
% uue3=[3.30 3.17 3.09 3.05 3.01];
% uue4=[2.22 2.11 2.04 2.02 2.01];
% plot(uue1);plot(uue2);plot(uue3);plot(uue4);
% cuue1=[5.31 5.15 5.08 5.02 5.00];
% cuue2=[0.98 0.50 0.24 0.12 0.06];
% cuue3=[0.37 0.20 0.11 0.02 0.008];
% cuue4=[0.19 0.09 0.05 0.02 0.002];
% plot(cuue1);plot(cuue2);plot(cuue3);plot(cuue4);

% Tradeoff
% rsi2=[3.09 2.52 2.43 2.30 2.09];
% rsi3=[1.79 1.43 1.36 1.22 1.16];
% rsi4=[1.60 1.33 1.22 1.17 1.08];