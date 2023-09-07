function [C, J_i, R_e, R_n, tilde_R_e, tilde_R_n] = gen_mat(N_r, N_t, L, H_RSI, a, b, n_l, tgt_index, w)

% C(\theta_l)
C=zeros(N_r, N_t, L);
for l=1:L
    C(:,:,l)=a(:,l)*b(:,l).';
end

% J_i
J_i = @(i) [zeros(N_t, N_t*(i-1)), eye(N_t), zeros(N_t, N_t*(N_t-i))];

% R
R=zeros(N_t,N_t);
for i=1:N_t
    R=R+J_i(i)*(w*w')*J_i(i)';
end

%% R_e & R_n & R'_e & R'_n

% R_e
R_e=zeros(N_r,N_r);
for l=1:L
    R_e=R_e+n_l(l)*(C(:,:,l)*R*C(:,:,l)');
end

% R_n
R_n=H_RSI*R*H_RSI';

% R'_e
n_l_temp=n_l;
n_l_temp(tgt_index)=[];
C_temp=C;
C_temp(:,:,tgt_index)=[];
tilde_R_e=zeros(N_r*N_t,N_r*N_t);
for tilde_l=1:L-1
    tilde_R_e=tilde_R_e+n_l_temp(tilde_l)*kron(eye(N_t),C_temp(:,:,tilde_l))*(w*w')*kron(eye(N_t),C_temp(:,:,tilde_l))';
end

% R'_n
tilde_R_n=kron(eye(N_t),H_RSI)*(w*w')*kron(eye(N_t),H_RSI)';