%Build matrix

%STANDARDMATRISER
I = i_t*i_f;
J_1 = J_e + J_t/i_t^2 + J_f/I^2;
J_2 = J_w + m*r_w^2;
b_1 = b_t/i_t^2 + b_f/I^2;
A = [0 1/I -1 0;...
    -k/(I*J_1) -(b_1+c)/(I^2*J_1) c/(I*J_1) 0;...
    k/J_2 c/(I*J_2) -c/J_2 0;...
    0 0 0 -1/tau_e];
B = [0 0 0 1/tau_e]';
C = [0 1 0 0];
H = [0 0 -1/J_2 0]';

l = r_w*m*g*c_r1;
%Lös ricatti

W = eye(size(A));
V = 0.05;
S = zeros(size(B));
E = eye(size(A));
P_f = care(A,B,W,V,S,E); 

%Observatörmatriser
K_f = P_f*C'*inv(V);

%temp 
%K_f = zeros(size(K_f));

A_obs = A - K_f*C;
B_obs = [B K_f H];
C_obs = eye(4);
