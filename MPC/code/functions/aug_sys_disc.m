function [sys_d, LTI, LTIe, Qe, Re, Pe, dim] = aug_sys_disc(params, Ts, w_on_states, w_on_inputs)
[A, B, C, D,Bd] = two_area_ss(params);

dim.nx = size(A,1);
dim.nu = size(B,2);
dim.nd = size(Bd,2);

Cd = [0.01 0.00;
      0.00 0.03];

sys_d = c2d(ss(A, B, C, D), Ts);
LTI.A = sys_d.A;
LTI.B = sys_d.B;
LTI.C = sys_d.C;
LTI.D = sys_d.D;

% Augmented state-space system for disturbance rejection
A_aug = [A, Bd;
         zeros(dim.nu,dim.nx), zeros(dim.nu)];
B_aug = [B; zeros(dim.nu)];
C_aug = [C, Cd];
D_aug = zeros(dim.nu);

% Check if the augmented system is observable
Obs_mat = [eye(dim.nx) - A, -Bd;
           C, Cd];

if rank(Obs_mat) == dim.nx + dim.nd
    disp('Continuous augmented system is observable!');
end

% Augmented system discretized 
sys_ce = ss(A_aug, B_aug, C_aug, D_aug);

% discretized augmented system
sys_de = c2d(sys_ce, Ts);

LTIe.A = sys_de.A;
LTIe.B = sys_de.B;
LTIe.C = sys_de.C;
LTIe.D = sys_de.D;

% Check if the discrete augmented system is observable
if rank(ctrb(LTIe.A', LTIe.C')) == dim.nx + dim.nd
disp("Discrete Augmented System is observable!")
end

% Defining weights
Q = w_on_states.*eye(dim.nx);
Qe = blkdiag(Q, zeros(dim.nd));

R = w_on_inputs.* eye(dim.nu);
Re = R;

P = idare(LTI.A,LTI.B,Q,R);
Pe = blkdiag(P, zeros(dim.nd));
end