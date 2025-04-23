%% FREQUENCY REGULATOR OF A TWO AREA SYSTEM WITH OUTPUT FEEDBACK AND DISTURBANCE REJECTION
clc, clear, close all

import casadi.*
addpath('functions\')


% defining the parameters
params.D = [0.015, 0.016];
params.H = [0.1667, 0.2017]./(2);
params.R = [3, 2.73];
params.Tg = [0.08, 0.06];
params.Tt = [0.4, 0.44];
params.beta = [0.3483, 0.3827];
params.Tij = [0, 0.2;
    0.2, 0];

Ts = 0.1;
w_on_states = 10;
w_on_inputs = 1;

% return a discretized augmented system with sample time Ts
% and weight matrices for the extended system and state dimension
[sysd, LTI, LTIe, weighte.Q, weighte.R, weighte.P, dim] = aug_sys_disc(params, Ts, w_on_states, w_on_inputs);


% defining the system dimension
dime.nx = size(LTIe.A,1);
dime.nu = size(LTIe.B,2);
dime.ny = size(LTIe.C,1);
% horizon
dime.N = 5;

% Prediction model
[T_pred, S_pred] = predmodgenX(LTIe, dime);
predmode.T = T_pred;
predmode.S = S_pred;

[He, he] = costgenO(predmode, weighte, dime);

% Observer gain 
L = place(LTIe.A', LTIe.C', [0.41;0.52;0.43;0.64;0.42;0.51;0.57;0.49;0.59])';

% Check if A - LC is stable
if max(eig(LTIe.A - L*LTIe.C)) < 1
    disp("A-LC is stable. Error dynamics are asymptotically stable!");
end

% Obtaining the terminal set
grc_con = 0.2; % generation rate constraint
tie_con = 0.03; % tie-line constraint
other_x = 1;    % other states constraints
d_dist = 0.5;
xlb = [-other_x;    % del_f1 loose 
       -grc_con;    % del_pg1 GRC
       -other_x;    % del_pm1 loose   
       -tie_con;    % del_ptie12 tie-line 
       -other_x;    % del_f2 loose 
       -grc_con;    % del_pg2 GRC
       -other_x];   % del_pm2 loose

xub = -1.*xlb;

u_con = 0.25;
ulb = [-u_con;    % del_pc1
       -u_con];   % del_pc2

uub = -1*ulb;

[Xf_H, Xf_h] = calcLQRXf(sysd, xlb, xub, ulb, uub, weighte.Q(1:dim.nx, 1:dim.nx), weighte.R(1:dime.nu, 1:dime.nu));

% Symbolic inputs
U = SX.sym('U', dime.N*dime.nu, 1);

% Symbolic parameters
% 7 + 2 = 9 xhat(k) or x(k)??
% 7 + 2 = 9 x_ref from OTS
% 2 u_ref from OTS
P_param = SX.sym('P', dime.nx*2 + dime.nu);

% State constraints as a function of symbolic input
g = [];

for j = 1:dim.nx
    for i = 1:(dime.N-1)
    g = [g;
         T_pred(i*dime.nx + j, :)*P_param(1:dime.nx) + S_pred(i*dime.nx + j, :)*U];
    end
end

% Terminal constraint Xf: H*x(N) <= h
constraints.Ae = S_pred((end - dime.nx + 1):end - dime.nu, :);
constraints.be = T_pred((end - dime.nx + 1):end - dime.nu, :);

g = [g;
     Xf_H * (constraints.Ae*U + constraints.be*P_param(1:dime.nx)) - Xf_h];



% Declaring the objective function in terms of symbolic inputs and
% parameters
obj = 0.5 * U' * He * U + (he*P_param)' * U;

OPT_variable = U;
qp_prob = struct('f', obj, 'x', OPT_variable, 'g', g, 'p', P_param);
solver = qpsol('solver', 'qpoases', qp_prob);

args = struct;
args.lbg = zeros(size(g));
args.ubg = zeros(size(g));

% Loose constraints
          % x1       % x3                % x5                   % x7
args.lbg([1:(dime.N-1), (2*(dime.N-1) +1):3*(dime.N-1), (4*(dime.N-1)+1):5*(dime.N-1), (6*(dime.N-1)+1):7*(dime.N-1)],:) = -other_x;
args.ubg([1:(dime.N-1), (2*(dime.N-1) +1):3*(dime.N-1), (4*(dime.N-1)+1):5*(dime.N-1), (6*(dime.N-1)+1):7*(dime.N-1)],:) = other_x;

% Generation-rate constraints
          % x2       % x6
args.lbg([((dime.N-1)+1):2*(dime.N-1), (5*(dime.N-1)+1):6*(dime.N-1)],:) = -grc_con;
args.ubg([((dime.N-1)+1):2*(dime.N-1), (5*(dime.N-1)+1):6*(dime.N-1)],:) = grc_con;

% Tie-line constraints
args.lbg((3*(dime.N-1)+1):4*(dime.N-1), :) = -tie_con;
args.ubg((3*(dime.N-1)+1):4*(dime.N-1), :) = tie_con;
%
% Xf inequalities
% H x(N) - h <= 0
args.lbg((7*(dime.N-1)+1):end,:) = -inf;
args.ubg((7*(dime.N-1)+1):end,:) = 0;

% Constraints on Input: Decision variables for optimization
args.lbx = -u_con;
args.ubx = u_con;

% THE SIMULATION SHOULD START FROM HERE
% -------------------------------------

% Simulation Time
sim_tim = 20;
T = sim_tim/Ts + 1;

% initialize the reference values
xref = zeros(dim.nx, T);  % state without disturbance
uref = zeros(dime.nu, T); % inputs

% MPC Control Law
% For all iterations capturing first input only
u_cl = zeros(T, dime.nu);

% input stack for full horizon used at every iteration
u_ncl = zeros(dime.N, dime.nu, T);

% extended state-estimate history
x_hat_d = zeros(dime.nx, T);
% extended output-estimate history
y_hat_d = zeros(dime.ny, T);

% extended state history
x_d = zeros(dime.nx, T);
y_d = zeros(dime.ny, T);

% Setting up...
% defining initial conditions for
% some bad initial estimates
x_hat_d(:,1) = [ 0.01         % d_f1
                 0.01         % d_pg1
                 0.01         % d_m1
                 0.01         % d_ptie1
                -0.01         % d_f2
                 0.01         % d_pg2
                 0.01         % d_pm2
                 0.01         % d_pl1
                 0.01];       % d_pl2
y_hat_d(:,1) = LTIe.C * x_hat_d(:,1);

% X_N initial state estimates
x_d(:,1) = [-0.03         % d_f1
             0.00         % d_pg1
             0.00         % d_m1
             0.03         % d_ptie1
            -0.03         % d_f2
             0.00         % d_pg2
             0.00         % d_pm2
             0.001         % d_pl1
             0.00];       % d_pl2
y_d(:,1) = LTIe.C * x_d(:,1);

% Constant Output reference
LTIe.yref = [0.00;      % ACE_1
             0.00];     % ACE_2

% At every iteration
% [xref, uref] = OTS(LTI, Ax=b, d_hat(from last iter))
% u* = MPC_Regulator(LTI, solver, xhat(from last iter), x_ref, u_ref)
% x(k+1) = LTIe.A x(k) + LTEe.B u*(k)
% y(k+1) = LTIe.C x(k+1) 
% xhat(k+1) = LTIe.A xhat(k) + LTIe.B u*(k)
% yhat(k+1) = LTIe.C xhat(k+1)

k = 1;
while (norm(y_d(:,k) - LTIe.yref, 2) > 1e-4 && k <= T)

    % OTS 
    % extract d_hat from last iteration
    d_hat = x_d((dim.nx +1):end,k); % last two rows
    eqconstraints = eqconstraintsgen(LTIe, dim, dime, d_hat);
    [xref(:,k), uref(:,k)] = calcOTS(dim, dime, eqconstraints);

    % MPC
    % Setup the parameters P_params for MPC
    args.p = [x_d(:,k); xref(:,k); d_hat; uref(:,k)];

    % use the initial guess of last iteration for MPC
    u0t = u_ncl(:,:,k)';
    % since matlab looks at things row-wise
    args.x0 = u0t(:);

    % calling the solver to find the optimal input sequence
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, ...
        'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);

    u_ncl(:, :, k+1) = reshape(full(sol.x)', dime.nu, dime.N)';
    
    % extract the first input and store it 
    u_cl(k,:) = u_ncl(1,:, k+1) + uref(:,k)';

    % Use this first input to calculate the state and output
    x_d(:, k+1) = LTIe.A * x_d(:,k) + LTIe.B * u_cl(k,:)';
    y_d(:, k+1) = LTIe.C * x_d(:,k+1);

    % Use the first input to calculate the state and output estimates
    x_hat_d(:, k+1) = LTIe.A * x_hat_d(:,k) + LTIe.B * u_cl(k,:)' + L * (y_d(:,k) - LTIe.C * x_hat_d(:,k));
    y_hat_d(:, k+1) = LTIe.C * x_hat_d(:,k+1);
    k = k + 1;
end

%%
stairs(x_hat_d(1,:))
hold on
stairs(x_d(1,:))


