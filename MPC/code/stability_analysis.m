%% FREQUENCY REGULATOR OF A TWO AREA SYSTEM
clc, clear, close all

import casadi.*
addpath('functions\')

params.D    = [0.015, 0.016];
params.H    = [0.1667, 0.2017]./(2);
params.R    = [3, 2.73];
params.Tg   = [0.08, 0.06];
params.Tt   = [0.4, 0.44];
params.beta = [0.3483, 0.3827];
params.Tij  = [0.0, 0.2;
              0.2, 0.0];
Ts = 0.1;
% Definition of the LTI system

[A, B, C, D, F] = two_area_ss(params);

if rank(ctrb(A,B)) == size(A,1)
    disp("The continuous system is controllable!");
end

if rank(ctrb(A',C')) == size(A,1)
    disp("The continuous system is observable!");
end

sys = ss(A,B,C,D);

% Discretizing the system
sysd = c2d(sys,Ts);
LTI = struct;
LTI.A = sysd.A;
LTI.B = sysd.B;
LTI.C = sysd.C;
LTI.D = sysd.D;
LTI.F = F;
LTI.x0 = [-0.03         % d_f1
           0.01         % d_pg1
           0.01         % d_m1
           0.03         % d_ptie1
          -0.01         % d_f2
           0.01         % d_pg2
          -0.01];       % dpm2

LTI.xref = zeros(size(LTI.x0,1),1);

% Defining the system dimension
dim.nx = size(LTI.A, 1); % state dimension
dim.ny = size(LTI.C, 1); % output dimension
dim.nu = size(LTI.B, 2); % input dimension
dim.N = 10;

% Defining the weights 
w_on_states = 10;
weights.Q = w_on_states.*eye(dim.nx);

w_on_inputs = 1;
weights.R = w_on_inputs.*eye(dim.nu);

% DARE solution
[Pdare, Kdare, ~] = idare(LTI.A, LTI.B, weights.Q, weights.R);

Kdare = -Kdare;

% Prediction model
[P_pred, S_pred] = predmodgenX(LTI, dim);

[H, h, const] = newcostgen(P_pred, S_pred, weights.Q, weights.R, dim, Pdare);


% Symbolic inputs
U = SX.sym('U', dim.N*dim.nu, 1);

% first 7 initial conditions for x
% last seven states as reference setpoint for regulation
P_param = SX.sym('P',dim.nx*2,1);


% State constraints as a function of symbolic input
g = [];

for j = 1:dim.nx
    for i = 1:(dim.N-1) % constraints on x(0):x(N-1)
    g = [g;
         P_pred(i*dim.nx + j, :)*P_param(1:dim.nx) + S_pred(i*dim.nx + j, :)*U];
    end
end

% Obtaining the terminal set x(N)
grc_con = 0.2; % generation rate constraint
tie_con = 0.03; % tie-line constraint
other_x = 0.3;    % other states constraints
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

[Xf_H, Xf_h] = calcLQRXf(sysd, xlb, xub, ulb, uub, weights.Q, weights.R);

% Generating the Polytope for Xf
P = polytope(Xf_H, Xf_h);                            
flag = isbounded(P);

% Check if the terminal set in bounded
if flag == true
    disp('The terminal set is bounded.')
    pause(1);
end

% Check if the Xf contain origin
is_origin = (Xf_H * zeros(dim.nx,1) - Xf_h) < 0;
if sum(is_origin) == size(Xf_h,1)
    disp('The origin is in the interior of the terminal set.')
    pause(1);
end

% Check if Xf contains initial condition
is_init = (Xf_H * LTI.x0 - Xf_h) < 0;
if sum(is_init) == size(Xf_h, 1)
    disp("The initial state x0 is in the interior of the terminal set.")
    pause(1);
else 
    disp("The initial state x0 is not in the interior of the terminal set.")
    pause(1);
end
% Terminal constraint Xf:  H*x(N) <= h
constraints.Ae = S_pred(end - dim.nx + 1:end,:);
constraints.be = P_pred(end - dim.nx + 1:end,:);

g = [];
g = P_pred(dim.nx +1:(end-dim.nx), :)*P_param(1:dim.nx) + ...
    S_pred(dim.nx +1:(end-dim.nx), :)*U;

% constraints on terminal x(N)
% Terminal constraint Xf:  H*x(N) <= h
constraints.Ae = S_pred(end - dim.nx + 1:end,:);
constraints.be = P_pred(end - dim.nx + 1:end,:);

g = [g; 
    Xf_H * (constraints.Ae*U  + constraints.be*P_param(1:dim.nx)) - Xf_h]; 

% Objective function
obj = 0.5 * U' * H * U + (h*P_param(1:dim.nx))' * U; %+ P_param(1:dim.nx)' * const * P_param(1:dim.nx);

options = struct;

OPT_variables = U;
qp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P_param);

% qpsol
solver = qpsol('solver', 'qpoases', qp_prob, options);
%
args = struct;
args.lbg = zeros(size(g));
args.ubg = zeros(size(g));

% constraining x(1) to x(N-1)
xlowlimit = kron(ones(dim.N-1,1),xlb);

args.lbg(1:(dim.N-1)*dim.nx,:) = xlowlimit;
args.ubg(1:(dim.N-1)*dim.nx,:) = -xlowlimit;

% constraining x(N)
args.lbg(((dim.N-1)*dim.nx + 1):end,:) = -inf;
args.ubg(((dim.N-1)*dim.nx + 1):end,:) = 0;

% Constraints on Input: Decision variables for optimization
args.lbx = -u_con;
args.ubx = u_con;


%
% THE SIMULATION SHOULD START FROM HERE
% ------------------------------------

t0 = 0;

x0 = LTI.x0;
xref = zeros(dim.nx,1);


xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(dim.N,dim.nu); % two control inputs
sim_tim = 50;

% start MPC
mpciter = 0; % counter for the loop
xx1 = [];
u_cl = [];


while ( norm( (x0 - xref), 2 ) > 1e-4 && mpciter < sim_tim / Ts )
% while mpciter < sim_tim / Ts 
    args.p = [x0; xref]; % set the values of the parameters

    % initial value of the decision variables
    u0t = u0';
    args.x0 = u0t(:);

    % calling the solver to find the optimal input sequence
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, ...
        'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);

    % reshaping u to get back a matrix
    u = reshape(full(sol.x)', dim.nu,dim.N)'; % u = rows: N, columns: inputs

    % STEP HERE
    x_pred = disc_sim(LTI, u, x0);

    % 3D PAGE HERE
    xx1(:,1:dim.nx, mpciter + 1) = x_pred;
    %
    u_cl = [u_cl; u(1, :)];

    t(mpciter + 1) = t0;
    
    [t0, x0, u0] = shift_disc(LTI, Ts, t0, x0, u);

    xx(:, mpciter + 2) = x0;
    mpciter = mpciter + 1;
end

% CLF Decrease

% u_mpc
u_mpc = u_cl';

% stage-cost 
l_xu = zeros(1, size(u_mpc,2));

% terminal-cost
vf_xu = zeros(1, size(xx,2));

for count = 1:size(u_mpc,2)
    l_xu(1, count) = 0.5 * xx(:,count)' * weights.Q * xx(:,count) + ...
        0.5 * u_mpc(:,count)' * weights.R * u_mpc(:, count);
    vf_xu(1, count) = 0.5 * xx(:,count)' * Pdare * xx(:, count);
end

% for last element
vf_xu(1, size(xx,2)) = 0.5 * xx(:,size(xx,2))' * Pdare * xx(:, size(xx,2));

% vf_diff should be the same size as l_xu
vf_diff = diff(vf_xu);

% vf_diff + l_xu <= 0 implies stability
clf_diff = vf_diff + l_xu;

figure 
stairs(t, clf_diff,'linewidth',2);
legend('$$V_{f}(x^+) - V_{f}(x) + l(x,u)$$','Interpreter','LaTex','Fontsize',12);
xlabel('Time [s]')
grid on
fontname(gcf,"Garamond")
fontsize(gcf,20,"pixels")
set(gcf,'Position',[150 50 500 500])














