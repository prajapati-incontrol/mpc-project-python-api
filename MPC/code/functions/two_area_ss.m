function [A, B, C, Dd, F] = two_area_ss(params)
D = params.D;
H = params.H;
Tt = params.Tt;
Tij = params.Tij;
R = params.R;
Tg = params.Tg;
beta = params.beta;

%     delf1                         delPg1          
% A = [-D(1)/(2 * H(1)),                 0,       1/(2 * H(1)),      -1/(2 * H(1)), zeros(1,8);
%      -1/(R(1) * Tg(1)),         -1/Tg(1),                   0,                 0, zeros(1,8);
%                      0,          1/Tt(1),            -1/Tt(1),                 0, zeros(1,8);
%       2*pi*(Tij(1,2) + Tij(1,3)),      0,                   0,                 0, -2*pi*Tij(1,2), zeros(1,3), -2*pi*Tij(1,3), zeros(1,3);
%       zeros(1,4), -D(2)/(2 * H(2)),                 0,       1/(2 * H(2)),      -1/(2 * H(2)), zeros(1,4);
%       zeros(1,4), -1/(R(2) * Tg(2)),         -1/Tg(2),                   0,                 0, zeros(1,4);
%       zeros(1,4),                 0,          1/Tt(2),            -1/Tt(2),                 0, zeros(1,4);
%       -2*pi*Tij(2,1), zeros(1,3), 2*pi*(Tij(2,1) + Tij(2,3)), zeros(1,3), -2*pi*(Tij(2,3)), zeros(1,3);
%       zeros(1,8), -D(3)/(2 * H(3)),                 0,       1/(2 * H(3)),      -1/(2 * H(3));
%       zeros(1,8), -1/(R(3) * Tg(3)),         -1/Tg(3),                   0,                 0;
%       zeros(1,8),                 0,          1/Tt(3),            -1/Tt(3),                 0;
%       -2*pi*Tij(3,1), zeros(1,3), -2*pi*Tij(3,2), zeros(1,3), 2*pi*(Tij(3,1) + Tij(3,2)), zeros(1,3)];
A = [-D(1)/(2 * H(1)),                 0,       1/(2 * H(1)),      -1/(2 * H(1)), zeros(1,3);
     -1/(R(1) * Tg(1)),         -1/Tg(1),                   0,                 0, zeros(1,3);
                     0,          1/Tt(1),            -1/Tt(1),                 0, zeros(1,3);
      2*pi*(Tij(1,2)),      0,                   0,                 0, -2*pi*Tij(1,2), zeros(1,2);
      zeros(1,4), -D(2)/(2 * H(2)),                 0,       1/(2 * H(2));
      zeros(1,4), -1/(R(2) * Tg(2)),         -1/Tg(2),                   0;
      zeros(1,4),                 0,          1/Tt(2),            -1/Tt(2)];

B = [zeros(1,2);
     1/(Tg(1)), 0;
     zeros(2,2);
     zeros(1,2);
     0, 1/(Tg(2));
     zeros(1,2)];

% C = [beta(1), 0, 0, 1, zeros(1,8);
%      zeros(1,4), beta(2), 0, 0, 1, zeros(1,4);
%      zeros(1,8), beta(3), 0, 0, 1];

C = [beta(1), 0, 0, 1, zeros(1,3);
     zeros(1,4), beta(2), 0, 0];
      

Dd = 0;

F = [-1/(2 * H(1)), 0;
     zeros(3,2);
     0, -1/(2 * H(2));
     zeros(2,2)];

end