function [t0, x0, u0] = shift_disc(LTI, Ts, t0, x0, u)
% shift the time instant to new instant
t0 = t0 + Ts;

% using first control input 
x0 = LTI.A * x0 + LTI.B * u(1,:)';

% guess for the next iteration
u0 = [u(2:end,:); u(end, :)];




end