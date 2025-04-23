function x_pred = disc_sim(LTI, u, x0)
n = size(LTI.A,1);
N = size(u,1);
x_pred = zeros(N,n);
x_pred(1,:) = x0';
for i = 1:(N-1)
    x_pred(i+1,:) = LTI.A * x_pred(i,:)' + LTI.B * u(i,:)';
end
end