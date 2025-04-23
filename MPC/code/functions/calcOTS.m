function [xref, uref] = calcOTS(dim, dime, eqcon) 

Q = 10*eye(dim.nx);
R = eye(dime.nu);

H = blkdiag(Q, R);

h = zeros(dim.nx + dime.nu,1);

options1 = optimoptions(@quadprog);
options1.OptimalityTolerance = 1e-20;
options1.ConstraintTolerance = 1.0000e-15;
options1.Display = 'off';
xur = quadprog(H, h, [], [], eqcon.A, eqcon.b, [],[],[],options1);
xref = xur(1:dim.nx);
uref = xur(dim.nx+1:end);
end
