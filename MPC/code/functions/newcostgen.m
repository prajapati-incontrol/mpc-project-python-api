function [H, h, const] = newcostgen(T, S, Q, R, dim, Pdare)
Qbar=blkdiag(kron(eye(dim.N),Q), Pdare);


H =S'*Qbar*S+kron(eye(dim.N),R);   
h =S'*Qbar*T;
const =T'*Qbar*T;
end