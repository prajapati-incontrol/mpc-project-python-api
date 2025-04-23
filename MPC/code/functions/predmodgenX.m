function [P,S]=predmodgenX(LTI,dim)

%Prediction matrices generation
%This function computes the prediction matrices to be used in the
%optimization problem

%Prediction matrix from initial state
P=zeros(dim.nx*(dim.N),dim.nx);
for k=0:dim.N
    P(k*dim.nx+1:(k+1)*dim.nx,:) = LTI.A^k;
end

%Prediction matrix from input
S=zeros(dim.nx*(dim.N),dim.nu*(dim.N-1));
for k=1:dim.N
    for i=0:k-1
        S(k*dim.nx+1:(k+1)*dim.nx,i*dim.nu+1:(i+1)*dim.nu) = LTI.A^(k-1-i)*LTI.B;
    end
end
