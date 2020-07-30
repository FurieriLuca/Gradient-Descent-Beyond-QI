%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION :
%       Numerical example of Section III.B attached to the paper: 
%       
%       "First Order Methods For Globally Optimal Distributed Controllers Beyond Quadratic Invariance"
%        by Luca Furieri (furieril@control.ethz.ch) and Maryam Kamgarpour
%        (mkamgar@control.ee.ethz.ch)

% This file validates the result of "QI_descent.m" by solving the
% corresponding convex program in the Youla parameter.
%%%%%%%%%%%%%%%%%%



clear all;
clc;


A=[1.6 0 0 0 0;0.5 1.6 0 0 0;2.5 2.5 -1.4 0 0;-2 1 -2 0.1 0;0 2 0 -0.5 1.1];
n=size(A,1);
B=eye(n);
C=eye(n);

m=size(B,2);
p=size(C,1);

N=3; 
x0=1*[1;-1;2;-3;3];

A_tilde=kron(eye(N+1),A);
B_tilde=[kron(eye(N),B);zeros(n,m*N)];
C_tilde=[kron(eye(N),C) zeros(p*N,n)];
Z=zeros(n*(N+1),n*(N+1));
for(i=1:N+1)
    for(j=1:N+1)
        if(i==j+1)
            Z([(i-1)*n+1:i*n],[(j-1)*n+1:j*n])=eye(n);
        end
    end
end


P11=round(inv(eye(n*(N+1))-Z*A_tilde),3);
P12=round(inv(eye(n*(N+1))-Z*A_tilde)*Z*B_tilde,3);

M=eye(n);
R=eye(m);
Sigmaw=eye(n);
Sigmav=eye(p);

M_tilde=blkdiag(eye(n),kron(eye(N),M));
R_tilde=kron(eye(N),R);
Sigmaw_tilde=kron(eye(N+1),Sigmaw);
Sigmav_tilde=kron(eye(N),Sigmav);
mu_w=[x0;zeros(N*n,1)];



S = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 1];
struct_small=ones(N,N);
struct_small=tril(struct_small);
struct=kron(struct_small,S);
cardinality=sum(sum(struct));

Q=sdpvar(m*N,p*N,'full');
Q=Q.*struct;


%COST in Q
w_x_cost=M_tilde^0.5*(eye(size(P12,1))+P12*Q*C_tilde)*P11*Sigmaw_tilde^0.5;
w_u_cost=R_tilde^0.5*Q*C_tilde*P11*Sigmaw_tilde^0.5;
v_x_cost=M_tilde^0.5*P12*Q*Sigmav_tilde^0.5;
v_u_cost=R_tilde^0.5*Q*Sigmav_tilde^0.5;
x0_x_cost=M_tilde^0.5*(eye(size(P12,1))+P12*Q*C_tilde)*P11*mu_w;
x0_u_cost=R_tilde^0.5*Q*C_tilde*P11*mu_w;


cost=trace(w_x_cost'*w_x_cost)+trace(w_u_cost'*w_u_cost)+trace(v_x_cost'*v_x_cost)+trace(v_u_cost'*v_u_cost)+x0_x_cost'*x0_x_cost+x0_u_cost'*x0_u_cost;

ops=sdpsettings('solver','mosek'); %also works with quadprog
sol=optimize([], cost, ops) 

optimal_value=value(cost)
