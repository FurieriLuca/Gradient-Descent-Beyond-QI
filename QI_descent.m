%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION :
%       Numerical example of Section III.B attached to the paper: 
%       
%       "First Order Methods For Globally Optimal Distributed Controllers Beyond Quadratic Invariance"
%        by Luca Furieri (furieril@control.ethz.ch) and Maryam Kamgarpour
%        (mkamgar@control.ee.ethz.ch)
%%%%%%%%%%%%%%%%%%


clear all;
clc;


%%Linear system definition
A=[1.6 0 0 0 0;0.5 1.6 0 0 0;2.5 2.5 -1.4 0 0;-2 1 -2 0.1 0;0 2 0 -0.5 1.1];
n=size(A,1);
B=eye(n);
C=eye(n);
m=size(B,2);
p=size(C,1);

N=3; 
mu0=[1;-1;2;-3;3];

Sigmaw=eye(n);
Sigmav=eye(p);



%%Stacked operator definition. "_b" stands for bold notation in the paper 
A_b=kron(eye(N+1),A);
B_b=[kron(eye(N),B);zeros(n,m*N)];
C_b=[kron(eye(N),C) zeros(p*N,n)];
Z=zeros(n*(N+1),n*(N+1));
for(i=1:N+1)
        for(j=1:N+1)
                if(i==j+1)
                        Z([(i-1)*n+1:i*n],[(j-1)*n+1:j*n])=eye(n);
                end
        end
end


P11=round(inv(eye(n*(N+1))-Z*A_b),3);
P12=round(inv(eye(n*(N+1))-Z*A_b)*Z*B_b,3);

Sigmaw_b=kron(eye(N+1),Sigmaw);
Sigmav_b=kron(eye(N),Sigmav);
mu_w=[mu0;zeros(N*n,1)];



%%Cost function parameters
M=eye(n);
R=eye(m);
M_b=blkdiag(eye(n),kron(eye(N),M));
R_b=kron(eye(N),R);



%%Information structure definition
S = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 1];
struct_small=ones(N,N);
struct_small=tril(struct_small);
struct=kron(struct_small,S); %%%% this struct is \textbf{S} in the paper. QI can be verified by checking that the structure of struct*C_b*P12*struct is included in struct itself.
cardinality=sum(sum(struct));
K=sym('K',[m*(N),p*(N)]);
assume(K,'real');
K=K.*struct;


%COST FUNCTION DEFINITION as per equation (6) of the paper

fprintf('Encoding the cost function...\n')
create_cost;



%%%stacking the non-zero decision variables into a single vector
vec_K=sym(zeros(cardinality,1));
count=0;
positions=zeros(cardinality,2);
for(i=1:m*N)
        for(j=1:p*N)
                if(struct(i,j)==1)
                        count=count+1;
                        vec_K(count)=K(i,j);
                        positions(count,:)=[i,j];
                end
        end
end


%%Compute the gradient symbolically
fprintf('Encoding the gradient...\n')
grad=simplify(gradient(cost,vec_K));


%% Encode the gradient and the cost as efficient matlab functions
fprintf('Encoding eval_cost...\n')
eval_cost=matlabFunction(cost,'Vars',{vec_K});

fprintf('Encoding eval_grad...\n')
eval_gradient=matlabFunction(grad,'Vars',{vec_K});



%%% GRADIENT DESCENT
parameters=10*(rand(cardinality,1)-rand(cardinality,1)); %% Initial controller is chosen at random
parameters_initial=parameters;
fprintf('Starting the gradient descent...\n')
exit=0;
t=0;
while(1)
        t=t+1;
        
        gradient_step=eval_gradient(parameters);
        
        eta=stepsize_bisection(eval_cost,eval_gradient,parameters); %%%Using the bisection algorithm of Proposition 5.5, "Nonlinear Optimization", F.J. Aragon and others, 2019
        
        parameters=parameters-eta*gradient_step;
        
        if(mod(t,100)==1)
                new_cost=eval_cost(parameters)
        end
        if(max(abs(eval_gradient(parameters)))<5e-5)
                exit=1;
        end
        if(exit==1)
                break;
        end
        
end

fprintf('DONE!\n')
final_cost=eval_cost(parameters)
initial_cost=eval_cost(parameters_initial)


%%rebuild K
Ksol=zeros(m*N,p*N);
for(i=1:cardinality)
        coord_x=positions(i,1);
        coord_y=positions(i,2);
        Ksol(coord_x,coord_y)=parameters(i);
end


