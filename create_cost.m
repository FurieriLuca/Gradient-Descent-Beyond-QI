w_x_cost_K=M_b^0.5*inv(eye(size(P12,1))-P12*K*C_b)*P11*Sigmaw_b^0.5;
v_x_cost_K=M_b^0.5*P12*K*inv(eye(size(C_b,1))-C_b*P12*K)*Sigmav_b^0.5;
w_u_cost_K=R_b^0.5*K*inv(eye(size(C_b,1))-C_b*P12*K)*C_b*P11*Sigmaw_b^0.5;
v_u_cost_K=R_b^0.5*K*inv(eye(size(C_b,1))-C_b*P12*K)*Sigmav_b^0.5;
x0_x_cost_K=M_b^0.5*inv(eye(size(P12,1))-P12*K*C_b)*P11*mu_w;
x0_u_cost_K=R_b^0.5*K*inv(eye(size(C_b,1))-C_b*P12*K)*C_b*P11*mu_w;

cost=trace(w_x_cost_K'*w_x_cost_K)+trace(w_u_cost_K'*w_u_cost_K)+trace(v_x_cost_K'*v_x_cost_K)+trace(v_u_cost_K'*v_u_cost_K)+trace(x0_x_cost_K'*x0_x_cost_K)+trace(x0_u_cost_K'*x0_u_cost_K);
