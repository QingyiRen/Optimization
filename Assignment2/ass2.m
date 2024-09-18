clear all 
clc 

global func;
% parameters for link (u,d)
c = 60;
N_u_d = 3;
v_u_d = 60; 
l_veh = 7;
l_u_d = 4000;
beta_u_d_o1 = 0.4;
beta_u_d_o2 = 0.3;
beta_u_d_o3 = 0.3;
mu_u_d = 1800;

% parameters for link (o1,d)
N_o1_d = 3;
v_o1_d = 50; 
l_o1_d = 4000;
beta_o1_d_u = 0.3;
beta_o1_d_o2 = 0.4;
beta_o1_d_o3 = 0.3;
mu_o1_d = 1700;

E_1 = 8+6;
E_2 = 7+0;
E_3 = 9+3;
    
C_u_d = (N_u_d*l_u_d)/l_veh;
C_o1_d = (N_o1_d*l_o1_d)/l_veh;

%------------- alpha enter (u,d) -------------%
alpha_enter_u_d(1:20) = 2300 + 10*E_1; 
alpha_enter_u_d(21:40) = 1800 + 10*E_2; 
alpha_enter_u_d(41:c) = 2100 + 10*E_3; 

alpha_enter_o1_d(1:c) = 1800+10*E_1;

%------------- Cd,u-------------%
C_d_u(1:40) = 20-E_3/2; 
C_d_u(41:c) = 20+E_3/2; 

%------------- Cd,o2-------------%
C_d_o2(1:20) = 10+E_2/2;
C_d_o2(21:35) = 10+E_2/2-2*(21:35)/(E_1+1);
C_d_o2(36:45) = 20+E_2/2;
C_d_o2(46:c) = 20+E_2/2-2*(46:c)/(E_1+1);
   
%------------- Cd,o1-------------%
C_d_o1(1:c)=C_d_o2(1:c)-2*(1:c)/(E_1+1);
   
%------------- Cd,o3-------------%
C_d_o3(1:20) = 10-E_3/2;
C_d_o3(21:c) = 10+E_3/2;

for k = 1:c
    if C_d_o3(k) < 0
        C_d_o3(k) = 0;
    end

    if C_d_o2(k) < 0
        C_d_o2(k) = 0;
    end

    if C_d_o1(k) < 0
        C_d_o1(k) = 0;
    end

    if C_d_u(k) < 0
        C_d_u(k) = 0;
    end
end

x0= [zeros(1,c) zeros(1,c) zeros(1,c) zeros(1,c)];
lb = ones(1,4*c)*15/3600;
ub = ones(1,4*c)*45/3600;

Aeq = [ones(1,c), -ones(1,c), zeros(1,c), zeros(1,c); zeros(1,c), zeros(1,c), ones(1,c), -ones(1,c); zeros(1,c), ones(1,c), ones(1,c) zeros(1,c)];
beq = [0; 0; c/3600];

options= optimoptions("fmincon","Algorithm","sqp");
[x,fval,exitflag,output] = fmincon(@(x)myfun(x,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d),x0,[],[],Aeq,beq,lb,ub,[],options)

function f = myfun(x,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d)

    global func;

    g_1 = x(1:60);
    g_2 = x(61:120);
    g_3 = x(121:180);
    g_4 = x(181:end);

    % parameters for link (u,d)
    c = 60;
    N_u_d = 3;
    v_u_d = 60; 
    l_veh = 7;
    l_u_d = 4000;
    beta_u_d_o1 = 0.4;
    beta_u_d_o2 = 0.3;
    beta_u_d_o3 = 0.3;
    mu_u_d = 1800;
    
    % parameters for link (o1,d)
    N_o1_d = 3;
    v_o1_d = 50; 
    l_o1_d = 4000;
    beta_o1_d_u = 0.3;
    beta_o1_d_o2 = 0.4;
    beta_o1_d_o3 = 0.3;
    mu_o1_d = 1700;
    

    for k = 1:c

        if k==1
          
            n_u_d(k)=0;
            n_o1_d(k)=0;
            q_u_d_o1(k) = 0;
            q_u_d_o2(k) = 0;
            q_u_d_o3(k) = 0;
    
            q_o1_d_u(k) = 0;
            q_o1_d_o2(k) = 0;
            q_o1_d_o3(k) = 0;


            tau_u_d(k) = floor((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh/(N_u_d*v_u_d*c/3600));
            tau_o1_d(k) = floor((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh/(N_o1_d*v_o1_d*c/3600));            
            
            gamma_u_d(k) = rem((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh, N_u_d*v_u_d*c/3600)/(N_u_d*v_u_d*c/3600);
            gamma_o1_d(k) = rem((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh, N_o1_d*v_o1_d*c/3600)/(N_u_d*v_u_d*c/3600);
        
            alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k)*(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(k)*(k-tau_u_d(k)-1);
            alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k)*(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(k)*(k-tau_o1_d(k)-1);
        
            alpha_leave_u_d_o1(k) = min([beta_u_d_o1*mu_u_d*g_1(k)/c/3600, q_u_d_o1(k)/c/3600 + beta_u_d_o1*alpha_arrive_u_d(k), C_d_o1(k)/c/3600]);
            alpha_leave_u_d_o2(k) = min([beta_u_d_o2*mu_u_d*g_2(k)/c/3600, q_u_d_o2(k)/c/3600 + beta_u_d_o2*alpha_arrive_u_d(k), C_d_o2(k)/c/3600]);
            alpha_leave_u_d_o3(k) = min([beta_u_d_o3*mu_u_d*(c/3600)/c/3600, q_u_d_o3(k)/c/3600 + beta_u_d_o3*alpha_arrive_u_d(k), C_d_o3(k)/c/3600]);
        
            alpha_leave_o1_d_u(k) = min([beta_o1_d_u*mu_o1_d*(c/3600)/(c/3600), q_o1_d_u(k)/c/3600 + beta_o1_d_u*alpha_arrive_o1_d(k), C_d_u(k)/c/3600]);
            alpha_leave_o1_d_o2(k) = min([beta_o1_d_o2*mu_o1_d*g_3(k)/c/3600, q_o1_d_o2(k)/c/3600 + beta_o1_d_o2*alpha_arrive_o1_d(k), C_d_o2(k)/c/3600]);
            alpha_leave_o1_d_o3(k) = min([beta_o1_d_o3*mu_o1_d*g_4(k)/c/3600, q_o1_d_o3(k)/c/3600 + beta_o1_d_o3*alpha_arrive_o1_d(k), C_d_o3(k)/c/3600]);

        else
                    
            n_u_d(k)=n_u_d(k-1)+(alpha_enter_u_d(k-1)-(alpha_leave_u_d_o1(k-1)+alpha_leave_u_d_o2(k-1)+alpha_leave_u_d_o3(k-1)))*c/3600;
            n_o1_d(k)=n_o1_d(k-1)+(alpha_enter_o1_d(k-1)-(alpha_leave_o1_d_u(k-1)+alpha_leave_o1_d_o2(k-1)+alpha_leave_o1_d_o3(k-1)))*c/3600;

            q_u_d_o1(k) = q_u_d_o1(k-1)+(beta_u_d_o1*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o1(k-1))*c/3600;
            q_u_d_o2(k) = q_u_d_o2(k-1)+(beta_u_d_o2*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o2(k-1))*c/3600;
            q_u_d_o3(k) = q_u_d_o3(k-1)+(beta_u_d_o3*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o3(k-1))*c/3600;
    
            q_o1_d_u(k) = q_o1_d_u(k-1)+(beta_o1_d_u*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_u(k-1))*c/3600;
            q_o1_d_o2(k) = q_o1_d_o2(k-1)+(beta_o1_d_o2*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_o2(k-1))*c/3600;
            q_o1_d_o3(k) = q_o1_d_o3(k-1)+(beta_o1_d_o3*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_o3(k-1))*c/3600;      

            tau_u_d(k) = floor((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh/(N_u_d*v_u_d*c/3600));
            tau_o1_d(k) = floor((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh/(N_o1_d*v_o1_d*c/3600));
            
            gamma_u_d(k) = rem((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh, N_u_d*v_u_d*c/3600)/(N_u_d*v_u_d*c/3600);
            gamma_o1_d(k) = rem((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh, N_o1_d*v_o1_d*c/3600)/(N_u_d*v_u_d*c/3600);
        
            alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k)*(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(k)*(k-tau_u_d(k)-1);
            alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k)*(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(k)*(k-tau_o1_d(k)-1);
        
            alpha_leave_u_d_o1(k) = min([beta_u_d_o1*mu_u_d*g_1(k-1)/(c/3600), q_u_d_o1(k)/(c/3600) + beta_u_d_o1*alpha_arrive_u_d(k), C_d_o1(k)/(c/3600)]);
            alpha_leave_u_d_o2(k) = min([beta_u_d_o2*mu_u_d*g_2(k-1)/(c/3600), q_u_d_o2(k)/(c/3600) + beta_u_d_o2*alpha_arrive_u_d(k), C_d_o2(k)/(c/3600)]);
            alpha_leave_u_d_o3(k) = min([beta_u_d_o3*mu_u_d*(c/3600)/(c/3600), q_u_d_o3(k)/(c/3600) + beta_u_d_o3*alpha_arrive_u_d(k), C_d_o3(k)/(c/3600)]);
        
            alpha_leave_o1_d_u(k) = min([beta_o1_d_u*mu_o1_d*(c/3600)/c/3600, q_o1_d_u(k)/c/3600 + beta_o1_d_u*alpha_arrive_o1_d(k), C_d_u(k)/c/3600]);
            alpha_leave_o1_d_o2(k) = min([beta_o1_d_o2*mu_o1_d*g_3(k-1)/c/3600, q_o1_d_o2(k)/c/3600 + beta_o1_d_o2*alpha_arrive_o1_d(k), C_d_o2(k)/c/3600]);
            alpha_leave_o1_d_o3(k) = min([beta_o1_d_o3*mu_o1_d*g_4(k-1)/c/3600, q_o1_d_o3(k)/c/3600 + beta_o1_d_o3*alpha_arrive_o1_d(k), C_d_o3(k)/c/3600]);
         
        end
        func(k) = (n_u_d(k)+n_o1_d(k))*c/3600;
    end
    f = sum(func);
end 