clear all 
close all

%% parameters
% parameters for link (u,d)
c = 60/3600;
N_u_d = 3;
v_u_d = 60; 
l_veh = 7/1000;
l_u_d = 4;
beta_u_d_o1 = 0.4;
beta_u_d_o2 = 0.3;
beta_u_d_o3 = 0.3;
mu_u_d = 1800;

% parameters for link (o1,d)
N_o1_d = 3;
v_o1_d = 50; 
l_o1_d = 4;
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
alpha_enter_u_d(1:21) = 2300 + 10*E_1; 
alpha_enter_u_d(22:41) = 1800 + 10*E_2; 
alpha_enter_u_d(42:61) = 2100 + 10*E_3; 

alpha_enter_o1_d(1:61) = 1800+10*E_1;

%------------- Cd,u-------------%
C_d_u(1:41) = 20-E_3/2; 
C_d_u(42:61) = 20+E_3/2; 

%------------- Cd,o2-------------%
C_d_o2(1:21) = 10+E_2/2;
C_d_o2(22:36) = 10+E_2/2-2*(21:35)/(E_1+1);
C_d_o2(37:46) = 20+E_2/2;
C_d_o2(47:61) = 20+E_3/2+2*(46:60)/(E_1+1);
   
%------------- Cd,o1-------------%
C_d_o1(1:61)=C_d_o2(1:61)-2*(0:60)/(E_1+1);
   
%------------- Cd,o3-------------%
C_d_o3(1:21) = 10-E_3/2;
C_d_o3(22:61) = 10+E_3/2;

for k = 1:61
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

%% question 3
lb = 15/3600;
ub = 45/3600;

rng default % For reproducibility
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',8e4,'TolFun',5e-2);

for k = 1:60

    n_u_d_(1) = 0;
    n_o1_d_(1) = 0;
 
    q_u_d_o1(1) = 0;
    q_u_d_o2(1) = 0;
    q_u_d_o3(1) = 0;
    
    q_o1_d_u(1) = 0;
    q_o1_d_o2(1) = 0;
    q_o1_d_o3(1) = 0;
    
    iteration(1) = 0;
    alpha_enter_u_d(k==1) = 0;
    alpha_enter_o1_d(k==1) = 0;
    value(1) = 0;
    x1 = 15/3600;
    green_light(1) = x1;

    tau_u_d(k) = floor((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh/(N_u_d*v_u_d*c));
    tau_o1_d(k) = floor((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh/(N_o1_d*v_o1_d*c));     

    gamma_u_d(k) = rem((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh, N_u_d*v_u_d*c)/(N_u_d*v_u_d*c);
    gamma_o1_d(k) = rem((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh, N_o1_d*v_o1_d*c)/(N_o1_d*v_o1_d*c);

    if k-tau_u_d(k)<=0
        alpha_arrive_u_d(k) = 0;
    else
        if k-tau_u_d(k)-1<=1
           alpha_arrive_u_d(k) = 0;
        elseif  k-tau_u_d(k)==2
           alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(2);
        else
           alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(k-tau_u_d(k)-1);
        end
    end
    
    if k-tau_o1_d(k)<=0
        alpha_arrive_o1_d(k) = 0;
    else
        if k-tau_o1_d(k)-1<1
           alpha_arrive_o1_d(k) = 0;
        elseif  k-tau_o1_d(k)==2
           alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(2);
        else 
           alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(k-tau_o1_d(k)-1);
        end
    end

    alpha_leave_u_d_o1 = @(x1)min([beta_u_d_o1*mu_u_d*x1/c, q_u_d_o1(k)/c + beta_u_d_o1*alpha_arrive_u_d(k), C_d_o1(k)/c]);
    alpha_leave_u_d_o2 =@(x1)min([beta_u_d_o2*mu_u_d*x1/c, q_u_d_o2(k)/c + beta_u_d_o2*alpha_arrive_u_d(k), C_d_o2(k)/c]);
    alpha_leave_u_d_o3 = min([beta_u_d_o3*mu_u_d, q_u_d_o3(k)/c + beta_u_d_o3*alpha_arrive_u_d(k), C_d_o3(k)/c]);

    alpha_leave_o1_d_u = min([beta_o1_d_u*mu_o1_d, q_o1_d_u(k)/c + beta_o1_d_u*alpha_arrive_o1_d(k), C_d_u(k)/c]);
    alpha_leave_o1_d_o2 = @(x1)min([beta_o1_d_o2*mu_o1_d*(c-x1)/c, q_o1_d_o2(k)/c + beta_o1_d_o2*alpha_arrive_o1_d(k), C_d_o2(k)/c]);
    alpha_leave_o1_d_o3 = @(x1)min([beta_o1_d_o3*mu_o1_d*(c-x1)/c, q_o1_d_o3(k)/c + beta_o1_d_o3*alpha_arrive_o1_d(k), C_d_o3(k)/c]);

    n_u_d=@(x1)(n_u_d_(k)+(alpha_enter_u_d(k)-(alpha_leave_u_d_o1(x1)+alpha_leave_u_d_o2(x1)+alpha_leave_u_d_o3))*c);
    n_o1_d=@(x1)(n_o1_d_(k)+(alpha_enter_o1_d(k)-(alpha_leave_o1_d_u+alpha_leave_o1_d_o2(x1)+alpha_leave_o1_d_o3(x1)))*c);

    func = @(x1)((n_u_d(x1)+n_o1_d(x1))*c);

    [x1,fval1,exitflag1,output1] = fmincon(func,15/3600,[],[],[],[],lb,ub,[],options);

    n_u_d_(k+1)=n_u_d(x1);
    n_o1_d_(k+1)=n_o1_d(x1);

    iteration(k+1)=(n_u_d_(k+1)+n_o1_d_(k+1))*c;
   
    q_u_d_o1(k+1) = q_u_d_o1(k)+(beta_u_d_o1*alpha_arrive_u_d(k)-alpha_leave_u_d_o1(x1))*c;
    q_u_d_o2(k+1) = q_u_d_o2(k)+(beta_u_d_o2*alpha_arrive_u_d(k)-alpha_leave_u_d_o2(x1))*c;
    q_u_d_o3(k+1) = q_u_d_o3(k)+(beta_u_d_o3*alpha_arrive_u_d(k)-alpha_leave_u_d_o3)*c;

    q_o1_d_u(k+1) = q_o1_d_u(k)+(beta_o1_d_u*alpha_arrive_o1_d(k)-alpha_leave_o1_d_u)*c;
    q_o1_d_o2(k+1) = q_o1_d_o2(k)+(beta_o1_d_o2*alpha_arrive_o1_d(k)-alpha_leave_o1_d_o2(x1))*c;
    q_o1_d_o3(k+1) = q_o1_d_o3(k)+(beta_o1_d_o3*alpha_arrive_o1_d(k)-alpha_leave_o1_d_o3(x1))*c;   

    value(k+1) = fval1;
    green_light(k+1) = x1;
end

f = sum(value)

figure()
hold on
grid on
plot(green_light(1:61)*3600);
plot(green_light(1:61)*3600);
plot(ones(1,61)*c*3600);
plot(ones(1,61)*c*3600);
plot(c*3600-green_light(1:61)*3600);
plot(c*3600-green_light(1:61)*3600);
xlim([1 61]);
ylim([0 65]);
legend('$g_{u,d,o1}$','$g_{u,d,o2}$', '$g_{u,d,o3}$', '$g_{o1,d,u}$','$g_{o1,d,o2}$','$g_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Length of the green time (sec)');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=15s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;

figure()
plot(iteration(1:61),"Color",'m')
xlim([1 61]);
ylabel('$y(k)$','interpreter','latex')
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=15s \forall k$, for $i\in O_{u,d}$','interpreter','latex')

figure();
hold on
grid on
plot(n_u_d_(1:61)*c);
plot(n_o1_d_(1:61)*c);
xlim([1 61]);
legend('$n_{u,d}$', '$n_{o1,d}$','interpreter','latex','location',"best")
ylabel('Number of vehicles on the links');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=15s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;

figure();
hold on
grid on
plot(q_u_d_o1(1:61));
plot(q_u_d_o2(1:61));
plot(q_u_d_o3(1:61));
plot(q_o1_d_u(1:61));
plot(q_o1_d_o2(1:61));
plot(q_o1_d_o3(1:61));
xlim([1 61]);
legend('$q_{u,d,o1}$', '$q_{u,d,o2}$', '$q_{u,d,o3}$', '$q_{o1,d,u}$', '$q_{o1,d,o2}$', '$q_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Queue length of each lane (veh)');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=15s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;
