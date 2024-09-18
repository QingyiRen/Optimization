clear all 
close all

%% parameters
global func;
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
lb = ones(1,61)*15/3600;
ub = ones(1,61)*45/3600;

rng default % For reproducibility
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',8e4,'TolFun',5e-2);
[x1,fval1,exitflag1,output1] = fmincon(@(x1)myfun(x1,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d),lb,[],[],[],[],lb,ub,[],options)

figure()
hold on
grid on
plot(x1(1:61)*3600);
plot(x1(1:61)*3600);
plot(ones(1,61)*c*3600);
plot(ones(1,61)*c*3600);
plot(c*3600-x1(1:61)*3600);
plot(c*3600-x1(1:61)*3600);
xlim([0 61]);
ylim([0 65]);
legend('$g_{u,d,o1}$','$g_{u,d,o2}$', '$g_{u,d,o3}$', '$g_{o1,d,u}$','$g_{o1,d,o2}$','$g_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Length of the green time (sec)');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=15s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;

figure()
plot(func,"Color",'m')
xlim([0 61]);
ylabel('$y(k)$','interpreter','latex')
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=15s \forall k$, for $i\in O_{u,d}$','interpreter','latex')

global n_u_d n_o1_d q_u_d_o1 q_u_d_o2 q_u_d_o3 q_o1_d_u q_o1_d_o2 q_o1_d_o3;
figure();
hold on
grid on
plot(n_u_d(1:end)*c);
plot(n_o1_d(1:end)*c);
xlim([0 61]);
legend('$n_{u,d}$', '$n_{o1,d}$','interpreter','latex','location',"best")
ylabel('Number of vehicles on the links');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=15s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;

figure();
hold on
grid on
plot(q_u_d_o1(1:end));
plot(q_u_d_o2(1:end));
plot(q_u_d_o3(1:end));
plot(q_o1_d_u(1:end));
plot(q_o1_d_o2(1:end));
plot(q_o1_d_o3(1:end));
xlim([0 61]);
legend('$q_{u,d,o1}$', '$q_{u,d,o2}$', '$q_{u,d,o3}$', '$q_{o1,d,u}$', '$q_{o1,d,o2}$', '$q_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Queue length of each lane (veh)');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=15s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;

rng default % For reproducibility
[x2,fval2,exitflag2,output2] = fmincon(@(x2)myfun(x2,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d),ub,[],[],[],[],lb,ub,[],options)

figure()
hold on
grid on
plot(x2(1:61)*3600);
plot(x2(1:61)*3600);
plot(ones(1,61)*c*3600);
plot(ones(1,61)*c*3600);
plot(c*3600-x2(1:61)*3600);
plot(c*3600-x2(1:61)*3600);
xlim([0 61]);
ylim([0 65]);
legend('$g_{u,d,o1}$','$g_{u,d,o2}$', '$g_{u,d,o3}$', '$g_{o1,d,u}$','$g_{o1,d,o2}$','$g_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Length of the green time (sec)');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=45s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;

figure()
plot(func,"Color",'m')
xlim([0 61]);
ylabel('$y(k)$','interpreter','latex')
xlabel('simulation cycle')
title('Starting point with $g_{u,d,i(k)}=45s \forall k$, for $i\in O_{u,d}$','interpreter','latex')

global n_u_d n_o1_d q_u_d_o1 q_u_d_o2 q_u_d_o3 q_o1_d_u q_o1_d_o2 q_o1_d_o3;
figure();
hold on
grid on
plot(n_u_d(1:end)*c);
plot(n_o1_d(1:end)*c);
xlim([0 61]);
legend('$n_{u,d}$', '$n_{o1,d}$','interpreter','latex','location',"best")
ylabel('Number of vehicles on the links');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=45s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;

figure();
hold on
grid on
plot(q_u_d_o1(1:end));
plot(q_u_d_o2(1:end));
plot(q_u_d_o3(1:end));
plot(q_o1_d_u(1:end));
plot(q_o1_d_o2(1:end));
plot(q_o1_d_o3(1:end));
xlim([0 61]);
legend('$q_{u,d,o1}$', '$q_{u,d,o2}$', '$q_{u,d,o3}$', '$q_{o1,d,u}$', '$q_{o1,d,o2}$', '$q_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Queue length of each lane (veh)');
xlabel('simulation cycle');
title('Starting point with $g_{u,d,i(k)}=45s \forall k$, for $i\in O_{u,d}$','interpreter','latex')
hold off;

%% question 4
global nocont;
rng default % For reproducibility
[x3,fval3,exitflag3,output3] = fmincon(@(x3)nocontrol(x3,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d),ones(1,61)*30/3600,[],[],[],[],lb,ub,[],options)

figure()
hold on
grid on
plot(x3(1:61)*3600);
plot(x3(1:61)*3600);
plot(ones(1,61)*c*3600);
plot(ones(1,61)*c*3600);
plot(c*3600-x3(1:61)*3600);
plot(c*3600-x3(1:61)*3600);
xlim([0 61]);
ylim([0 65]);
legend('$g_{u,d,o1}$','$g_{u,d,o2}$', '$g_{u,d,o3}$', '$g_{o1,d,u}$','$g_{o1,d,o2}$','$g_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Length of the green time (sec)');
xlabel('simulation cycle');
title('No-Control case ($g_{u,d,i(k)}=g_{o_1,d,j(k)}=30s \forall k$, for $i\in O_{u,d},j \in O_{o_1,d}$)','interpreter','latex')
hold off;

figure()
plot(nocont,"Color",'m')
xlim([0 61]);
title('No-Control case ($g_{u,d,i(k)}=g_{o_1,d,j(k)}=30s \forall k$, for $i\in O_{u,d},j \in O_{o_1,d}$)','interpreter','latex')
ylabel('$y(k)$','interpreter','latex')
xlabel('simulation cycle');

global n_u_d n_o1_d q_u_d_o1 q_u_d_o2 q_u_d_o3 q_o1_d_u q_o1_d_o2 q_o1_d_o3;
figure();
hold on
grid on
plot(n_u_d(1:end)*c);
plot(n_o1_d(1:end)*c);
xlim([0 61]);
legend('$n_{u,d}$', '$n_{o1,d}$','interpreter','latex','location',"best")
ylabel('Number of vehicles on the links');
xlabel('simulation cycle');
title('No-Control case ($g_{u,d,i(k)}=g_{o_1,d,j(k)}=30s \forall k$, for $i\in O_{u,d},j \in O_{o_1,d}$)','interpreter','latex')
hold off;

figure();
hold on
grid on
plot(q_u_d_o1(1:end));
plot(q_u_d_o2(1:end));
plot(q_u_d_o3(1:end));
plot(q_o1_d_u(1:end));
plot(q_o1_d_o2(1:end));
plot(q_o1_d_o3(1:end));
xlim([0 61]);
legend('$q_{u,d,o1}$', '$q_{u,d,o2}$', '$q_{u,d,o3}$', '$q_{o1,d,u}$', '$q_{o1,d,o2}$', '$q_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Queue length of each lane (veh)');
xlabel('simulation cycle');
title('No-Control case ($g_{u,d,i(k)}=g_{o_1,d,j(k)}=30s \forall k$, for $i\in O_{u,d},j \in O_{o_1,d}$)','interpreter','latex')
hold off;

%% question 5
lb=3*ones(1,61);
ub=9*ones(1,61);
opts = optimoptions('ga','ConstraintTolerance',1e-6, 'MaxStallGenerations',1000,'MaxGenerations',1000,'InitialPopulationMatrix', lb, 'PlotFcn', @gaplotbestf);
rng default % For reproducibility
[x4,fval4,eflag4,outpt4] = ga(@(x4)integer_sol(x4,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d),61,[],[],[],[],lb,ub,[],1:61,opts)

global integ n_u_d n_o1_d q_u_d_o1 q_u_d_o2 q_u_d_o3 q_o1_d_u q_o1_d_o2 q_o1_d_o3;

figure()
hold on
grid on
plot(x4(1:61)*5);
plot(x4(1:61)*5);
plot(ones(1,61)*c*3600);
plot(ones(1,61)*c*3600);
plot(c*3600-x4(1:61)*5);
plot(c*3600-x4(1:61)*5);
xlim([0 61]);
ylim([0 65]);
legend('$g_{u,d,o1}$','$g_{u,d,o2}$', '$g_{u,d,o3}$', '$g_{o1,d,u}$','$g_{o1,d,o2}$','$g_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Length of the green time (sec)');
xlabel('simulation cycle');
title('Genetic algorithm')
hold off;

figure()
plot(integ,"Color",'m')
xlim([0 61]);
title('Genetic algorithm')
ylabel('$y(k)$','interpreter','latex')
xlabel('simulation cycle');

figure();
hold on
grid on
plot(n_u_d(1:end)*c);
plot(n_o1_d(1:end)*c);
xlim([0 61]);
legend('$n_{u,d}$', '$n_{o1,d}$','interpreter','latex','location',"best")
ylabel('Number of vehicles on the links');
xlabel('simulation cycle');
title('Genetic algorithm')
hold off;

figure();
hold on
grid on
plot(q_u_d_o1(1:end));
plot(q_u_d_o2(1:end));
plot(q_u_d_o3(1:end));
plot(q_o1_d_u(1:end));
plot(q_o1_d_o2(1:end));
plot(q_o1_d_o3(1:end));
xlim([0 61]);
legend('$q_{u,d,o1}$', '$q_{u,d,o2}$', '$q_{u,d,o3}$', '$q_{o1,d,u}$', '$q_{o1,d,o2}$', '$q_{o1,d,o3}$','interpreter','latex','location',"best")
ylabel('Queue length of each lane (veh)');
xlabel('simulation cycle');
title('Genetic algorithm')
hold off;
%% functions used

function i = integer_sol(x,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d)
    
    global integ n_u_d n_o1_d q_u_d_o1 q_u_d_o2 q_u_d_o3 q_o1_d_u q_o1_d_o2 q_o1_d_o3;

    % parameters for link (u,d)
    c = 60/3600;
    N_u_d = 3;
    v_u_d = 60; 
    l_veh = 7/1000;
    beta_u_d_o1 = 0.4;
    beta_u_d_o2 = 0.3;
    beta_u_d_o3 = 0.3;
    mu_u_d = 1800;
    
    % parameters for link (o1,d)
    N_o1_d = 3;
    v_o1_d = 50; 
    beta_o1_d_u = 0.3;
    beta_o1_d_o2 = 0.4;
    beta_o1_d_o3 = 0.3;
    mu_o1_d = 1700;
    
    x = x*5/3600; % for values between 15s and 45s

    g_1 = x(1:61);
    g_2 = x(1:61);
    g_3 = ones(1,61)*c;
    g_4 = ones(1,61)*c;
    g_5 = c-x(1:61);
    g_6 = c-x(1:61);

    for k = 1:61

        if k==1
            n_u_d(k) = 0;
            n_o1_d(k) = 0;

            q_u_d_o1(k) = 0;
            q_u_d_o2(k) = 0;
            q_u_d_o3(k) = 0;
    
            q_o1_d_u(k) = 0;
            q_o1_d_o2(k) = 0;
            q_o1_d_o3(k) = 0;

        else
            n_u_d(k)=n_u_d(k-1)+(alpha_enter_u_d(k-1)-(alpha_leave_u_d_o1(k-1)+alpha_leave_u_d_o2(k-1)+alpha_leave_u_d_o3(k-1)))*c;
            n_o1_d(k)=n_o1_d(k-1)+(alpha_enter_o1_d(k-1)-(alpha_leave_o1_d_u(k-1)+alpha_leave_o1_d_o2(k-1)+alpha_leave_o1_d_o3(k-1)))*c;

            q_u_d_o1(k) = q_u_d_o1(k-1)+(beta_u_d_o1*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o1(k-1))*c;
            q_u_d_o2(k) = q_u_d_o2(k-1)+(beta_u_d_o2*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o2(k-1))*c;
            q_u_d_o3(k) = q_u_d_o3(k-1)+(beta_u_d_o3*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o3(k-1))*c;
    
            q_o1_d_u(k) = q_o1_d_u(k-1)+(beta_o1_d_u*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_u(k-1))*c;
            q_o1_d_o2(k) = q_o1_d_o2(k-1)+(beta_o1_d_o2*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_o2(k-1))*c;
            q_o1_d_o3(k) = q_o1_d_o3(k-1)+(beta_o1_d_o3*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_o3(k-1))*c;      

        end

        tau_u_d(k) = floor((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh/(N_u_d*v_u_d*c));
        tau_o1_d(k) = floor((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh/(N_o1_d*v_o1_d*c));     

        gamma_u_d(k) = rem((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh, N_u_d*v_u_d*c)/(N_u_d*v_u_d*c);
        gamma_o1_d(k) = rem((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh, N_o1_d*v_o1_d*c)/(N_o1_d*v_o1_d*c);

        if k-tau_u_d(k)<1
            alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(1)+gamma_u_d(k)*alpha_enter_u_d(1);
        else
            if k-tau_u_d(k)-1<1
               alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(1);
            else
               alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(k-tau_u_d(k)-1);
            end
        end

        if k-tau_o1_d(k)<1
            alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(1)+gamma_o1_d(k)*alpha_enter_o1_d(1);
        else
            if k-tau_o1_d(k)-1<1
               alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(1);
            else
               alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(k-tau_o1_d(k)-1);
            end
        end

        alpha_leave_u_d_o1(k) = min([beta_u_d_o1*mu_u_d*g_1(k)/c, q_u_d_o1(k)/c + beta_u_d_o1*alpha_arrive_u_d(k), C_d_o1(k)/c]);
        alpha_leave_u_d_o2(k) = min([beta_u_d_o2*mu_u_d*g_2(k)/c, q_u_d_o2(k)/c + beta_u_d_o2*alpha_arrive_u_d(k), C_d_o2(k)/c]);
        alpha_leave_u_d_o3(k) = min([beta_u_d_o3*mu_u_d*g_3(k)/c, q_u_d_o3(k)/c + beta_u_d_o3*alpha_arrive_u_d(k), C_d_o3(k)/c]);
    
        alpha_leave_o1_d_u(k) = min([beta_o1_d_u*mu_o1_d*g_4(k)/c, q_o1_d_u(k)/c + beta_o1_d_u*alpha_arrive_o1_d(k), C_d_u(k)/c]);
        alpha_leave_o1_d_o2(k) = min([beta_o1_d_o2*mu_o1_d*g_5(k)/c, q_o1_d_o2(k)/c + beta_o1_d_o2*alpha_arrive_o1_d(k), C_d_o2(k)/c]);
        alpha_leave_o1_d_o3(k) = min([beta_o1_d_o3*mu_o1_d*g_6(k)/c, q_o1_d_o3(k)/c + beta_o1_d_o3*alpha_arrive_o1_d(k), C_d_o3(k)/c]);

        integ(k) = (n_u_d(k)+n_o1_d(k))*c;
    end
    i = sum(integ);
end

function m = nocontrol(x,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d)
    
    global nocont n_u_d n_o1_d q_u_d_o1 q_u_d_o2 q_u_d_o3 q_o1_d_u q_o1_d_o2 q_o1_d_o3;

    % parameters for link (u,d)
    c = 60/3600;
    N_u_d = 3;
    v_u_d = 60; 
    l_veh = 7/1000;
    beta_u_d_o1 = 0.4;
    beta_u_d_o2 = 0.3;
    beta_u_d_o3 = 0.3;
    mu_u_d = 1800;
    
    % parameters for link (o1,d)
    N_o1_d = 3;
    v_o1_d = 50; 
    beta_o1_d_u = 0.3;
    beta_o1_d_o2 = 0.4;
    beta_o1_d_o3 = 0.3;
    mu_o1_d = 1700;
    
    g_1 = ones(1,61)*30/3600;
    g_2 = ones(1,61)*30/3600;
    g_3 = ones(1,61)*30/3600;
    g_4 = ones(1,61)*30/3600;
    g_5 = ones(1,61)*30/3600;
    g_6 = ones(1,61)*30/3600;
    
    for k = 1:61

        if k==1
            n_u_d(k) = 0;
            n_o1_d(k) = 0;

            q_u_d_o1(k) = 0;
            q_u_d_o2(k) = 0;
            q_u_d_o3(k) = 0;
    
            q_o1_d_u(k) = 0;
            q_o1_d_o2(k) = 0;
            q_o1_d_o3(k) = 0;

        else
            n_u_d(k)=n_u_d(k-1)+(alpha_enter_u_d(k-1)-(alpha_leave_u_d_o1(k-1)+alpha_leave_u_d_o2(k-1)+alpha_leave_u_d_o3(k-1)))*c;
            n_o1_d(k)=n_o1_d(k-1)+(alpha_enter_o1_d(k-1)-(alpha_leave_o1_d_u(k-1)+alpha_leave_o1_d_o2(k-1)+alpha_leave_o1_d_o3(k-1)))*c;

            q_u_d_o1(k) = q_u_d_o1(k-1)+(beta_u_d_o1*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o1(k-1))*c;
            q_u_d_o2(k) = q_u_d_o2(k-1)+(beta_u_d_o2*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o2(k-1))*c;
            q_u_d_o3(k) = q_u_d_o3(k-1)+(beta_u_d_o3*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o3(k-1))*c;
    
            q_o1_d_u(k) = q_o1_d_u(k-1)+(beta_o1_d_u*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_u(k-1))*c;
            q_o1_d_o2(k) = q_o1_d_o2(k-1)+(beta_o1_d_o2*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_o2(k-1))*c;
            q_o1_d_o3(k) = q_o1_d_o3(k-1)+(beta_o1_d_o3*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_o3(k-1))*c;      

        end

        tau_u_d(k) = floor((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh/(N_u_d*v_u_d*c));
        tau_o1_d(k) = floor((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh/(N_o1_d*v_o1_d*c));     

        gamma_u_d(k) = rem((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh, N_u_d*v_u_d*c)/(N_u_d*v_u_d*c);
        gamma_o1_d(k) = rem((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh, N_o1_d*v_o1_d*c)/(N_o1_d*v_o1_d*c);

        if k-tau_u_d(k)<1
            alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(1)+gamma_u_d(k)*alpha_enter_u_d(1);
        else
            if k-tau_u_d(k)-1<1
               alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(1);
            else
               alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(k-tau_u_d(k)-1);
            end
        end

        if k-tau_o1_d(k)<1
            alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(1)+gamma_o1_d(k)*alpha_enter_o1_d(1);
        else
            if k-tau_o1_d(k)-1<1
               alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(1);
            else
               alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(k-tau_o1_d(k)-1);
            end
        end

        alpha_leave_u_d_o1(k) = min([beta_u_d_o1*mu_u_d*g_1(k)/c, q_u_d_o1(k)/c + beta_u_d_o1*alpha_arrive_u_d(k), C_d_o1(k)/c]);
        alpha_leave_u_d_o2(k) = min([beta_u_d_o2*mu_u_d*g_2(k)/c, q_u_d_o2(k)/c + beta_u_d_o2*alpha_arrive_u_d(k), C_d_o2(k)/c]);
        alpha_leave_u_d_o3(k) = min([beta_u_d_o3*mu_u_d*g_3(k)/c, q_u_d_o3(k)/c + beta_u_d_o3*alpha_arrive_u_d(k), C_d_o3(k)/c]);
    
        alpha_leave_o1_d_u(k) = min([beta_o1_d_u*mu_o1_d*g_4(k)/c, q_o1_d_u(k)/c + beta_o1_d_u*alpha_arrive_o1_d(k), C_d_u(k)/c]);
        alpha_leave_o1_d_o2(k) = min([beta_o1_d_o2*mu_o1_d*g_5(k)/c, q_o1_d_o2(k)/c + beta_o1_d_o2*alpha_arrive_o1_d(k), C_d_o2(k)/c]);
        alpha_leave_o1_d_o3(k) = min([beta_o1_d_o3*mu_o1_d*g_6(k)/c, q_o1_d_o3(k)/c + beta_o1_d_o3*alpha_arrive_o1_d(k), C_d_o3(k)/c]);

        nocont(k) = (n_u_d(k)+n_o1_d(k))*c;
    end
    m = sum(nocont);
end

function f = myfun(x,C_u_d,C_o1_d,C_d_u,C_d_o3,C_d_o2,C_d_o1,alpha_enter_o1_d,alpha_enter_u_d)

    global func n_u_d n_o1_d q_u_d_o1 q_u_d_o2 q_u_d_o3 q_o1_d_u q_o1_d_o2 q_o1_d_o3;

    % parameters for link (u,d)
    c = 60/3600;
    N_u_d = 3;
    v_u_d = 60; 
    l_veh = 7/1000;
    beta_u_d_o1 = 0.4;
    beta_u_d_o2 = 0.3;
    beta_u_d_o3 = 0.3;
    mu_u_d = 1800;
    
    % parameters for link (o1,d)
    N_o1_d = 3;
    v_o1_d = 50; 
    beta_o1_d_u = 0.3;
    beta_o1_d_o2 = 0.4;
    beta_o1_d_o3 = 0.3;
    mu_o1_d = 1700;
    
    g_1 = x(1:61);
    g_2 = x(1:61);
    g_3 = ones(1,61)*c;
    g_4 = ones(1,61)*c;
    g_5 = c-x(1:61);
    g_6 = c-x(1:61);
    
    for k = 1:61

        if k==1
            n_u_d(k) = 0;
            n_o1_d(k) = 0;

            q_u_d_o1(k) = 0;
            q_u_d_o2(k) = 0;
            q_u_d_o3(k) = 0;
    
            q_o1_d_u(k) = 0;
            q_o1_d_o2(k) = 0;
            q_o1_d_o3(k) = 0;

        else
            n_u_d(k)=n_u_d(k-1)+(alpha_enter_u_d(k-1)-(alpha_leave_u_d_o1(k-1)+alpha_leave_u_d_o2(k-1)+alpha_leave_u_d_o3(k-1)))*c;
            n_o1_d(k)=n_o1_d(k-1)+(alpha_enter_o1_d(k-1)-(alpha_leave_o1_d_u(k-1)+alpha_leave_o1_d_o2(k-1)+alpha_leave_o1_d_o3(k-1)))*c;

            q_u_d_o1(k) = q_u_d_o1(k-1)+(beta_u_d_o1*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o1(k-1))*c;
            q_u_d_o2(k) = q_u_d_o2(k-1)+(beta_u_d_o2*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o2(k-1))*c;
            q_u_d_o3(k) = q_u_d_o3(k-1)+(beta_u_d_o3*alpha_arrive_u_d(k-1)-alpha_leave_u_d_o3(k-1))*c;
    
            q_o1_d_u(k) = q_o1_d_u(k-1)+(beta_o1_d_u*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_u(k-1))*c;
            q_o1_d_o2(k) = q_o1_d_o2(k-1)+(beta_o1_d_o2*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_o2(k-1))*c;
            q_o1_d_o3(k) = q_o1_d_o3(k-1)+(beta_o1_d_o3*alpha_arrive_o1_d(k-1)-alpha_leave_o1_d_o3(k-1))*c;      

        end

        tau_u_d(k) = floor((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh/(N_u_d*v_u_d*c));
        tau_o1_d(k) = floor((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh/(N_o1_d*v_o1_d*c));     

        gamma_u_d(k) = rem((C_u_d-q_u_d_o1(k)-q_u_d_o2(k)-q_u_d_o3(k))*l_veh, N_u_d*v_u_d*c)/(N_u_d*v_u_d*c);
        gamma_o1_d(k) = rem((C_o1_d-q_o1_d_u(k)-q_o1_d_o2(k)-q_o1_d_o3(k))*l_veh, N_o1_d*v_o1_d*c)/(N_o1_d*v_o1_d*c);

        if k-tau_u_d(k)<1
            alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(1)+gamma_u_d(k)*alpha_enter_u_d(1);
        else
            if k-tau_u_d(k)-1<1
               alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(1);
            else
               alpha_arrive_u_d(k) = (1-gamma_u_d(k))*alpha_enter_u_d(k-tau_u_d(k))+gamma_u_d(k)*alpha_enter_u_d(k-tau_u_d(k)-1);
            end
        end

        if k-tau_o1_d(k)<1
            alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(1)+gamma_o1_d(k)*alpha_enter_o1_d(1);
        else
            if k-tau_o1_d(k)-1<1
               alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(1);
            else
               alpha_arrive_o1_d(k) = (1-gamma_o1_d(k))*alpha_enter_o1_d(k-tau_o1_d(k))+gamma_o1_d(k)*alpha_enter_o1_d(k-tau_o1_d(k)-1);
            end
        end

        alpha_leave_u_d_o1(k) = min([beta_u_d_o1*mu_u_d*g_1(k)/c, q_u_d_o1(k)/c + beta_u_d_o1*alpha_arrive_u_d(k), C_d_o1(k)/c]);
        alpha_leave_u_d_o2(k) = min([beta_u_d_o2*mu_u_d*g_2(k)/c, q_u_d_o2(k)/c + beta_u_d_o2*alpha_arrive_u_d(k), C_d_o2(k)/c]);
        alpha_leave_u_d_o3(k) = min([beta_u_d_o3*mu_u_d*g_3(k)/c, q_u_d_o3(k)/c + beta_u_d_o3*alpha_arrive_u_d(k), C_d_o3(k)/c]);
    
        alpha_leave_o1_d_u(k) = min([beta_o1_d_u*mu_o1_d*g_4(k)/c, q_o1_d_u(k)/c + beta_o1_d_u*alpha_arrive_o1_d(k), C_d_u(k)/c]);
        alpha_leave_o1_d_o2(k) = min([beta_o1_d_o2*mu_o1_d*g_5(k)/c, q_o1_d_o2(k)/c + beta_o1_d_o2*alpha_arrive_o1_d(k), C_d_o2(k)/c]);
        alpha_leave_o1_d_o3(k) = min([beta_o1_d_o3*mu_o1_d*g_6(k)/c, q_o1_d_o3(k)/c + beta_o1_d_o3*alpha_arrive_o1_d(k), C_d_o3(k)/c]);

        func(k) = (n_u_d(k)+n_o1_d(k))*c;
    end

    f = sum(func);
end 