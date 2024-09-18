%--------- LP and QP programming assignment--------%
% Maria de Neves de Fonseca, 5611679
% Qingyi Ren, 5684803

%Note: to not have any kind of interferance in the results from one
%question to the other, please run one section at the time, even though the code runs either way. In every
%question all variables are cleared before starting the code.

%% Question 1.1

clear all
close all

c=[-3.6 -2.2];
A=[250 120 
   1 1];
b=[2000+20*14 10]';

lb=zeros(2,1);
ub=[];             % No upper bound
[x,fval,exitflag]=linprog(c,A,b,[],[],lb,ub)
rate=-c*x 

%There might be four pairs of integer values of (G,B):(8,1),(8,2),(9,1),(9,2)
%First we have to check whether or not the pair of integer values of (G,B) satisfies the contraints
p=size(A);
constr_num=p(1);%constr_num means how many constraints inequality we have
pair=[8 8 9 9;
      1 2 1 2];
for j=1:4 % we have to check four pairs of integer values
  for i=1:constr_num
    if A(i,:)*pair(:,j)>b(i)% the jth pair doesn't satisfy the ith constraint inequality
      disp('pair');
      disp(j);
      disp('does not satisfy constraint');
      disp(i);
    end
  end
end
% After the constraint check, only the first pair and the second pair
% satisfy the constraints inequality

%First pair (8,1)
G=8;
B=1;
rate1=-c*[G;B];

%Second pair (8,2)
G=8;
B=2;
rate2=-c*[G;B];
poss_rate=[rate1 rate2];
[max_rate,max_pair_number]=max(poss_rate);

%max_pair_number=2, we should choose (G,B)=(8,2)
%The fastest mining rate of DelftCoin is max_rate=33.2000

%% Question 1.3

clear all
close all

c=[250 120];
A=[250 120 
   1 1];
b=[2000+20*14 10]';

Aeq=[3.6*365 2.2*365];
beq=8000+20*7;

lb=zeros(2,1);
ub=[];             % No upper bound
[x,fval,exitflag]=linprog(c,A,b,Aeq,beq,lb,ub)
spent=c*x

%% Question 3

clear all
close all

T = readtable("measurements22.csv");

delta_t = 60;
T_k = T.Var1;
q_in = T.Var3;
q_out = T.Var4;
T_amb = T.Var2;
T_k_1 = T_k(2:end);

Y = -T_k_1+T_k(1:end-1);
phi = [(T_k(1:end-1)-T_amb(1:end-1))*delta_t (-q_in(1:end-1)+q_out(1:end-1))*delta_t];
H = 2*phi'*phi;
c = -2*phi'*Y;

% -------check if H is symmetric positive definite------%
try chol(H)
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end

%since there are no constraints we can obtain x as follows
x = -inv(H)*c

%% Question 4

clear all
close all

T = readtable("measurements22.csv");

%variables
delta_t = 60;
T_k = T.Var1;
q_in = T.Var3;
q_out = T.Var4;
T_amb = T.Var2;
T_k_1 = T_k(2:end);
phi = T.Var5/3600000;
E3 = 12;
psi = (0.38+0.01*E3)/3600000;
N = 1440;
T_max = 95;
q_in_max = 125;
T_1 = 23.3;
a_1 = 0.0005;
a_2 =  0.0020;

%optimization problem
c=[zeros(1,N), delta_t*(phi(1:N)'-psi)];

D=zeros(1,N-1);
for i=1:numel(D)
    D(i)=-(1-a_1*delta_t);
end

J=zeros(1,N);
for i=1:numel(J)
    J(i)=-a_2*delta_t;
end

A_eq_1 = diag(ones(1,N))+diag(D,-1);
A_eq_2 = diag(J);

A_eq = [A_eq_1 A_eq_2];

b_eq_1=zeros(1,N);
for i=1:numel(b_eq_1)
    b_eq_1(i)=-a_2*delta_t*q_out(i)+a_1*delta_t*T_amb(i);
end

b_eq_2 = [(1-a_1*delta_t)*T_1, zeros(1,N-1)]; 

b_eq = (b_eq_1+b_eq_2)';

A = [eye(N) zeros(N,N);
    zeros(N,N) eye(N);
    zeros(N,N) -eye(N)];

b = [T_max*ones(N,1);q_in_max*ones(N,1); zeros(N,1)];        

[x,fval,exitflag]=linprog(c,A,b,A_eq,b_eq,[],[]);
total_cost = c*x + sum(phi(1:N).*q_out(1:N)*delta_t)

% Behaviour of the internal temperature and the power used by the mining device for computation 
% figure(1);
% plot(x(1:N),'b-.')
% hold on
% plot(x(N+1:2*N),'r')
% hold off
% 
% legend('T','$\dot{q}^{in}$','Interpreter','latex','Location','southwest')