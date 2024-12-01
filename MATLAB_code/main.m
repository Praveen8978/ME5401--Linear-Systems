clear all;
clc;

%% Given Information
% Given matrices
A = [-1.7, -0.25, 0;
    23, -30, 20;
    0, -660, -860];
B = [7, 0;
    -118, 0;
    0, -1300];
C = [0, 1, 0;
    0, 0, 1];
D = [0 0;
    0 0];
% initial conditions
x0 = [1 100 200]';

% calculated zeta and omega based on time domain specifications
zeta = 0.8;
omega = 0.2;

%% pole placement 
[p1, p2] = pole_estimator(zeta, omega); % 3rd pole is assumed to be -1

P = [p1 p2 -1];
% K = pole_placement(A,B,P);

% K is used to simulate the system using zero_input_pole_placement.slx

%% LQR
% define Q and R matrices

Q = [15 0 0;
    0 900 0;
    0 0 150];
R = [200 0;
    0 160];

[K, V, U, P, M] = LQR_own(A,B,Q,R);

% K is used to simulate the system using step_input_lqr.slx and
% zero_input_lqr.slx

%% LQR controller with observer

% check for observability of pair(A,C)
obsr = [C; C*A; C*A*A];
if rank(obsr) < 3
        error('it is unobservable!');
end

% We use pole placement to find observer L
% observer poles:
P = [-0.5 -0.5 -1];
% L = pole_placement(A,B,P);

% L and K from LQR are used to simulate using observer_controller.slx

%% Decoupling control

[F,K] = decoupling(A,B,C);

% these F,K values combined with K from LQR are used to simulate using
% decoupler_ss.slx

%% Intergal control + state observer

y_sp=[100;150]; % settling point
disturbance=[-2;5]; % disturbance load applied at t=10s

Contr = rank([A B;C zeros(2,2)]);
if Contr ~= size(A,1) + size(C,1)
    error("The matrix is un-controllable");
end

% defining augmented matrices
A_bar=[A zeros(3,2);-C zeros(2,2)];
B_bar=[B;zeros(2,2)];
B_w_bar=[B;zeros(2,2)];
B_r_bar=[zeros(3,2);eye(2)];
C_bar=[C zeros(2,2)];

% defining Q and R to solve the augmented system
Q=[1 0 0 0 0;
   0 1 0 0 0;
   0 0 1 0 0; 
   0 0 0 1 0; 
   0 0 0 0 1]*10;

R=[1 0
   0 1]*0.5;

% solving using LQR
gamma=[A_bar -B_bar/R*B_bar'
    -Q -A_bar'];
[vector,value]=eig(gamma);

value=sum(value);
v=vector(:,find(real(value)<0));

P=v(6:10,:)/v(1:5,:);
K=real(inv(R)*B_bar'*P);

K1 = K(:,1:3);
K2 = K(:,4:5);

% to design the observer
% controller poles
desired_poles = [-20000 -5 -2000];
L = place(A',C',desired_poles)';

% using K1, K2, L the system is simulated using servos.slx

%% Question 6

% calculating given weight matrix
a=4;
b=6;
c=8;
d=8;
w=[a+b+1 0 0
    0 c+4 0
    0 0 d+5];

syms s u1 u2
U=[u1;u2];
xsp=[5 250 300]'; % settling point
xs=inv(s*eye(3)-A)*B*U;
xs=subs(xs,s,0); % steady state value of the system

% solving the jacobian to find the u1,u2 that minimize the loss function.
J=1/2*(xs-xsp)'*w*(xs-xsp);
jacob = jacobian(J, [u1 u2]);
ans = vpasolve(jacob==0);
u1=double(ans.u1);
u2=double(ans.u2);
u=[u1;u2];

% This u is used to simulate the system using state_control_last_model.slx










