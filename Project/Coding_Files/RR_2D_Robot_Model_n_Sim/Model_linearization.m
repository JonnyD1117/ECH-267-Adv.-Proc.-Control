clear all 

% close all 


syms m1 m2 m3 L1 L2 L3 q1 q2 q3 q1_dot q2_dot g_val t1 t2

tau = [t1;t2];

eps = .1 ;

m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(q2) + (L2)^2) + eps, m2*(L1*L2*cos(q2)+ (L2)^2); ...
     m2*(L1*L2*cos(q2) + (L2)^2), m2*(L2)^2+eps];

v = [-m2*L1*L2*sin(q2)*(2*q1_dot*q2_dot + (q2_dot)^2);...
     m2*L1*L2*(q1_dot)^2*sin(q2)] ; 

g = [(m1 + m2)*L1*g_val*cos(q1)+ m2*g_val*L2*cos(q1+ q2); ...
     m2*g_val*L2*cos(q1+q2)];
 

f = [q1_dot*.25;q2_dot*.25];
 
 q_ddot = m\(tau - v -g - f);
 
 
 f1 = q1_dot; 
 f2 = q2_dot;
 f3 = q_ddot(1);
 f4 = q_ddot(2);
 
 d_f1_dx1 = diff(f1,q1);
 d_f1_dx2 = diff(f1,q2);
 d_f1_dx3 = diff(f1,q1_dot);
 d_f1_dx4 = diff(f1,q2_dot);
 

 d_f2_dx1 = diff(f2,q1);
 d_f2_dx2 = diff(f2,q2);
 d_f2_dx3 = diff(f2,q1_dot);
 d_f2_dx4 = diff(f2,q2_dot);
 
 
 d_f3_dx1 = diff(f3,q1);
 d_f3_dx2 = diff(f3,q2);
 d_f3_dx3 = diff(f3,q1_dot);
 d_f3_dx4 = diff(f3,q2_dot);


 d_f4_dx1 = diff(f4,q1);
 d_f4_dx2 = diff(f4,q2);
 d_f4_dx3 = diff(f4,q1_dot);
 d_f4_dx4 = diff(f4,q2_dot);
 
 
 
 A = [d_f1_dx1, d_f1_dx2, d_f1_dx3, d_f1_dx4;...
      d_f2_dx1, d_f2_dx2, d_f2_dx3, d_f2_dx4;...
      d_f3_dx1, d_f3_dx2, d_f3_dx3, d_f3_dx4; ...
      d_f4_dx1, d_f4_dx2, d_f4_dx3, d_f4_dx4];
  
  A = simplify(A); 
  
  
  
 d_f1_du1 = diff(f1,t1);
 d_f1_du2 = diff(f1,t2);

 d_f2_du1 = diff(f2,t1);
 d_f2_du2 = diff(f2,t2);
 
 d_f3_du1 = diff(f3,t1);
 d_f3_du2 = diff(f3,t2);

 d_f4_du1 = diff(f4,t1);
 d_f4_du2 = diff(f4,t2);
 
 B = [d_f1_du1,d_f1_du2;...
      d_f2_du1,d_f2_du2;...
      d_f3_du1,d_f3_du2;...
      d_f4_du1,d_f4_du2];
  
  C = eye(4);


%% Full-State Feedback Regulator

x0 = [0,0,0,0];         % Equilibrium Point
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
states_sym = [q1, q2, q1_dot, q2_dot, 9.82]; % State Symbols

params_val = [.25,.25, .5, .25]; 
% states_val = []


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);


G_ss = ss(A_mat, B_mat, C, zeros(4,2));


K = place(A_mat, B_mat, [-5, -4.5, -1, -1]);

% X_stable = (A_mat-B_mat*K);
% 
% eig(X_stable)

r = [-.5;0];




x_init = [pi;-pi/4; 0; 0] ;
X_stable = ss((A_mat-B_mat*K) + B_mat*r, eye(4,2) , C, zeros(4,2));

figure(5)
initial(X_stable, x_init)
title('Full State Feedback Stabilizing Controller X_0=[pi,-.25*pi ,0,0], X_ref=[0,0,0,0]')




%% Full-State Feedback Steady State Reference 
x0 = [0,0,0,0];         % Equilibrium Point
states_sym = [q1, q2, q1_dot, q2_dot]; % State Symbols
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
params_val = [.25,.25, .5, .25, 9.84];


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);


G_ss = ss(A_mat, B_mat, C, zeros(4,2));


% K = place(A_mat, B_mat, [-5, -4.5, -1, -1]);
K = place(A_mat, B_mat, [-7, -7, -2, -2]);


x_targ = [pi; deg2rad(-45);0;0];

x_init = [0;-pi/2; 0; 0] ;

x_new = [];
u_new = []; 
dt = .1;

for k = 1:1:100
    u = -K*(x_init-x_targ);
    x_dot = A_mat*x_init + B_mat*(u);
    
    x_new(:,k) = x_init + dt.*x_dot;
    u_new(:,k) = u; 
    x_init = x_new(:,k);
end 

figure()
subplot(4,1,1)
plot(x_new(1,:))
ylabel('Q1')
title('Full State Feedback Stabilizing Controller X_0=[0,-.5*pi,0,0], X_{ref}=[pi, -.25*pi,0,0]')

subplot(4,1,2)
plot(x_new(2,:))
ylabel('Q2')

subplot(4,1,3)
plot(x_new(3,:))
ylabel('Q1_{dot}')

subplot(4,1,4)
plot(x_new(4,:))
ylabel('Q2_{dot}')
xlabel('Time')



figure()
subplot(2,1,1)
plot(u_new(1,:))
ylabel('Tau_1')
title('Full State Feedback Stabilizing Controller Torque Inputs')

subplot(2,1,2)
plot(u_new(2,:))
ylabel('Tau_2')
 



%% Linear Quadratic Regulator

x0 = [0,0,0,0];         % Equilibrium Point
states_sym = [q1, q2, q1_dot, q2_dot]; % State Symbols
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
params_val = [.25,.25, .5, .25, 9.84];


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);


G_ss = ss(A_mat, B_mat, C, zeros(4,2));


Q = .5.*eye(4); 
R = 1.*eye(2); 

[K,S,CLP] = lqr(G_ss,Q,R); 

x_targ = [0; 0;0;0];

x_init = [pi;-pi/2; 0; 0] ;

x_new = [];
u_new = []; 
dt = .1;

for k = 1:1:100
    u = -K*(x_init-x_targ);
    x_dot = A_mat*x_init + B_mat*(u);
    
    x_new(:,k) = x_init + dt.*x_dot;
    u_new(:,k) = u; 
    x_init = x_new(:,k);
end 

figure()
subplot(4,1,1)
plot(x_new(1,:))
ylabel('Q1')
title('LQR X_0=[0,-.5*pi,0,0]')

subplot(4,1,2)
plot(x_new(2,:))
ylabel('Q2')

subplot(4,1,3)
plot(x_new(3,:))
ylabel('Q1_{dot}')

subplot(4,1,4)
plot(x_new(4,:))
ylabel('Q2_{dot}')
xlabel('Time')

figure()
subplot(2,1,1)
plot(u_new(1,:))
ylabel('Tau_1')
title('LQR Torque Inputs')

subplot(2,1,2)
plot(u_new(2,:))
ylabel('Tau_2')
 
 


%% LQR State Following 

x0 = [0,0,0,0];         % Equilibrium Point
states_sym = [q1, q2, q1_dot, q2_dot]; % State Symbols
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
params_val = [.25,.25, .5, .25, 9.84];


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);


G_ss = ss(A_mat, B_mat, C, zeros(4,2));


Q = .5.*eye(4); 
R = 1.*eye(2); 

[K,S,CLP] = lqr(G_ss,Q,R); 

x_targ = [pi; deg2rad(-45);0;0];

x_init = [0;-pi/2; 0; 0] ;

x_new = [];
u_new = []; 
dt = .1;

for k = 1:1:100
    u = -K*(x_init-x_targ);
    x_dot = A_mat*x_init + B_mat*(u);
    
    x_new(:,k) = x_init + dt.*x_dot;
    u_new(:,k) = u; 
    x_init = x_new(:,k);
end 

figure()
subplot(4,1,1)
plot(x_new(1,:))
ylabel('Q1')
title('LQR X_0=[0,-.5*pi,0,0], X_{ref}=[pi, -.25*pi,0,0]')

subplot(4,1,2)
plot(x_new(2,:))
ylabel('Q2')

subplot(4,1,3)
plot(x_new(3,:))
ylabel('Q1_{dot}')

subplot(4,1,4)
plot(x_new(4,:))
ylabel('Q2_{dot}')
xlabel('Time')



figure()
subplot(2,1,1)
plot(u_new(1,:))
ylabel('Tau_1')
title('LQR Torque Inputs')

subplot(2,1,2)
plot(u_new(2,:))
ylabel('Tau_2')
 
