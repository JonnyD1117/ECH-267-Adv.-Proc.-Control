%% 3DOF RRR Robot Dynamics Derivation 
clear all 
clc 

syms t q1(t) q2(t) q3(t) L1 L2 L3 g m1 m2 m3 I1 I2 I3 


q1_dot = diff(q1,t) ; 
q2_dot = diff(q2,t); 
q3_dot = diff(q3,t); 



% Lagrange Equation for 

% Link Mass Carteasian Coordinates
y1 = L1; 
y2 = y1 + (.5*L2)*sin(q2); 
y3 = y2 + (.5*L3)*sin(q2 + q3); 

x1 = 0;
x2 = (.5*L2)*cos(q2);
x3 = L2*cos(q2) + (.5*L3)*cos(q2 + q3);

% Link Mass Carteasian Velocities
y1_dot = 0; 
y2_dot = .5*L2*cos(q2)*q2_dot ; 
y3_dot = L2*cos(q2)*q2_dot + (.5*L3)*cos(q2 + q3)*(q2_dot + q3_dot); 

x1_dot = 0; 
x2_dot = -0.5*L2*sin(q2)*q2_dot; 
x3_dot = -L2*sin(q2)*q2_dot -(.5*L3)*sin(q2+q3)*(q2_dot + q3_dot); 

v1_sqr = x1_dot^2 + y1_dot^2; 
v2_sqr = x2_dot^2 + y2_dot^2; 
v3_sqr = x3_dot^2 + y3_dot^2; 

omega_1 = q1_dot; 
omega_2 = q2_dot;
omega_3 = q2_dot + q3_dot;

% Define Potential Energy for Each Link 
U1 = m1*g*y1; 
U2 = m2*g*y2; 
U3 = m3*g*y3; 
U_tot = U1 + U2 + U3;


% Define Kinetic Energy for Each Link
k1 = .5*m1*(v1_sqr) + .5*I1*(omega_1)^2;
k2 = .5*m1*(v1_sqr) + .5*I2*(omega_2)^2;
k3 = .5*m1*(v1_sqr) + .5*I3*(omega_3)^2;

K_tot = k1 + k2 + k3; 


L = K_tot - U_tot;

latex(L)





