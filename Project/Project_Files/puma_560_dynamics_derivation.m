%% PUMA 560 Dynamics Derivation 

Puma_560_parameters


w0 = [0; 0; 0]; 
v0 = [0; 0; 0];

syms q1 q2 q3 q1_dot q2_dot q3_dot g_val

% Gravity Vector 

g_t = [0,0,g_val];


% Linear & Angular Velocities

P_01 = [0;0;0] ;  
v1 = v0 + cross(w0, P_01);
w1 = w0 + q1_dot*R_01(q1)*[0;0;1];

P_12 = [0; .2435; 0];
P_02 = vpa(R_01(q1)*P_12, 4);
v2 = v1 + cross(w1,P_02);
w2 = w1 + q2_dot*R_02(q1, q2)*[0;0;1];

%
P_23 = [.4318; 0 ; -.093];
P_03 = vpa(R_02(q1, q2)*P_23, 4);
v3 = v2 + cross(w2,P_03); 
w3 = simplify(w2 + q3_dot*R_03(q1, q2, q3)*[0;0;1]);

% Linear Velocities about Link CG 
r2 = vpa(R_02(q1, q2)*S2,4);
Vc2 = vpa(simplify(v2 + cross(w2, r2)),4);

r3 = R_03(q1, q2, q3)*S3;
Vc3 = vpa(simplify(v3 + cross(w3, r3)),4) ;

% Velocity Magnitudes 

Vc2_mag =  norm(Vc2);
Vc3_mag = norm(Vc3);
w1_mag = norm(w1);
w2_mag = norm(w2);
w3_mag = norm(w3);

%% Kinetic & Potential Energy Functions 

k1 = .5*(Izz_1*(w1_mag).^2 ); 
k2 = .5*(Izz_2*(w2_mag).^2 + m2*(Vc2_mag).^2);
k3 = .5*(Izz_3*(w3_mag).^2 + m3*(Vc3_mag).^2);

K = vpa(simplify(k1 + k2 + k3),4);

U = vpa(m2*dot(g_t',r2),4) + vpa(m3*dot(g_t',r3),4) ;


L = K - U; 


dL_dq1 = diff(L,q1)
dL_dq2 = diff(L,q2)
dL_dq3 = diff(L,q3)










%% Rotation Matrices


function out = R_01(theta1)

out =[cos(theta1), -sin(theta1), 0; ...
      sin(theta1), cos(theta1) , 0; ... 
           0     ,      0      , 1]; 
end 


function out = R_02(theta1, theta2)

out =[cos(theta1)*cos(theta2), -cos(theta1)*sin(theta2), -sin(theta1); ...
      sin(theta1)*cos(theta2), -sin(theta1)*sin(theta2),cos(theta1); ... 
      -sin(theta2)         ,      -cos(theta2)      , 0]; 
end 

function out = R_03(theta1, theta2, theta3)

out =[cos(theta1)*cos(theta2+ theta3), -cos(theta1)*sin(theta2 + theta3), -sin(theta1); ...
      sin(theta1)*cos(theta2 + theta3), -sin(theta1)*sin(theta2 + theta3) , cos(theta1); ... 
      -sin(theta2+theta3),      -cos(theta2 + theta3)      , 0]; 
end 













