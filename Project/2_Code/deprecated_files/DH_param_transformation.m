%% 3 Degree of Freedom: DH Parameter Homogeneous Transformation Derivation  
clear all 
clc

syms a0 a1 a2 alpha0 alpha1 alpha2 d1 d2 d3 theta1 theta2 theta3 


% Derivation 

% Frame to Frame Transforms
T_01 =  T_h(0, 0, 0, theta1);
T_12 =  T_h(-pi/2, 0, .2435, theta2);
T_23 =  T_h(0, .4318, -.0934, theta3);
T_34 = T_h(pi/2, -.0203, -.0934, 0);

% Frame to Origin Transforms
T_02 = vpa(simplify(T_01*T_12),4);
T_03 = vpa(simplify(T_01*T_12*T_23),4);
T_04 = vpa(simplify(T_01*T_12*T_23*T_34),4);


% size(T_04)
P_04_x = vpa(simplify(T_04(1,4)),4)
P_04_y = vpa(simplify(T_04(2,4)),4)
P_04_z = vpa(simplify(T_04(3,4)),4)





% T1 = T_h(alpha0, a0, d1, theta1); 
% T2 = T_h(alpha1, a1, d2, theta2); 
% T3 = T_h(alpha2, a2, d3, theta3);
% 
% T_tot = T1*T2*T3;



%%

function output = T_Rz(theta)

R_z = [cos(theta) -sin(theta) 0; 
       sin(theta) cos(theta) 0; 
           0          0      1 ]; 
       
    output = [R_z [0;0;0]; 
              [0 0 0] 1];  
end 

function output = T_z(dist)

output = [1 0 0 0; 
          0 1 0 0; 
          0 0 1 dist;
          0 0 0 1];
end 

function output = T_Rx(theta)

R_x = [ 1 0 0 ;
        0 cos(theta) -sin(theta); 
       0 sin(theta) cos(theta) ]; 
       
    output = [R_x [0;0;0]; 
              [0 0 0] 1];  
end 

function output = T_x(dist)

output = [1 0 0 dist; 
          0 1 0 0; 
          0 0 1 0;
          0 0 0 1];
end 

function output = T_h(alpha, a, d, theta)


% output = T_Rz(theta)*T_z(d)*T_x(a)*T_Rx(alpha);
% output = T_Rz(theta)*T_z(d)*T_Rx(alpha)*T_x(a);

output = T_Rx(alpha)*T_x(a)*T_Rz(theta)*T_z(d);



end 