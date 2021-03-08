syms q1 q2 q3 q1_dot q2_dot q3_dot q1_ddot q2_ddot q3_ddot m1 m2 m3


q0 = [0;0; 0];
q0_dot = [0;0;0];

tau = [1;0;0];
dt = .1;

p = struct ;

p.m1 = .5;
p.m2 = .25; 
p.L1 = 1.4;
p.L2 = 1.2; 
p.g = 9.81;

   q0_ddot(:,1) = PUMA_560( q0(:,1), q0_dot(:,1), tau);
%  q0_ddot(:,1) = Robot_Stuff( q0(:,1), q0_dot(:,1), tau);
%  q0_ddot(:,1) = PUMA_560( q0(:,1), q0_dot(:,1), tau);
%  q0_ddot(:,1) = RobotModel_q_ddot(p, q0(:,1), q0_dot(:,1), tau);


 pos =  forward_kinematics(q0(:,1), q0_dot(:,1));
    
    x_vec = pos.x;
    y_vec = pos.y;
    z_vec = pos.z;

 plot_robot(x_vec, y_vec, z_vec)






%%





for i=2:1:100
    

     q0_ddot(:,i) = PUMA_560( q0(:,i-1), q0_dot(:,i-1), tau);
%    q0_ddot(:,i-1) = Robot_Stuff( q0(:,i-1), q0_dot(:,i-1), tau);
%    q0_ddot(:,i) = RobotModel_q_ddot(p, q0(:,i-1), q0_dot(:,i-1), tau);
%    q0_ddot(:,i) = [.1;0;0];
    
    q0_ddot(3,i) = 0; 
    q0_dot(:,i) = q0_dot(:,i-1) + dt*q0_ddot(:,i); 
    q0(:,i) = q0(:,i-1) + dt*q0_dot(:,i) + .5*(dt)^2*(q0_ddot(:,i));
    
    
    pos =  forward_kinematics(q0(:,i), q0_dot(:,i));
    
    x_vec = pos.x;
    y_vec = pos.y;
    z_vec = pos.z;

    plot_robot(x_vec, y_vec, z_vec)
    
%     bot.plot([q0(:,i)' 0])
    % p560.plot([q0(:,i)'0])
    % bot.plot([q0(:,i)'])
end 






%%


function q_ddot = RobotModel_q_ddot(p, theta, theta_dot, tau_input)

q1 = theta(1);
q2 = theta(2); 
q3 = theta(3); 

q1_dot = theta_dot(1);
q2_dot = theta_dot(2); 
q3_dot = theta_dot(3); 

tau1 = tau_input(1);
tau2 = tau_input(2);
tau3 = tau_input(3);

tau = [tau1;tau2];

I_3 = 1;


m1 = p.m1; 
m2 = p.m2; 
L1 = p.L1; 
L2 = p.L2; 
g_val = p.g; 

eps = .1 ;


m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(q2) + (L2)^2) + eps, m2*(L1*L2*cos(q2)+ (L2)^2); ...
     m2*(L1*L2*cos(q2) + (L2)^2), m2*(L2)^2+eps];

v = [-m2*L1*L2*sin(q2)*(2*q1_dot*q2_dot + (q2_dot)^2);...
     m2*L1*L2*(q1_dot)^2*sin(q2)] ; 

g = [(m1 + m2)*L1*g_val*cos(q1)+ m2*g_val*L2*cos(q1+ q2); ...
     m2*g_val*L2*cos(q1+q2)];
 

 f_const= .00001;
f = [q1_dot*f_const;q2_dot*f_const];
 
 q_ddot = m\(tau -v -g - f); 
 
I_3 = .25+.524+.125;


q3_ddot = (1/I_3)*(tau3) ;

q_ddot = [q3_ddot;q_ddot];

end 

function q_ddot = Robot_Stuff(theta, theta_dot, tau)

m1 = .2; 
m2 = .2; 
m3 = .15; 
m4 = .5; 

d1 = .05; 
a2 = .140; 
a3 = .120; 

eps = .1;

g = 9.81;

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);

q1_dot = theta_dot(1);
q2_dot = theta_dot(2); 
q3_dot = theta_dot(3); 

M11 = eps + (m4+m3+.25*m2)*(a2)^2*(cos(q2))^2 + (m4+.25*m3)*(a3)^2*(sin(q2+q3))^2 + (m4+m3)*a2*a3*cos(q2)*sin(q2+q3)+1;
M12 = 0; 
M13 = 0; 

M21 = 0;
M22 = eps + (m4+m3+.25*m2)*(a2)^2 + (m4+.25*m3)*(a3)^2 + (2*m4+m3)*a2*a3*sin(q2) + 1;
M23 = (m4+m3+.25*m2)*(a3)^2 + (m4+.5*m3)*a2*a3*sin(q2) +1;

M31 = 0 ;
M32 = M23;
M33 = eps + (m4+.25*m3)*(a3)^2 +1 ;

C1 = q1_dot*q2_dot*(m4+.25*m3)*(a3)^2*sin(2*(q2+q3)) -(m4+m3+.25*m2)*(a2)^2*sin(2*q2)+ q1_dot*q3_dot*((m4+.25*m3)*(a3)^2*sin(2*(q2+q3)) + (m4+m3)*a2*a3*cos(q2)*cos(q2+q3)) + (m4+m3)*a2*a3*cos(2*(q2+q3));
C2 =  q2_dot*q3_dot*(-(2*m4+m3)*a2*a3*cos(q3)) + (q3_dot)^2*(-(m4+.5*m3)*a2*a3*cos(q3)) + (m4+.25*m3)*(a3)^2*sin(2*(q2+q3)) + (m4+.5*m3)*a2*a3*cos(2*(q2+q3)) + .5*(q1_dot)^2*((m4+m3+.25*m2)*(a2)^2*sin(2*q2));
C3 = .5*(q1_dot)^2*((m4 + .25*m3)*(a3)^2*sin(2*(q2+q3)) + (m4+m3)*a2*a3*cos(q2)*cos(q2+q3)) -(.5*(q2_dot)^2 + q2_dot*q3_dot)*(2*m4 + m3)*a2*a3*cos(q3);

G1 = 0;
G2 = (m4+m3+.5*m2)*g*a2*cos(q2) + (m4+.5*m3)*g*a3*sin(q2+q3);
G3 = (m4+.5*m3)*g*a3*sin(q2+q3);

M_m = [M11, M12, M13;...
       M21, M22, M23;... 
       M31, M32, M33]; 
   
C_m = [C1; C2 ;C3];

G_m = [G1; G2; G3];


F = .5*theta_dot; 


% q_ddot = M_m\(tau - G_m - C_m -F);
q_ddot = M_m\(tau - G_m - C_m -F);



end 

function plot_robot(x_vec, y_vec, z_vec)

x1 = x_vec(1); 
x2 = x_vec(2); 
x3 = x_vec(3);
x4 = x_vec(4);

y1 = y_vec(1); 
y2 = y_vec(2); 
y3 = y_vec(3);
y4 = y_vec(4);

z1 = z_vec(1); 
z2 = z_vec(2); 
z3 = z_vec(3);
z4 = z_vec(4);

hold on
% fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background

plot3([ 0, x1 ], [0, y1], [0, z1], "b");
plot3([ x1, x2], [y1, y2], [z1, z2], "r");
plot3([ x2, x3], [y2, y3], [z2, z3], "g");
plot3([ x3, x4], [y3, y4], [z3, z4], "k");

hold off

pause(.2)


title('Planar 2D Robot in Workspace')
xlim([-2,2])
ylim([-2,2])
zlim([-2,2])

xlabel("X-Axis")
ylabel("Y-Axis") 
zlabel("Z-Axis") 

az = 45;
el = 45;
view(az, el);


end

function pos = forward_kinematics(theta, theta_dot)
 
L1 = .5; 
L2 = 1.4; 
L3 = 1.2; 

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);

q1_dot = theta_dot(1); 
q2_dot = theta_dot(2); 
q3_dot = theta_dot(3); 


% Frame to Frame Transforms
T_01 =  T_h(0, 0, .125, q1);
T_12 =  T_h(-pi/2, 0, .2435, q2);
T_23 =  T_h(0, .4318, -.0934, q3);
% T_34 = T_h(pi/2, -.0203, .4331,0);
T_34 = T_h(pi/2, -.0203, .4331,0);


% % Frame to Origin Transforms
T_02 = vpa((T_01*T_12),4);
T_03 = vpa((T_01*T_12*T_23),4);
T_04 = vpa((T_01*T_12*T_23*T_34),4);

P_01_x = vpa((T_01(1,4)),4);
P_01_y = vpa((T_01(2,4)),4);
P_01_z = vpa((T_01(3,4)),4);

P_02_x = vpa((T_02(1,4)),4);
P_02_y = vpa((T_02(2,4)),4);
P_02_z = vpa((T_02(3,4)),4);

P_03_x = vpa((T_03(1,4)),4);
P_03_y = vpa((T_03(2,4)),4);
P_03_z = vpa((T_03(3,4)),4);

P_04_x = vpa((T_04(1,4)),4);
P_04_y = vpa((T_04(2,4)),4);
P_04_z = vpa((T_04(3,4)),4);

x = [P_01_x;P_02_x;P_03_x; P_04_x]; 
y = [P_01_y;P_02_y;P_03_y;P_04_y]; 
z = [P_01_z;P_02_z;P_03_z;P_04_z];



% P_01_x = vpa((T_01(1,4)),4);
% P_01_y = vpa((T_01(2,4)),4);
% P_01_z = vpa((T_01(3,4)),4);
% 
% P_12_x = vpa((T_12(1,4)),4);
% P_12_y = vpa((T_12(2,4)),4);
% P_12_z = vpa((T_12(3,4)),4);
% 
% P_23_x = vpa((T_23(1,4)),4);
% P_23_y = vpa((T_23(2,4)),4);
% P_23_z = vpa((T_23(3,4)),4);
% 
% P_34_x = vpa((T_34(1,4)),4);
% P_34_y = vpa((T_34(2,4)),4);
% P_34_z = vpa((T_34(3,4)),4);
% 
% x = [P_01_x;P_12_x;P_23_x; P_34_x]; 
% y = [P_01_y;P_12_y;P_23_y;P_34_y]; 
% z = [P_01_z;P_12_z;P_23_z;P_34_z];

% x1 = 0;
% x2 = L2*cos(q1)*cos(q3); 
% x3 = (L2*cos(q1) + L3*cos(q1+q2))*cos(q3); 
% 
% y1 = 0;
% y2 = L2*cos(q1)*sin(q3); 
% y3 = (L2*sin(q1) + L3*sin(q1+q2))*sin(q3); 
% 
% z1 = L1;
% z2 = L1 + L2*sin(q1); 
% z3 = L1 + L2*sin(q1) + L3*sin(q2+q1);

% x = [x1;x2;x3]; 
% y = [y1;y2;y3]; 
% z = [z1;z2;z3];

pos = struct ;

pos.x = x; 
pos.y = y; 
pos.z = z;

end 

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





function out = M(theta)

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);

eps = .1;

m11 =  (22+.9*cos(q2)^2) + eps; 
m12 =  (1.92*sin(q2)+1.17);
m13 =  -.3*cos(q2+q3);

m21 =  (1.92*sin(q2)+1.17);
m22 =  (1.66)+ eps;
m23 =  (-.29);

m31 = -.3*cos(q2+q3);
m32 = (-.29);
m33 = (.11+eps);

out = [m11, m12, m13; 
       m21, m22, m23;
       m31, m32, m33];

end 

function out = V(theta, theta_dot)

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);

q1_dot = theta_dot(1);
q2_dot = theta_dot(2);
q3_dot = theta_dot(3);


v1 = 1.92*cos(q2)*(q2_dot)^2 -1.8*cos(q2)*sin(q2)*q1_dot*q2_dot + .3*sin(q2+q3)*q2_dot*q3_dot + .3*sin(q2+q3)*(q3_dot)^2;  
v2 = 1.92*cos(q2)*q1_dot*q2_dot + .9*cos(q2)*sin(q2)*(q1_dot)^2 -1.92*cos(q2)*q1_dot*q2_dot -.3*sin(q2+q3)*q1_dot*q3_dot;
v3 = .3*sin(q2+q3)*q1_dot*q2_dot + .3*sin(q2+q3)*q1_dot*q3_dot -.3*sin(q2+q3)*q1_dot*q3_dot; 
out = [v1; v2; v3];

end 

function out = G(theta)

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);

g1 = 0;
g2 = 91.9*sin(q2) -3.3*sin(q2)*cos(q3) -3.3*cos(q2)*sin(q3);
g3 = -3.3*sin(q3)*cos(q2) -3.3*sin(q2)*cos(q3);

out = [g1; g2; g3];
end 

function q_ddot = PUMA_560(q0, q0_dot, tau)
    m = M(q0); 
    v = V(q0, q0_dot);
    g = G(q0);

    f = .5*q0_dot;
    q_ddot = m\(tau - v - g -f);
end 

