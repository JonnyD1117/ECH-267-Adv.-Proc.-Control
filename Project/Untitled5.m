
syms q1 q2 q3 q1_dot q2_dot q3_dot q1_ddot q2_ddot q3_ddot m1 m2 m3


q0 = [0;0;0];
q0_dot = [0;0;0];

tau = [0;0;0];
dt = .1;


L1 = Link('a', 0, 'd', 0, 'alpha', 0);
L2 = Link('a', 0, 'd', .2435, 'alpha', -pi/2);
L3 = Link('a', .4318, 'd', -.0934, 'alpha',0);
L4 = Link('a',.4331 , 'd',-.0203 , 'alpha', pi/2);


% bot = SerialLink([L1 L2 L3], 'name', 'PUMA 560');
bot = SerialLink([L1 L2 L3 L4 ], 'name', 'PUMA 560');


bot.plot([0 0 0 0])


%%





for i=2:1:100
    
    q0_ddot(:,i) = Robot_Stuff(q0(:,i-1), q0_dot(:,i-1), tau);
    q0_dot(:,i) = q0_dot(:,i-1) + dt*q0_ddot(:,i);
    q0(:,i) = q0(:,i-1) + dt*q0_dot(:,i);
    
    
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

y1 = y_vec(1); 
y2 = y_vec(2); 
y3 = y_vec(3);

z1 = z_vec(1); 
z2 = z_vec(2); 
z3 = z_vec(3);




hold on
fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background

plot3([ 0, x1 ], [0, y1], [0, z1], "b");
plot3([ x1, x2], [y1, y2], [z1, z2], "r");
plot3([ x2, x3], [y2, y3], [z2, z3], "g");
hold off

pause(.2)


title('Planar 2D Robot in Workspace')
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])

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



x1 = 0;
x2 = L1*cos(q2)*cos(q1); 
x3 = (L2*cos(q2) + L3*cos(q2+q3))*cos(q1); 

y1 = 0;
y2 = L1*cos(q2)*sin(q1); 
y3 = (L2*cos(q2) + L3*cos(q2+q3))*sin(q1); 

z1 = L1;
z2 = L1 + L2*sin(q2); 
z3 = L1 + L2*sin(q2) + L3*sin(q2+q3);

x = [x1;x2;x3]; 
y = [y1;y2;y3]; 
z = [z1;z2;z3];

pos = struct ;

pos.x = x; 
pos.y = y; 
pos.z = z;

end 



