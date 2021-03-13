
L1 = Link('a', 0, 'd', 0, 'alpha', 0);
L2 = Link('a', 0, 'd', .2435, 'alpha', -pi/2);
L3 = Link('a', .4318, 'd', -.0934, 'alpha',0);
L4 = Link('d',.4331 , 'a',-.0203 , 'alpha', pi/2);


% bot = SerialLink([L1 L2 L3], 'name', 'PUMA 560');
bot = SerialLink([L1 L2 L3 L4], 'name', 'PUMA 560');


bot.plot([0 0 0 0 ])
%%
tau = [0;0;0];
dt =.1;

q0 = [0;0;pi/4];
q0_dot = [0;0;0];
q0_ddot = PUMA_560(q0, q0_dot, tau);


pos =  forward_kinematics(q0(:,1), q0_dot(:,1));
    
    x_vec = pos.x;
    y_vec = pos.y;
    z_vec = pos.z;

    plot_robot(x_vec, y_vec, z_vec)



%%


for i=2:1:100
    
        q0_ddot(:,i) = PUMA_560(q0(:,i-1), q0_dot(:,i-1), tau);

%     q0_ddot(:,i) = PUMA_560(q0(:,i-1), q0_dot(:,i-1), tau);
    q0_dot(:,i) = q0_dot(:,i-1) + dt*q0_ddot(:,i);
    q0(:,i) = q0(:,i-1) + dt*q0_dot(:,i);
    
    
    pos =  forward_kinematics(q0(:,i), q0_dot(:,i));
    
    x_vec = pos.x;
    y_vec = pos.y;
    z_vec = pos.z;

    plot_robot(x_vec, y_vec, z_vec)
    
    
    %p560.plot([q0(:,i)'0])
%     bot.plot([q0(:,i)'])
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
 
L1 = .25 
L2 = .35; 
L3 = .35;

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);

q1_dot = theta_dot(1); 
q2_dot = theta_dot(2); 
q3_dot = theta_dot(3); 


x1 = 0;
x2 = L1*cos(q2)*sin(q1); 
x3 = (L2*cos(q2) + L3*cos(q2+q3))*sin(q1); 

y1 = 0;
y2 = L1*cos(q2)*cos(q1); 
y3 = (L2*cos(q2) + L3*cos(q2+q3))*cos(q1); 

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




