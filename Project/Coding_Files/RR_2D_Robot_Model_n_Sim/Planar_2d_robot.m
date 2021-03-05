
clear all 

close all 
p = struct ;

p.m1 = .5;
p.m2 = .25; 
p.L1 = .25;
p.L2 = .25; 
p.g = 9.81;

q0 = [deg2rad(0);deg2rad(25)];
q0_dot = [0;0];

dt = .1;

t1 = 2*ones(1,150);
t2 = zeros(1,150);

% tau = [ [t1;t2] ,zeros(2,200)];
tau = [2;0];



q0_ddot = Planar2D_Robot(p,q0, q0_dot, [0;0]);

for i=2:1:350
    figure(1);

    q0_ddot(:,i) = Planar2D_Robot(p,q0(:,i-1), q0_dot(:,i-1), tau);
    q0_dot(:,i) = q0_dot(:,i-1) + dt*q0_ddot(:,i); 
    q0(:,i) = q0(:,i-1) + dt*q0_dot(:,i) + .5*(dt)^2*(q0_ddot(:,i));
    
    pos = forward_kinematics(p, q0(:,i), q0_dot(:,i));
    
    x = pos.x; 
    y = pos.y; 
    
  
    plot_robot(x,y)
   
end 


theta_1 = q0(1,:);
theta_2 = q0(2,:);

% plot(theta_1)

% plot_robot([1,2],[1,3])

function plot_robot(x_vec, y_vec)

x1 = x_vec(1); 
x2 = x_vec(2); 
y1 = y_vec(1); 
y2 = y_vec(2); 



hold on
fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background

plot([ 0, x1], [0, y1], "b");
plot([ x1, x2], [y1, y2], "r");
hold off

title('Planar 2D Robot in Workspace')
xlim([-.6,.6])
ylim([-.6,.6])

xlabel("X-Axis")
ylabel("Y-Axis") 



end

function pos = forward_kinematics(p,theta, theta_dot)
 
L1 = p.L1; 
L2 = p.L2; 
g = p.g; 

q1 = theta(1);
q2 = theta(2);

q1_dot = theta_dot(1); 
q2_dot = theta_dot(2); 


x = [L1*cos(q1) ;L1*cos(q1) + L2*cos(q1 + q2)];
y = [L1*sin(q1) ;L1*sin(q1) + L2*sin(q1 + q2)];

pos = struct ;

pos.x = x; 
pos.y = y; 

end 

function q_ddot = Planar2D_Robot(p, theta, theta_dot, tau)

m1 = p.m1; 
m2 = p.m2; 
L1 = p.L1; 
L2 = p.L2; 
g = p.g; 

eps = .1 ;

m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(theta(2)) + (L2)^2) + eps, m2*(L1*L2*cos(theta(2))+ (L2)^2); ...
     m2*(L1*L2*cos(theta(2)) + (L2)^2), m2*(L2)^2+eps];

v = [-m2*L1*L2*sin(theta(2))*(2*theta_dot(1)*theta_dot(2) + (theta_dot(2))^2);...
     m2*L1*L2*(theta_dot(1))^2*sin(theta(2))] ; 

g = [(m1 + m2)*L1*g*cos(theta(1))+ m2*g*L2*cos(theta(1)+ theta(2)); ...
     m2*g*L2*cos(theta(1)+theta(2))];
 

f = [theta_dot(1)*.25;theta_dot(2)*.25];
 
 q_ddot = m\(tau - v -g - f); 


end 

