
L1 = Link('a', 0, 'd', 0, 'alpha', 0);
L2 = Link('a', 0, 'd', .2435, 'alpha', -pi/2);
L3 = Link('a', .4318, 'd', -.0934, 'alpha',0);
% L4 = Link('a',.4331 , 'd',-.0203 , 'alpha', pi/2);


bot = SerialLink([L1 L2 L3], 'name', 'PUMA 560');

bot.plot([0 0 0 ])
%%
tau = [0;0;0];
dt =.1;

q0 = [0;0;.1];
q0_dot = [0;0;0];
q0_ddot = PUMA_560(q0, q0_dot, tau);


for i=2:1:100
    
    q0_ddot(:,i) = PUMA_560(q0(:,i-1), q0_dot(:,i-1), tau);
    q0_dot(:,i) = q0_dot(:,i-1) + dt*q0_ddot(:,i);
    q0(:,i) = q0(:,i-1) + dt*q0_dot(:,i);
    
%     p560.plot([q0(:,i)'0])
    bot.plot([q0(:,i)'])

    
end 


function out = M(theta)

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);

eps = .1;

m11 =  (22+.9*cos(q2)^2) + eps; 
m12 =  (1.92*sin(q2)+1.17);
m13 =  0;

m21 =  (1.92*sin(q2)+1.17);
m22 =  (1.66)+ eps;
m23 =  (-.29);

m31 = 0;
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


v1 = 1.92*cos(q2)*(q2_dot)^2 -1.8*cos(q2)*sin(q2)*q1_dot*q2_dot;  
v2 = .45*sin(2*q2)*(q1_dot)^2;
v3 = 0; 
out = [v1; v2; v3];

end 

function out = G(theta)

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);


g1 = 0;
g2 = 31.9*sin(q2) -3.6*sin(q2+q3);
g3 = -3.6*sin(q2+q3);



out = [g1; g2; g3];

end 

function q_ddot = PUMA_560(q0, q0_dot, tau)
    m = M(q0); 
    v = V(q0, q0_dot);
    g = G(q0);

    f = .25*q0_dot;
    q_ddot = m\(tau - v - g -f);
end 
