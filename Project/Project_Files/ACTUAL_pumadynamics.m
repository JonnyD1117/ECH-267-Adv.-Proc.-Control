
L1 = Link('d', 0, 'a', 0, 'alpha', 0);
L2 = Link('d', .2435, 'a', 0, 'alpha', -pi/2);
L3 = Link('d', -.0934, 'a', .4318, 'alpha', 0);
L4 = Link('d',.4331 , 'a', -.0203, 'alpha', pi/2);

clear L
deg = pi/180;

% joint angle limits from 
% A combined optimization method for solving the inverse kinematics problem...
% Wang & Chen
% IEEE Trans. RA 7(4) 1991 pp 489-

% all parameters are in SI units: m, radians, kg, kg.m2, N.m, N.m.s etc.
L(1) = Revolute('d', 0, ...   % link length (Dennavit-Hartenberg notation)
    'a', 0, ...               % link offset (Dennavit-Hartenberg notation)
    'alpha', pi/2, ...        % link twist (Dennavit-Hartenberg notation)
    'I', [0, 0.35, 0, 0, 0, 0], ... % inertia tensor of link with respect to center of mass I = [L_xx, L_yy, L_zz, L_xy, L_yz, L_xz]
    'r', [0, 0, 0], ...       % distance of ith origin to center of mass [x,y,z] in link reference frame
    'm', 0, ...               % mass of link
    'Jm', 200e-6, ...         % actuator inertia 
    'G', -62.6111, ...        % gear ratio
    'B', 1.48e-3, ...         % actuator viscous friction coefficient (measured at the motor)
    'Tc', [0.395 -0.435], ... % actuator Coulomb friction coefficient for direction [-,+] (measured at the motor)
    'qlim', [-160 160]*deg ); % minimum and maximum joint angle

L(2) = Revolute('d', 0, 'a', 0.4318, 'alpha', 0, ...
    'I', [0.13, 0.524, 0.539, 0, 0, 0], ...
    'r', [-0.3638, 0.006, 0.2275], ...
    'm', 17.4, ...
    'Jm', 200e-6, ...
    'G', 107.815, ...
    'B', .817e-3, ...
    'Tc', [0.126 -0.071], ...
    'qlim', [-45 225]*deg );

L(3) = Revolute('d', 0.15005, 'a', 0.0203, 'alpha', -pi/2,  ...
    'I', [0.066, 0.086, 0.0125, 0, 0, 0], ...
    'r', [-0.0203, -0.0141, 0.070], ...
    'm', 4.8, ...
    'Jm', 200e-6, ...
    'G', -53.7063, ...
    'B', 1.38e-3, ...
    'Tc', [0.132, -0.105], ...
    'qlim', [-225 45]*deg );

L(4) = Revolute('d', 0.4318, 'a', 0, 'alpha', pi/2,  ...
    'I', [1.8e-3, 1.3e-3, 1.8e-3, 0, 0, 0], ...
    'r', [0, 0.019, 0], ...
    'm', 0.82, ...
    'Jm', 33e-6, ...
    'G', 76.0364, ...
    'B', 71.2e-6, ...
    'Tc', [11.2e-3, -16.9e-3], ...
    'qlim', [-110 170]*deg);

L(5) = Revolute('d', 0, 'a', 0, 'alpha', -pi/2,  ...
    'I', [0.3e-3, 0.4e-3, 0.3e-3, 0, 0, 0], ...
    'r', [0, 0, 0], ...
    'm', 0.34, ...
    'Jm', 33e-6, ...
    'G', 71.923, ...
    'B', 82.6e-6, ...
    'Tc', [9.26e-3, -14.5e-3], ...
    'qlim', [-100 100]*deg );


L(6) = Revolute('d', 0, 'a', 0, 'alpha', 0,  ...
    'I', [0.15e-3, 0.15e-3, 0.04e-3, 0, 0, 0], ...
    'r', [0, 0, 0.032], ...
    'm', 0.09, ...
    'Jm', 33e-6, ...
    'G', 76.686, ...
    'B', 36.7e-6, ...
    'Tc', [3.96e-3, -10.5e-3], ...
    'qlim', [-266 266]*deg );

% bot = SerialLink([L1 L2 L3 L4], 'name', 'PUMA 560');
% bot = SerialLink([L1 L2 L3], 'name', 'PUMA 560');

% bot = SerialLink([L2 L3 L4], 'name', 'PUMA 560');
% % 
% bot.plot([0 0 0])
% bot.plot([0 0 0])


p560 = SerialLink(L, 'name', 'Puma 560', ...
    'manufacturer', 'Unimation', 'ikine', 'puma', 'comment', 'viscous friction; params of 8/95');


p560.model3d = 'UNIMATE/puma560';


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
    
    p560.plot([q0(:,i)', 0 0 0])
    
end 

% theta_1 = q0(1,:);
% theta_2 = q0(2,:);
% theta_3 = q0(3,:);
% 
% figure()
% hold on 
% plot(theta_1)
% plot(theta_2)
% plot(theta_3)





function out = M(theta)

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);

m11 = (22 +.9*(cos(q2))^2); 
m12 = (1.17 + 1.92*sin(q2)); 
m13 = (-.3*cos(q2 + q3)); 

m21 = (1.17 + 1.92*sin(q2)); 
m22 = 1.66; 
m23 = (-.29);

m31 = (-.3*cos(q2 + q3)); 
m32 = (-.29); 
m33 = (.11);

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


v1 = -1.8*cos(q2)*sin(q2)*q1_dot*q2_dot + 1.92*cos(q2)*(q2_dot)^2 + .3*sin(q2+q3)*q2_dot*q3_dot +.3*sin(q2+q3)*(q3_dot)^2; 
v2 = 1.92*cos(q2)*q1_dot*q2_dot + .9*cos(q2)*sin(q2)*(q1_dot)^2 - 1.92*cos(q2)*q1_dot*q2_dot -.3*sin(q2+q3)*q1_dot*q2_dot; 
v3 = .3*sin(q2+q3)*q1_dot*q2_dot + .3*sin(q2+q3)*q1_dot*q3_dot - .3*sin(q2+q3)*q1_dot*q3_dot; 
out = [v1; v2; v3];

end 

function out = G(theta)

q1 = theta(1);
q2 = theta(2);
q3 = theta(3);


g1 = 0;
g2 = 91.9*sin(q2) - 3.3*sin(q2)*cos(q3) - 3.3*cos(q2)*sin(q3);
g3 = -3.3*sin(q3)*cos(q2) - 3.3*sin(q2)*cos(q3);



out = [g1; g2; g3];

end 

function q_ddot = PUMA_560(q0, q0_dot, tau)


m = M(q0); 
v = V(q0, q0_dot);
g = G(q0);

q_ddot = m\(tau - v - g);
end 
