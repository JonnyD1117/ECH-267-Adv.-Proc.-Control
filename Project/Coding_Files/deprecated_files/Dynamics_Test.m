
syms t q1(t) q2(t) q3(t) L1 L2 L3 g m1 m2 m3 I1 I2 I3 T1 T2 T3


q1_dot = diff(q1,t) ; 
q2_dot = diff(q2,t); 
q3_dot = diff(q3,t); 

q1_ddot = diff(q1_dot,t) ; 
q2_ddot = diff(q2_dot,t); 
q3_ddot = diff(q3_dot,t); 

M = [(I1 + .25*L2^2*m2+L2^2*m3) (.5*L2*L3*m3*cos(q1-q2)) 0; ...  
     (.5*L2*L3*m3*cos(q1-q2))  (I2 + .25*L3^2*m3) 0; ...
     0 0 I3];
 

V = [.5*L2*L3*m3*sin(q1-q2)*(q2_dot)^2; -.5*L2*L3*m3*sin(q1-q2)*(q1_dot)^2 ; 0];

K = [.5*L2*g*m2*cos(q1)+L2*g*m3*cos(q1); .5*L3*m3*cos(q2);0]; 

M_inv = inv(M);





theta_dot_dot = M_inv*([T1;T2;T3]- V -K)

% theta_dot_dot_sub = subs(theta_dot_dot,{T1, T2, T3, L1, L2, L3, I1, I2, I3, m1,m2, m3, g, q1, q2, q3, q1_dot, q2_dot,q3_dot,},{0,0,0,1,1,1,1,1,1,.5,.5,.5,9.81, 0,0,0, 0,0,0})
% 


%% 

params = struct;

params.L1 = .25;
params.L2 = .25;
params.L3 = .25;
params.I1 = 1;
params.I2 = 1;
params.I3 = .5;
params.m1 = 1;
params.m2 = 1;
params.m3 = 1;
params.g = 9.81;


q0 = [0;0;0];
q0_dot = [0;0;0];

u = [0, 0, 0]; 

q_ddot = RobotDynamics(q0, q0_dot, u, params);

i_param = struct;
i_param.dt = .1;


for t=1:1:100
    
    output(:, t) = RobotDynamics_with_integrators(q0, q0_dot, u, params, i_param) ;

    q0 = output(1:3,t);
    q0_dot = output(4:6,t);
end 


positions = output(1:3,:);

q1 = positions(1,:); 
q2 = positions(2,:) ;
q3 = positions(3,:) ; 

figure()

subplot(3,1,1)
plot(q1)
title('Robot Joint Angles')
ylabel('q1')
subplot(3,1,2)
plot(q2)
ylabel('q2')
subplot(3,1,3)
plot(q3)
ylabel('q3')
xlabel('Time Steps')




% output = RobotDynamics_with_integrators(q0, q0_dot, u, params, i_param) 




%% Animate the Robot:


L1 = .25; 
L2 = .5; 
L3 = .5; 


for k = 1:1:100
    
    % Point-Wise Kinematics 
x0 = 0; 
y0 = 0; 
z0 = 0;

x1 = 0;
x2 = L2*cos(q1(k));
x3 = L2*cos(q1(k)) + L3*cos(q2(k));

y1 = L1;
y2 = L1 + L2*sin(q1(k));
y3 = L1 + L2*sin(q1(k)) + L3*sin(q2(k));

z1 = 0; 
z2 = -(L2*cos(q1(k))*cos(q3(k))); 
z3 = -(L2*cos(q1(k))+L3*cos(q2(k)))*cos(q3(k)); 


p1 = [x0,x1,x2,x3]; % List of all X points in Sequence
p2 = [y0,y1, y2,y3]; % List of all Y points in Sequence
p3 = [z0,z1,z2,z3]; % List of all Z points in Sequence
    
    
    
    
    
plot3(p1, p3, p2)
xlabel("X Label")
ylabel("Y Label")
zlabel("Z Label")
xlim([-2,2])
ylim([0,2])
zlim([-2,2])

    pause(0.5);    
    
end 







%%

function accel = RobotDynamics(q_s, q_dots, Taus, params)

L1 = params.L1;
L2 = params.L2;
L3 = params.L3;
I1 = params.I1;
I2 = params.I2;
I3 = params.I3;
m1 = params.m1;
m2 = params.m2;
m3 = params.m3;
g = params.g;

q1 = q_s(1);
q2 = q_s(2);
q3 = q_s(3);

q1_dot = q_dots(1);
q2_dot = q_dots(2);
q3_dot = q_dots(3);

T1 = Taus(1);
T2 = Taus(2);
T3 = Taus(3);

 
 P1 =  - (4*(m3*L3^2 + 4*I2)*((L2*g*m2*cos(q1))/2 - T1 + L2*g*m3*cos(q1) + (L2*L3*m3*sin(q1 - q2)*q2_dot^2)/2))/(16*I1*I2 + 4*L2^2*L3^2*m3^2 + 4*I2*L2^2*m2 + 4*I1*L3^2*m3 + 16*I2*L2^2*m3 + L2^2*L3^2*m2*m3 - 4*L2^2*L3^2*m3^2*cos(q1 - q2)^2) - (8*L2*L3*m3*cos(q1 - q2)*(T2 - (L3*m3*cos(q2))/2 + (L2*L3*m3*sin(q1 - q2)*q1_dot^2)/2))/(16*I1*I2 + 4*L2^2*L3^2*m3^2 + 4*I2*L2^2*m2 + 4*I1*L3^2*m3 + 16*I2*L2^2*m3 + L2^2*L3^2*m2*m3 - 4*L2^2*L3^2*m3^2*cos(q1 - q2)^2);
 
 P2= 4*(T2 - (L3*m3*cos(q2))/2 + (L2*L3*m3*sin(q1 - q2)*q1_dot^2)/2)*(4*I1 + L2^2*m2 + 4*L2^2*m3)/(16*I1*I2 + 4*L2^2*L3^2*m3^2 + 4*I2*L2^2*m2 + 4*I1*L3^2*m3 + 16*I2*L2^2*m3 + L2^2*L3^2*m2*m3 - 4*L2^2*L3^2*m3^2*cos(q1 - q2)^2) + (8*L2*L3*m3*cos(q1 - q2)*((L2*g*m2*cos(q1))/2 - T1 + L2*g*m3*cos(q1) + (L2*L3*m3*sin(q1 - q2)*q2_dot^2)/2))/(16*I1*I2 + 4*L2^2*L3^2*m3^2 + 4*I2*L2^2*m2 + 4*I1*L3^2*m3 + 16*I2*L2^2*m3 + L2^2*L3^2*m2*m3 - 4*L2^2*L3^2*m3^2*cos(q1 - q2)^2);
 
 P3 = T3/I3;
 
 accel = [ P1 ; P2 ; P3]; 
 
end 



function out = RobotDynamics_with_integrators(q_s, q_dots, Taus, robot_params, int_params) 

% Define Parameters
dt = int_params.dt;

% Set Initial Values
q0 = q_s;
q0_dot = q_dots;
q0_ddot = RobotDynamics(q_s, q_dots, Taus, robot_params);


% Perform Euler-Integration
q_dot_next = q0_dot + q0_ddot*dt;
q_next = q0 + q0_dot*dt + .5*q0_ddot*(dt)^2;


pos = q_next; 
vel =q_dot_next;
accel = q0_ddot;

out = [pos; vel; accel;];

end 
   