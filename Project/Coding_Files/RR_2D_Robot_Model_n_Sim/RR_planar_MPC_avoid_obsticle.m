%% Planar Robot 2-Link Multiple Shooting MPC 

import casadi.*

clear all 
close all 
clc 


p = struct ;

p.m1 = .5;
p.m2 = .25; 
p.L1 = .25;
p.L2 = .25; 
p.g = 9.81;

N = 5;
T = .2;

% Control Input Saturation Limits
tau_1_max = 10; % N*m
tau_1_min = -10; %N*m

tau_2_max = 10; % N*m
tau_2_min = -10; %N*m

% Define Symbolic Variables
x1 = SX.sym('x1'); x2 = SX.sym('x2'); x3 = SX.sym('x3'); x4 = SX.sym('x4'); tau1 = SX.sym('tau1'); tau2 = SX.sym('tau2'); 

states = [x1 ; x2 ; x3 ; x4]; n_states = length(states);
ctrl_input = [tau1 ; tau2]; n_ctrl = length(ctrl_input); 

% Define Symbolic State Space RHS
rhs = [x3; x4; RobotModel_q1_ddot(p, x1, x2, x3, x4,tau1, tau2);RobotModel_q2_ddot(p, x1, x2, x3, x4,tau1, tau2)]; 

% Define Symbolic CasADi function for state dynamics 
f = Function('f', {states, ctrl_input}, {rhs}, {'state','input'}, {'dynamics'}); % This function enables us to numerically evaluate the symbolic graph given numerical inputs

% Define Symbolic Matrix for Inputs & Params for each step of Time Horizon 
U = SX.sym('U', n_ctrl, N); % Inputs at each step K in horizon N 
P = SX.sym('P', 2*n_states);% Params at each step K in horizon N, P = (x1_0, x2_0, x1_ref, x2_ref)
X = SX.sym('X', n_states, (N+1)); % States at each step K in horizon N

% Using this symbolic representation of the State Update process, create a function to implement this behavior 
ff = Function('ff',{U,P},{X}, {'Input','Reference '},{'New State'}); % This enables evaluation of symbolic graph give numerical inputs

% Define Objective & Constraints Vector 
obj = 0; % Initialize the Objective Func. to zero
g = [] ; % Initialize the Constraints as Empty 

% Define Loss function Weights on the States and Control Inputs
Q = zeros(4,4); Q(1,1)=10; Q(2,2) =10 ;Q(3,3) =4 ;Q(4,4) =4 ; % Weighting matrices (states)
R = .125*eye(2,2); % Weighting matrices (controls)
Term_cost = zeros(4,4); Term_cost(1,1)=10; Term_cost(2,2) =10 ;Term_cost(3,3) =1 ;Term_cost(4,4) =1 ; % Weighting matrices (states)

st = X(:,1); 
g = [g; (st-P(1:4))]; % Initial Constraint via Defined Initial Condition



% Compute Objective function from Stage Cost 
for k=1:N
    % Symbolically Computes Graph for the objective function
    st = X(:,k); con = U(:,k); ref = P(5:8); 

    obj = obj + (st - ref)'*Q*(st - ref) + con'*R*con; % Construct Symbolic representation of objective function at time step "K"   
  
    st_next = X(:,(k+1));
    f_value = f(st,con); 
    st_next_euler = st + (T*f_value); 
    
    if k < N 
      g = [g; st_next-st_next_euler];
      
    elseif k == N 
      g = [g; st_next-st_next_euler; st_next - ref];

    end 
    
    
end 


obj = obj + (st_next - ref)'*Term_cost*(st_next - ref); % Construct Symbolic representation of objective function at time step "K"   


obs_x = 0.25; % meters
obs_y = 0.25; % meters
obs_diam = 0.3; % meters
for k = 1:N+1   % box constraints due to the map margins
    g = [g ; -sqrt((X(1,k)-obs_x)^2+(X(2,k)-obs_y)^2) + (p.L1/2 + obs_diam/2)];
end

% Define Nonlinear Programming Structure 
OPT_variables = [reshape(X, 4*(N+1), 1) ;reshape(U, 2*(N),1)];  % Via Single Shooting only the Control Input set is used as the decision variable

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
opts = struct; 

% Define Solver Options 
opts.ipopt.max_iter =200; 
opts.ipopt.print_level=0; % 0 - 3
opts.print_time  = 0 ; 
opts.ipopt.acceptable_tol = 1e-6; 
opts.ipopt.acceptable_obj_change_tol = 1e-6; 

% Enstantiate nlpsol class using ipopt 
solver = nlpsol('solver','ipopt', nlp_prob, opts); 

% Define Problem Arguments 
args = struct;

% Equality Constraints, for Entire Trajectory over Horizon
% args.lbg(1:4*(N+1)) = 0; 
% args.ubg(1:4*(N+1)) = 0; 

% Model Prediction Constraint
args.lbg(1:4*(N+1)) = 0; 
args.ubg(1:4*(N+1)) = 0; 
% Terminal State Constraint
args.lbg(4*(N+1):4*(N+1)+4) = 0; 
args.ubg(4*(N+1):4*(N+1)+4) = 0; 
% Object Avoidance Constraint
args.lbg(3*(N+1)+1 : 3*(N+1)+ (N+1)) = -inf; % inequality constraints
args.ubg(3*(N+1)+1 : 3*(N+1)+ (N+1)) = 0; % inequality constraints


%Combined State & Input Constraints (due to multishooting formulation)  
args.lbx(1:4:4*(N+1),1) = -inf; % THETA1 Lower Bound
args.ubx(1:4:4*(N+1),1) = inf;  % THETA1 Upper Bound

args.lbx(2:4:4*(N+1),1) = -inf; % THETA2 Lower Bound
args.ubx(2:4:4*(N+1),1) = inf;  % THETA2 Upper Bound

args.lbx(3:4:4*(N+1),1) = -inf; % THETA1_dot Lower Bound
args.ubx(3:4:4*(N+1),1) = inf;  % THETA1_dot Upper Bound

args.lbx(4:4:4*(N+1),1) = -inf; % THETA2_dot Lower Bound
args.ubx(4:4:4*(N+1),1) = inf;  % THETA2_dot Upper Bound


args.lbx(4*(N+1):2: 4*(N+1) + 2*N,1) = tau_1_min;   % Input Lower Bound TAU1
args.ubx(4*(N+1):2: 4*(N+1) + 2*N,1) = tau_1_max;   % Input Upper Bound TAU1

args.lbx(4*(N+1)+1:2: 4*(N+1) + 2*N,1) = tau_2_min;   % Input Lower Bound TAU2
args.ubx(4*(N+1)+1:2: 4*(N+1) + 2*N,1) = tau_2_max;   % Input Upper Bound TAU2


% SIMULATION LOOP 
t0 = 0; 
x0 = [deg2rad(25); deg2rad(45); 0;0] ;
x_ref = [pi;-pi/4; 0; 0] ;

x_list(:,1) = x0; 
t(1) = t0; 
u0 = zeros(2,N);
X0 = repmat(x0,1, N+1)';
sim_time = 20; 

% Start MPC 
mpc_iter = 0; 
xxl = [ ]; 
u_cl = [ ];


% Main Loop 
main_loop = tic;

while(norm((x0 - x_ref),2) > 1e-4 && mpc_iter < sim_time/T)
    
    args.p = [x0;x_ref];
    args.x0 =  [reshape(X0', 4*(N+1),1); reshape(u0,2*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p); 
    
    u = reshape(full(sol.x(4*(N+1)+1:end))',2,N)'; 
    xxl(:, 1:4, mpc_iter +1) = reshape(full(sol.x(1:4*(N+1)))',4,N+1)';
%     ff_value = ff(u, args.p); 
%     xxl(:, 1:2, mpc_iter +1) = full(ff_value)';
%     u_cl = [u_cl; u(:)]; 
    u_cl = [u_cl; u(1,:)];
    
    t(mpc_iter +1) = t0; 
    [t0, x0, u0] = shift(T, t0, x0, u, f);
    x_list(:,mpc_iter+2) = x0; 
    mpc_iter = mpc_iter +1 ;     
    
end 

u_cl = u_cl';
main_loop_time = toc(main_loop)

figure(1)
subplot(411) 
plot(x_list(1,:));
ylabel('Theta [rads]'); 

subplot(412) 
plot(x_list(2,:));
ylabel('Theta [rads]');

subplot(413)
plot(x_list(3,:)); 
ylabel('Angular Velocity [rad/s]'); 
xlabel('Time [seconds]');

subplot(414)
plot(x_list(4,:)); 
ylabel('Angular Velocity [rad/s]'); 
xlabel('Time [seconds]');


figure(2)
grid on
subplot(211)
stairs(t,u_cl(1,:),'k','linewidth',1.5)
ylabel('Torque (N*m)')

subplot(212)
stairs(t,u_cl(2,:),'k','linewidth',1.5); axis([0 t(end) -0.35 0.75])
xlabel('time (seconds)')
ylabel('Torque (N*m)')
disp("DONE");


figure(3)
for i=1:1:(sim_time/T)
    
    
    theta = x_list(1:2,i); 
    theta_dot = x_list(3:4,i); 
    pos = forward_kinematics(p,theta, theta_dot);
    
    x_vec = pos.x;
    y_vec = pos.y;
    
    plot_robot(x_vec, y_vec);
    pause(.1)
end 

disp(rad2deg(x_list(1,end)))
close all


%% End Effector Trajectory 

x1_traj = []; 
y1_traj = [];

x2_traj = []; 
y2_traj = [];


for i=1:1:100
    
    
    theta = x_list(1:2,i); 
    theta_dot = x_list(3:4,i); 
    pos = forward_kinematics(p,theta, theta_dot);
    
    x_vec = pos.x;
    y_vec = pos.y;
    
    
    x1_traj = [x1_traj,x_vec(1)];
    y1_traj = [y1_traj,y_vec(1)];
    
    x2_traj = [x2_traj,x_vec(2)];
    y2_traj = [y2_traj,y_vec(2)];

    
   
end 

figure(4)
hold on
plot(x1_traj, y1_traj, 'b')
plot(x2_traj, y2_traj,'r')
title("Joint Trajectories")
xlabel("X-Axis")
ylabel("Y-Axis")
legend(["Joint 1"," Joint 2"])


%%

function [t0, x0, u0] = shift(T, t0, x0, u, f)

st = x0; 
con = u(1,:)'; 
f_value = f(st,con); 
st = st+ (T*f_value); 
x0 = full(st); 
t0 = t0 + T;
u0 = [u(2:size(u,1),:); u(size(u,1),:) ];  


end


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




% function Animate_pendulum(x_traj, dt, sim_time)
% 
% v = VideoWriter('pendulum.avi'); 
% v.FrameRate = 30;
% open(v); 
% 
% x1 = x_traj(1,:); 
% 
% sim_len = (sim_time / dt);
% 
% disp(sim_len)
% figure(3)
%     for k=2:length(x1)
%         %Plot for Video 
%         hold  on
%         fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background
%         plot([0 sin(x1(k-1))],[0 -cos(x1(k-1))], 'b', 'LineWidth', 3); % Plots Rod
%     %     plot(,'Marker','o', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); % Plots Bob 
% 
%         xlim([-1-0.2*1 1+0.2*1]); 
%         ylim([-1-0.2*1 1+0.2*1]);
%         title('MPC Inverted Pendulum using CasADi');
%         frame=getframe(gcf); 
%         writeVideo(v, frame);
%     end 
% end 





function q1_ddot = RobotModel_q1_ddot(p, q1, q2, q1_dot, q2_dot,tau1, tau2)

tau = [tau1;tau2] ;

m1 = p.m1; 
m2 = p.m2; 
L1 = p.L1; 
L2 = p.L2; 
g = p.g; 

eps = .1 ;

m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(q2) + (L2)^2) + eps, m2*(L1*L2*cos(q2)+ (L2)^2); ...
     m2*(L1*L2*cos(q2) + (L2)^2), m2*(L2)^2+eps];

v = [-m2*L1*L2*sin(q2)*(2*q1_dot*q2_dot + (q2_dot)^2);...
     m2*L1*L2*(q1_dot)^2*sin(q2)] ; 

g = [(m1 + m2)*L1*g*cos(q1)+ m2*g*L2*cos(q1+ q2); ...
     m2*g*L2*cos(q1+q2)];
 

f = [q1_dot*.25;q2_dot*.25];
 
 q_ddot = m\(tau - v -g - f); 
 
 q1_ddot = q_ddot(1);


end 

function q2_ddot = RobotModel_q2_ddot(p, q1, q2, q1_dot, q2_dot,tau1, tau2)

tau = [tau1;tau2] ;

m1 = p.m1; 
m2 = p.m2; 
L1 = p.L1; 
L2 = p.L2; 
g = p.g; 

eps = .1 ;

m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(q2) + (L2)^2) + eps, m2*(L1*L2*cos(q2)+ (L2)^2); ...
     m2*(L1*L2*cos(q2) + (L2)^2), m2*(L2)^2+eps];

v = [-m2*L1*L2*sin(q2)*(2*q1_dot*q2_dot + (q2_dot)^2);...
     m2*L1*L2*(q1_dot)^2*sin(q2)] ; 

g = [(m1 + m2)*L1*g*cos(q1)+ m2*g*L2*cos(q1+ q2); ...
     m2*g*L2*cos(q1+q2)];
 

f = [q1_dot*.25;q2_dot*.25];
 
 q_ddot = m\(tau - v -g - f); 
 
 q2_ddot = q_ddot(2);


end 
