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

N = 55;
T = .05;
fold_name = '\T18';


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
Q = zeros(4,4); Q(1,1)=10; Q(2,2) = 10;Q(3,3) =6 ;Q(4,4) =6 ; % Weighting matrices (states)
R = 1*eye(2,2); % Weighting matrices (controls)
Term_cost = zeros(4,4); Term_cost(1,1)=5; Term_cost(2,2) =5 ;Term_cost(3,3) =1 ;Term_cost(4,4) =1 ; % Weighting matrices (states)

st = X(:,1); 
g = [g; (st-P(1:4))]; % Initial Constraint via Defined Initial Condition

% con_ref = [0;0];
con_ref = Tau_Ref(P(5:6), p);

% obj_x = -.20;
% obj_y = .30;

obj_x = -.3;
obj_y = .35;

% obj_x = -.15;
% obj_y = .47;

% obj_x = .25;
% obj_y = .35;

obj_dia = .5*.125;

end_effect_dia = .5*.125;

% Compute Objective function from Stage Cost 
for k=1:N
    % Symbolically Computes Graph for the objective function
    st = X(:,k); con = U(:,k); ref = P(5:8); 
    
    L1 = .25; 
    L2 = .25;

    x = L1*cos(X(1,k)) + L2*cos(X(1,k) + X(2,k));
    y = L1*sin(X(1,k)) + L2*sin(X(1,k) + X(2,k));
    
    c = -sqrt((x - obj_x )^2 + (y - obj_y)^2) + (end_effect_dia + obj_dia +.01);
   
    obj = obj + (st - ref)'*Q*(st - ref) + (con- con_ref)'*R*(con-con_ref) + 50000*(norm(max(c,0))); % Construct Symbolic representation of objective function at time step "K"   
  
    st_next = X(:,(k+1));
    f_value = f(st,con); 
    st_next_euler = st + (T*f_value); 
    
    if k < N 
      g = [g; st_next-st_next_euler];
      
    elseif k == N 
      
      g = [g; st_next-st_next_euler; st_next - ref];
    end 
end 

% obj_x = -.35;
% obj_y = .35;
% obj_x = -0;
% obj_y = 0;

% obj_dia = .5*.125;
% 
% end_effect_dia = .5*.125;
for k = 1:N
    
    pos = forward_kinematics(p,[X(1,k);X(2,k)], [X(3,k);X(4,k)]);
    
    L1 = .25; 
    L2 = .25;

    x = L1*cos(X(1,k)) + L2*cos(X(1,k) + X(2,k));
    y = L1*sin(X(1,k)) + L2*sin(X(1,k) + X(2,k));
    
    c = -sqrt((x - obj_x )^2 + (y - obj_y)^2) ;%+ (end_effect_dia + obj_dia);
    c = 0; 
    g = [g; c];

end 


obj = obj + (st_next - ref)'*Term_cost*(st_next - ref); % Construct Symbolic representation of objective function at time step "K"   


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
args.lbg(1:4*(N+1)) = 0; 
args.ubg(1:4*(N+1)) = 0; 

args.lbg(4*(N+1):4*(N+1)+4) = 0; 
args.ubg(4*(N+1):4*(N+1)+4) = 0; 

% % Object Constraints
args.lbg(4*(N+1) + 5:2:4*(N+1) + 4 + N) = 0; % Link 1 Constraints 
args.ubg(4*(N+1) + 5:2:4*(N+1) + 4 + N) = 0; 

%Combined State & Input Constraints (due to multishooting formulation)  
args.lbx(1:4:4*(N+1),1) = -inf; % THETA1 Lower Bound
args.ubx(1:4:4*(N+1),1) = inf;  % THETA1 Upper Bound

args.lbx(2:4:4*(N+1),1) = -inf; % THETA2 Lower Bound
args.ubx(2:4:4*(N+1),1) = inf;  % THETA2 Upper Bound

% args.lbx(3:4:4*(N+1),1) = -inf; % THETA1_dot Lower Bound
% args.ubx(3:4:4*(N+1),1) = inf;  % THETA1_dot Upper Bound
% 
% args.lbx(4:4:4*(N+1),1) = -inf; % THETA2_dot Lower Bound
% args.ubx(4:4:4*(N+1),1) = inf;  % THETA2_dot Upper Bound

args.lbx(3:4:4*(N+1),1) = -3.5; % THETA1_dot Lower Bound
args.ubx(3:4:4*(N+1),1) = 3.5;  % THETA1_dot Upper Bound

args.lbx(4:4:4*(N+1),1) = -3.5; % THETA2_dot Lower Bound
args.ubx(4:4:4*(N+1),1) = 3.5;  % THETA2_dot Upper Bound

args.lbx(4*(N+1):2: 4*(N+1) + 2*N,1) = tau_1_min;   % Input Lower Bound TAU1
args.ubx(4*(N+1):2: 4*(N+1) + 2*N,1) = tau_1_max;   % Input Upper Bound TAU1

args.lbx(4*(N+1)+1:2: 4*(N+1) + 2*N,1) = tau_2_min;   % Input Lower Bound TAU2
args.ubx(4*(N+1)+1:2: 4*(N+1) + 2*N,1) = tau_2_max;   % Input Upper Bound TAU2

% Initial Goal -> Final Goal

x1_ref = .258; 
y1_ref = .3686;

% x1_ref = .2; 
% y1_ref = .45;

% x1_ref = .4; 
% y1_ref = .1;

% x1_ref = .43; 
% y1_ref = .15;
angles_0 = Inverse_Kinematics(x1_ref,y1_ref, .25, .25); 

q1_ref = angles_0(1);
q2_ref = angles_0(2);


x0_ref = [q1_ref;q2_ref; 0; 0] ; 


x2_ref = -.45;
y2_ref = .2; 

% x2_ref = 0;
% y2_ref = .25; 

% x2_ref = -.48;
% y2_ref = 0; 

% x2_ref = -.2;
% y2_ref = .43; 
angles_f = Inverse_Kinematics(x2_ref,y2_ref, .25, .25); 

q1f_ref = pi + angles_f(1); 
q2f_ref = angles_f(2);

xf_ref = [q1f_ref;q2f_ref; 0; 0] ;


% SIMULATION LOOP 
t0 = 0; 

x_list(:,1) = x0_ref; 
% x_list = [x_list, x0_ref]; 

t(1) = t0; 
u0 = zeros(2,N);
X0 = repmat(x0_ref,1, N+1)';
sim_time = 20; 

% Start MPC 
mpc_iter = 0; 
xxl = [ ]; 
u_cl = [ ];

pred_state_horizon = [] ; 


main_loop = tic;

while(norm((x0_ref - xf_ref),2) > 1e-3 && mpc_iter < sim_time/T)
    
    itter_loop = tic;
    
    args.p = [x0_ref;xf_ref];
    args.x0 =  [reshape(X0', 4*(N+1),1); reshape(u0,2*N,1)];
    
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p); 
    
    u = reshape(full(sol.x(4*(N+1)+1:end))',2,N)'; 
    xxl(:, 1:4, mpc_iter +1) = reshape(full(sol.x(1:4*(N+1)))',4,N+1)';
    u_cl = [u_cl; u(1,:)];

    t(mpc_iter +1) = t0; 
    [t0, x0_ref, u0] = shift(T, t0, x0_ref, u, f);
    x_list(:,mpc_iter+2) = x0_ref; 
%     x_list = [x_list, x0_ref]; 

    mpc_iter = mpc_iter +1 ;  
    
    itter_loop_time(mpc_iter) = toc(itter_loop);
    
end 
disp('Single Loop Stats')
mean(itter_loop_time)
disp('Sim Time')
length(t)*T


u_cl = u_cl';
main_loop_time = toc(main_loop)

%%
figure(3)
video_path = 'C:\Users\Indy-Windows\Documents\ECH-267-Adv.-Proc.-Control\Project\Report\images\OBS_Avoid_TEST';
v = VideoWriter(strcat(video_path, fold_name,'\robot_MPC.avi')); 
v.FrameRate = 20;
open(v); 
for i=1:1:length(x_list)-1
    
    
    theta = x_list(1:2,i); 
    theta_dot = x_list(3:4,i); 
    pos = forward_kinematics(p,theta, theta_dot);
    
    x_vec = pos.x;
    y_vec = pos.y;
%     plot_robot_with_object(x_vec, y_vec, obj_x, obj_y, x2_ref, y2_ref,x1_ref, y1_ref)
%     plot_robot(x_vec, y_vec, x1_ref, y1_ref, x2_ref, y2_ref);
    
    var = size(xxl);
    hold on
    for k = 1:1:var(1)

        q1 = xxl(k, 1, i);
        q2 = xxl(k, 2, i);
        
        
        pos = forward_kinematics(p, [q1,q2], [0,0]); 
        
        x = pos.x; 
        y = pos.y;
        
%         viscircles([x(2), y(2)], .0625,'LineStyle','--', 'Color', 'k');
          plot(x(2), y(2), '*k')
    end 
        plot_robot_with_object(x_vec, y_vec, obj_x, obj_y, x2_ref, y2_ref,x1_ref, y1_ref)

    hold off
      frame=getframe(gcf); 
        writeVideo(v, frame);
    pause(.25)
    fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background

    
  
end 

close(v)

%%

main_loop_time = toc(main_loop)


u_cl = u_cl;
 
figure(1)
subplot(211) 
plot(x_list(1,:));
ylabel('Theta_1 [rads]'); 
title(' MPC Obstacle Avoidance: Positions')


subplot(212) 
plot(x_list(2,:));
ylabel('Theta_2 [rads]');

figure(2)
subplot(211)
plot(x_list(3,:)); 
ylabel('Ang. Vel Theta_1 [rad/s]'); 
xlabel('Time [seconds]');
title('MPC Obstacle Avoidance: Velocities')


subplot(212)
plot(x_list(4,:)); 
ylabel('Ang. Vel Theta_2 [rad/s]'); 
xlabel('Time [seconds]');


figure(3)
grid on
subplot(211)
stairs(t,u_cl(1,:),'k','linewidth',1.5)
ylabel('Torque_1 (N*m)')
title('MPC Obstacle Avoidance: Torque Inputs')

subplot(212)
stairs(t,u_cl(2,:),'k','linewidth',1.5); axis([0 t(end) -0.35 0.75])
xlabel('time (seconds)')
ylabel('Torque_2 (N*m)')
disp("DONE");
disp(rad2deg(x_list(1,end)))


%% Graphs


fn_torque = '\torque.png';
fn_vel = '\vel.png';
fn_pos = '\pos.png';


path = 'C:\Users\Indy-Windows\Documents\ECH-267-Adv.-Proc.-Control\Project\Report\images\OBS_Avoid_TEST';

f_path_torque = strcat(path, fold_name,fn_torque);
f_path_vel = strcat(path, fold_name,fn_vel);
f_path_pos = strcat(path, fold_name,fn_pos);


% Save Torques 
saveas(figure(3),f_path_torque)

% Save Velocities 
 saveas(figure(2),f_path_vel)

% Save Positions
saveas(figure(1),f_path_pos)






%% End Effector Trajectory 

x1_traj = []; 
y1_traj = [];

x2_traj = []; 
y2_traj = [];


for i=1:1:length(x_list)
    
    
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

function plot_robot(x_vec, y_vec, x0_ref, y0_ref, xf_ref, yf_ref)

x1 = x_vec(1); 
x2 = x_vec(2); 
y1 = y_vec(1); 
y2 = y_vec(2); 



hold on
fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background

plot([ 0, x1], [0, y1], "b");
plot([ x1, x2], [y1, y2], "r");

plot(x0_ref, y0_ref, '-*r')
plot(xf_ref, yf_ref, '-*g')

hold off

title('Planar 2D Robot in Workspace')
xlim([-.6,.6])
ylim([-.6,.6])

xlabel("X-Axis")
ylabel("Y-Axis") 



end

function plot_robot_with_object(x_vec, y_vec, obj_x, obj_y, goal1, goal2, start1, start2)

x1 = x_vec(1); 
x2 = x_vec(2); 
y1 = y_vec(1); 
y2 = y_vec(2); 

x_goal = goal1;
y_goal = goal2;

x_start = start1; 
y_start = start2; 


hold on
% fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background


viscircles([obj_x,obj_y], .0625,'LineStyle','--');


plot([ 0, x1], [0, y1], "b");
plot([ x1, x2], [y1, y2], "r");

plot(x2, y2, '-*b')

viscircles([x2, y2], .0625,'LineStyle','--', 'Color', 'b');


plot(obj_x,obj_y, '*r') % Obstacle Position 
plot(x_start,y_start, '*k') % Start Position
plot(x_goal,y_goal, '*g') % End Position 

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


x = [L1*cos(q1) ; L1*cos(q1) + L2*cos(q1 + q2)];
y = [L1*sin(q1) ; L1*sin(q1) + L2*sin(q1 + q2)];

pos = struct ;

pos.x = x; 
pos.y = y; 

end 

function angles = Inverse_Kinematics(x,y, L1, L2)

q2 = acos(((x)^2 +(y)^2-(L1)^2-(L2)^2)/(2*L1*L2));
q1 = atan(y/x) - atan((L2*sin(q2))/(L1+L2*cos(q2)));

angles = [q1,q2];


end 

function pos_cg = CG_forward_kinematics(p,theta, theta_dot)
 
L1 = p.L1; 
L2 = p.L2; 
g = p.g; 

q1 = theta(1);
q2 = theta(2);

q1_dot = theta_dot(1); 
q2_dot = theta_dot(2); 


x = [.5*L1*cos(q1) ;L1*cos(q1) + .5*L2*cos(q1 + q2)];
y = [.5*L1*sin(q1) ;L1*sin(q1) + .5*L2*sin(q1 + q2)];

pos_cg = struct ;

pos_cg.x = x; 
pos_cg.y = y; 

end 

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

function out = Tau_Ref(theta, p)

q1 = theta(1);
q2 = theta(2);


m1 = p.m1; 
m2 = p.m2; 
L1 = p.L1; 
L2 = p.L2; 
g = p.g; 

out = [(m1 + m2)*L1*g*cos(q1)+ m2*g*L2*cos(q1+ q2); ...
     m2*g*L2*cos(q1+q2)];


end 



function [x0_2, y0_2] = T0_2(q1, q2)





end 


