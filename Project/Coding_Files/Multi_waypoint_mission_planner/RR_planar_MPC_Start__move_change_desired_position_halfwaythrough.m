%%

clear all 
close all
clc


MPC_config_file

% Initial Goal -> First Way Point

x1_ref = .258; 
y1_ref = .3686;

x2_ref = -.48;
y2_ref = 0; 



angles_0 = Inverse_Kinematics(x1_ref,y1_ref, .25, .25); 

q1_ref = angles_0(1);
q2_ref = angles_0(2);


x0_ref = [q1_ref;q2_ref; 0; 0] ; 




angles_f = Inverse_Kinematics(x2_ref,y2_ref, .25, .25); 

q1f_ref = pi + angles_f(1); 
q2f_ref = angles_f(2);

xf_ref = [q1f_ref;q2f_ref; 0; 0] ;


% SIMULATION LOOP 
t0 = 0; 

x_list(:,1) = x0_ref; 
% x_list = [x_list, x0_ref]; 

t(1) = t0; 

sim_time = 20; 

% Start MPC 
mpc_iter = 0; 
xxl = [ ]; 
u_cl = [ ];

pred_state_horizon = [] ; 


main_loop = tic;

while(norm((x0_ref - xf_ref),2) > 1e-4 && mpc_iter < sim_time/T)
    
    
    if mpc_iter == 10
        xf_ref = [deg2rad(200);deg2rad(96); 0; 0] ;

    end 
    
    


    mpc_iter = mpc_iter +1 ;     
    
end 

u_cl = u_cl';
main_loop_time = toc(main_loop)

figure(3)
for i=1:1:length(x_list)
    
    
    theta = x_list(1:2,i); 
    theta_dot = x_list(3:4,i); 
    pos = forward_kinematics(p,theta, theta_dot);
    
    x_vec = pos.x;
    y_vec = pos.y;
    plot_robot_with_object(x_vec, y_vec, obj_x, obj_y, x2_ref, y2_ref,x1_ref, y1_ref)
%     plot_robot(x_vec, y_vec, x1_ref, y1_ref, x2_ref, y2_ref);
    
    var = size(xxl);
    hold on
    for k = 1:3:var(1)

        q1 = xxl(k, 1, i);
        q2 = xxl(k, 2, i);
        
        
        pos = forward_kinematics(p, [q1,q2], [0,0]); 
        
        x = pos.x; 
        y = pos.y;
        
        viscircles([x(2), y(2)], .0625,'LineStyle','--', 'Color', 'k');
    end 
    
    hold off
    pause(.1)
end 

disp(rad2deg(x_list(1,end)))







%%



function [init_input_horizon, init_state_horizon]= MPC_init(init_state, pred_horizon)

    % Argument Translation 
    N = pred_horizon;

    % Mapping Initial Inputs and States onto the Initial Input & State
    % Horizons
    init_input_horizon = zeros(2,N);
    init_state_horizon = repmat(init_state,1, N+1)';

end 

function [opt_input, pred_state_horizon ]= MPC_step(x_init, x_goal,  solver, solver_args, time_step, pred_horizion, mpc_iter )


x_pred_horizon, u_pred_horizon,

    % Argument Translation 
    args = solver_args; 
    T = time_step;
    N = pred_horizion;
    
    
    if mpc_iter == 1
        [X0, u0] = MPC_init(x_init, N);   
    else
        [t0, x0_ref, u0] = shift(T, t0, x_goal, u, f);
    end
     
    % Initialize Arguments
    args.p = [x_init;x_goal];
    args.x0 =  [reshape(X0', 4*(N+1),1); reshape(u0,2*N,1)];
    
    
    % Solve MPC at Current Time Step
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p); 
    
    % Extract Optimal Input and State Prediction Horizion from Solution
    u = reshape(full(sol.x(4*(N+1)+1:end))',2,N)'; 
    xxl(:, 1:4, mpc_iter +1) = reshape(full(sol.x(1:4*(N+1)))',4,N+1)';
    u_cl = [u_cl; u(1,:)];
    t(mpc_iter +1) = t0; 
    
    x_list(:,mpc_iter+2) = x0_ref; 
%     x_list = [x_list, x0_ref]; 





end 







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
fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background


% viscircles([obj_x,obj_y], .0625,'LineStyle','--');


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

% function q1_ddot = RobotModel_q1_ddot(p, q1, q2, q1_dot, q2_dot,tau1, tau2)
% 
% tau = [tau1;tau2] ;
% 
% m1 = p.m1; 
% m2 = p.m2; 
% L1 = p.L1; 
% L2 = p.L2; 
% g = p.g; 
% 
% eps = .1 ;
% 
% m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(q2) + (L2)^2) + eps, m2*(L1*L2*cos(q2)+ (L2)^2); ...
%      m2*(L1*L2*cos(q2) + (L2)^2), m2*(L2)^2+eps];
% 
% v = [-m2*L1*L2*sin(q2)*(2*q1_dot*q2_dot + (q2_dot)^2);...
%      m2*L1*L2*(q1_dot)^2*sin(q2)] ; 
% 
% g = [(m1 + m2)*L1*g*cos(q1)+ m2*g*L2*cos(q1+ q2); ...
%      m2*g*L2*cos(q1+q2)];
%  
% 
% f = [q1_dot*.25;q2_dot*.25];
%  
%  q_ddot = m\(tau - v -g - f); 
%  
%  q1_ddot = q_ddot(1);
% 
% 
% end 
% 
% function q2_ddot = RobotModel_q2_ddot(p, q1, q2, q1_dot, q2_dot,tau1, tau2)
% 
% tau = [tau1;tau2] ;
% 
% m1 = p.m1; 
% m2 = p.m2; 
% L1 = p.L1; 
% L2 = p.L2; 
% g = p.g; 
% 
% eps = .1 ;
% 
% m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(q2) + (L2)^2) + eps, m2*(L1*L2*cos(q2)+ (L2)^2); ...
%      m2*(L1*L2*cos(q2) + (L2)^2), m2*(L2)^2+eps];
% 
% v = [-m2*L1*L2*sin(q2)*(2*q1_dot*q2_dot + (q2_dot)^2);...
%      m2*L1*L2*(q1_dot)^2*sin(q2)] ; 
% 
% g = [(m1 + m2)*L1*g*cos(q1)+ m2*g*L2*cos(q1+ q2); ...
%      m2*g*L2*cos(q1+q2)];
%  
% 
% f = [q1_dot*.25;q2_dot*.25];
%  
%  q_ddot = m\(tau - v -g - f); 
%  
%  q2_ddot = q_ddot(2);
% 
% 
% end 

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



