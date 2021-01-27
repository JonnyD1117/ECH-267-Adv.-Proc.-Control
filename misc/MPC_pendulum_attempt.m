%% Simple Pendulum MPC Attempt (Single Shooting) 

import casadi.*

clear all 
close all 
clc 


% System Parameters & Constants
T = .2; % Sampling Period 
N = 10; % Time Horizon 
g = -9.81; % Acceleration due to Gravity
b = .25; % Damping Const.
L = 1; % Pendulum Length 

% Control Input Saturation Limits
tau_max = 20; % N*m
tau_min = -20; %N*m

% Define Symbolic Variables
x1 = SX.sym('x1'); x2 = SX.sym('x2'); tau = SX.sym('tau'); 

states = [x1 ; x2]; n_states = length(states);
ctrl_input = [tau]; n_ctrl = length(ctrl_input); 

% Define Symbolic State Space RHS
rhs = [x2; (-(g/L)*sin(x1) - b*x2 +tau)]; 

% Define Symbolic CasADi function for state dynamics 
f = Function('f', {states, ctrl_input}, {rhs}, {'state','input'}, {'dynamics'}); % This function enables us to numerically evaluate the symbolic graph given numerical inputs

% Define Symbolic Matrix for Inputs & Params for each step of Time Horizon 
U = SX.sym('U', n_ctrl, N); % Inputs at each step K in horizon N 
P = SX.sym('P', 2*n_states);% Params at each step K in horizon N, P = (x1_0, x2_0, x1_ref, x2_ref)
X = SX.sym('X', n_states, (N+1)); % States at each step K in horizon N

% Initialize the State initial state vector using Parameter vector 
X(:,1) = P(1:2); 

% Using forward recursion, populate X matrix symbollically (Computational Graph)  
for k = 1:N
    st = X(:,k); con = U(:,k); % Define the current state and input    
    f_value = f(st, con); % Use 'st' & 'con' to compute the state derivatives
    st_new = st + (T*f_value); % Update Dynamics using Euler Integration
    X(:,k+1) = st_new;    % Update state at next time step K+1
end 

% Using this symbolic representation of the State Update process, create a function to implement this behavior 
ff = Function('ff',{U,P},{X}, {'Input','Reference '},{'New State'}); % This enables evaluation of symbolic graph give numerical inputs

% Define Objective & Constraints Vector 
obj = 0; % Initialize the Objective Func. to zero
g = [] ; % Initialize the Constraints as Empty 

% Define Loss function Weights on the States and Control Inputs
Q = zeros(2,2); Q(1,1)=1; Q(2,2) =2 ; % Weighting matrices (states)
R = .5; % Wieghting matrices (controls)

% Compute Objective function from Stage Cost 
for k=1:N
    % Symbolically Computes Graph for the objective function
    st = X(:,k); con = U(:,k);
    ref = P(3:4); 
    obj = obj + (st - ref)'*Q*(st - ref) + con*R*con;   
end 


% % Compute the Constraints Vector 
% for k=1:N+1
%     g = [g; X(1,k)]; % States X 
%     g = [g; X(2,k)];
% end 


% Define Nonlinear Programming Structure 
OPT_variables = U;  % Via Single Shooting only the Control Input set is used as the decision variable
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
opts = struct; 

% Define Solver Options 
opts.ipopt.max_iter =100; 
opts.ipopt.print_level=0; % 0 - 3
opts.print_time  = 0 ; 
opts.ipopt.acceptable_tol = 1e-8; 
opts.ipopt.acceptable_obj_change_tol = 1e-6; 

% Enstantiate nlpsol class using ipopt 
solver = nlpsol('solver','ipopt', nlp_prob, opts); 

% Define Problem Arguments 
args = struct;

% Inequality Constraints
args.lbg = 0; 
args.ubg = 0; 

%Input Constraints 
args.lbx(1:1:N,1) = tau_min;
args.ubx(1:1:N,1) = tau_max;


% SIMULATION LOOP 
t0 = 0; 
x0 = [(5*pi/4); 0] ;
x_ref = [pi;0] ;

x_list(:,1) = x0; 
t(1) = t0; 
u0 = zeros(N,1);
sim_time = 20; 

% Start MPC 
mpc_iter = 0; 
xxl = [ ]; 
u_cl = [ ];


% Main Loop 
main_loop = tic;

while(norm((x0 - x_ref),2) > 1e-2 && mpc_iter < sim_time/T)
    
    args.p = [x0;x_ref];
    args.x0 =  u0;
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p); 
    
    u = reshape(full(sol.x)',1,N)'; 
    ff_value = ff(u, args.p); 
    xxl(:, 1:2, mpc_iter +1) = full(ff_value)';
    u_cl = [u_cl; u(1)]; 
    
    t(mpc_iter +1) = t0; 
    [t0, x0, u0] = shift(T, t0, x0, u, f);
    x_list(:,mpc_iter+2) = x0; 
    mpc_iter = mpc_iter +1 ; 
    
    
end 
main_loop_time = toc(main_loop)


figure
grid on

stairs(t,u_cl(:),'k','linewidth',1.5); axis([0 t(end) -0.35 0.75])
xlabel('time (seconds)')
ylabel('Torque (N*m)')
disp("DONE");

Animate_pendulum(x_list, T, sim_time)


function [t0, x0, u0] = shift(T, t0, x0, u, f)

st = x0; 
con = u(1,:)'; 
f_value = f(st,con); 
st = st+ (T*f_value); 
x0 = full(st); 
t0 = t0 + T;
u0 = [u(2:size(u,1),:); u(size(u,1),:) ];  


end

function Animate_pendulum(x_traj, dt, sim_time)

v = VideoWriter('pendulum.avi'); 
v.FrameRate = 30;
open(v); 

x1 = x_traj(1,:); 

sim_len = (sim_time / dt);

disp(sim_len)

    for k=2:length(x1)
        %Plot for Video 
        hold  on
        fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background
        plot([0 sin(x1(k-1))],[0 -cos(x1(k-1))], 'b', 'LineWidth', 3); % Plots Rod
    %     plot(,'Marker','o', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); % Plots Bob 

        xlim([-1-0.2*1 1+0.2*1]); 
        ylim([-1-0.2*1 1+0.2*1]);
        title('Simple Pendulum');
        frame=getframe(gcf); 
        writeVideo(v, frame);

        disp(k)
    end 


end 
