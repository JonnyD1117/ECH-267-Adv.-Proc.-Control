import casadi.*


function [outputArg1,outputArg2] = MPC_Controller(state,state_reference, input, input_reference, set_up_flag)


if set_up_flag == 1
    

p = struct ;

p.m1 = .5;
p.m2 = .25; 
p.L1 = .25;
p.L2 = .25; 
p.g = 9.81;

N = 15;
T = .1;

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
q_ddot = RobotModel_q_ddot(p, x1, x2, x3, x4,tau1, tau2);
rhs = [x3; x4 ; q_ddot(1);q_ddot(2)]; 

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

% con_ref = [0;0];
con_ref = Tau_Ref(P(5:6), p);

phi = deg2rad(15); 
alpha = deg2rad(10); 


% Compute Objective function from Stage Cost 
for k=1:N
    % Symbolically Computes Graph for the objective function
    st = X(:,k); con = U(:,k); ref = P(5:8); 

    obj = obj + (st - ref)'*Q*(st - ref) + (con- con_ref)'*R*(con-con_ref); % Construct Symbolic representation of objective function at time step "K"   
  
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



%Combined State & Input Constraints (due to multishooting formulation)  
% args.lbx(1:4:4*(N+1),1) = -inf; % THETA1 Lower Bound
% args.ubx(1:4:4*(N+1),1) = inf;  % THETA1 Upper Bound
% 
% args.lbx(2:4:4*(N+1),1) = -inf; % THETA2 Lower Bound
% args.ubx(2:4:4*(N+1),1) = inf;  % THETA2 Upper Bound

args.lbx(1:4:4*(N+1),1) = -(.5*pi-phi); % THETA1 Lower Bound
args.ubx(1:4:4*(N+1),1) = (3/2)*pi - phi;  % THETA1 Upper Bound

args.ubx(2:4:4*(N+1),1) = pi - alpha; % THETA2 Lower Bound
args.lbx(2:4:4*(N+1),1) = alpha - pi;  % THETA2 Upper Bound

args.lbx(3:4:4*(N+1),1) = -inf; % THETA1_dot Lower Bound
args.ubx(3:4:4*(N+1),1) = inf;  % THETA1_dot Upper Bound

args.lbx(4:4:4*(N+1),1) = -inf; % THETA2_dot Lower Bound
args.ubx(4:4:4*(N+1),1) = inf;  % THETA2_dot Upper Bound


args.lbx(4*(N+1):2: 4*(N+1) + 2*N,1) = tau_1_min;   % Input Lower Bound TAU1
args.ubx(4*(N+1):2: 4*(N+1) + 2*N,1) = tau_1_max;   % Input Upper Bound TAU1

args.lbx(4*(N+1)+1:2: 4*(N+1) + 2*N,1) = tau_2_min;   % Input Lower Bound TAU2
args.ubx(4*(N+1)+1:2: 4*(N+1) + 2*N,1) = tau_2_max;   % Input Upper Bound TAU2

    
    
    
else 
    
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
    
    u_cl = [u_cl; u(1,:)];
    
    t(mpc_iter +1) = t0; 
    [t0, x0, u0] = shift(T, t0, x0, u, f);
    x_list(:,mpc_iter+2) = x0; 
    mpc_iter = mpc_iter +1 ;   
    
    
    
end 






end

