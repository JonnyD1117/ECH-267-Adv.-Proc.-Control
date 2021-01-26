import casadi.*

% Rosenbrock Problem

% Create Symbolic Variables
x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z');

% Populate NLP Structure (%Optimizer Args, Cost funct, Constraints)
nlp = struct('x',[x;y;z], 'f',x^2+100*z^2, 'g',z+(1-x)^2-y);

% Define the NLP Solver ( Func handle, optimizer algo, problem structure )
S = nlpsol('S', 'ipopt', nlp);

% Display the Solver Seetings
disp(S)

% Pass Initial Conditions into the Solver, along with parameters such as
% bounds. 
r = S('x0',[2.5,3.0,0.75],...
      'lbg',0,'ubg',0);

% Extract Argmax  
x_opt = r.x;

% Displace Argmax
disp(x_opt)

% NOTE: r.x is NOT the x component of the multivariate vector arguments,
% 'x' is the function handle for the optimization argument defined in LINE
% 9 for the NLP struct. In the same vein, 'x0' is NOT the initial value of
% the x variable, it is the initial value(S) for the optimizer arguments 