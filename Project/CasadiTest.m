% addpath('C:\Users\Indy-Windows\Desktop\Winter 2021\Casadi\Matlab/casadi-matlabR2016a-v3.5.5')
import casadi.*

clear all
close all
clc


T = 0.2;  % Sampling Time 
N = 3; % Prediction horizon 
rob_diam = 0.3; 

v_max = .6; v_min=-v_max; 
omega_max = pi/4; omega_min = -omega_max; 

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states); 

v=SX.sym('v'); omega = SX.sym('omega'); 
controls = [v;omega]; n_controls = length(controls); 
rhs = [v*cos(theta);v*sin(theta); omega]; 


f = Function('f', {states, controls}, {rhs}); 
U = SX.sym('U', n_controls, N);
P = SX.sym('P',n_states + n_states);

X = SX.sym('X', n_states, (N+1)); 


% Compute Solution Symbolically
X(:,1) = P(1:3); 

for k=1:N
    
    st = X(:,k); con = U(:,k); 
    f_value = f(st,con); 
    st_next = st + (T*f_value);
    X(:,k+1) = st_next;
    
end 

% This function to get the optimal trajectory knowning the optimal solution
ff = Function('ff', {U,P},{X});


obj = 0;  % objective function 
g = []; % constraints vector 


Q = zeros(3,3); Q(1,1)=1; Q(2,2) =5 ; Q(3,3) = 0.1; % Wieghting matrices (states)
R = zeros(2,2); Q(1,1)=0.5; Q(2,2) = 0.05; % Wieghting matrices (controls)

% compute Objective 

for k=1:N
    
    st = X(:,k); con = U(:,k); 
    obj = obj + (st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con;
    
end 

% Compute Constraints 

for k=1:N+1
    
    g = [g; X(1,k)]; % States X 
    g = [g; X(2,k)];
end 


OPT_variables = reshape(U,2*N,1); % Reshaped matrix to a single vector 
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P); 


opts = struct;
opts.ipopt.max_iter = 100; 
opts.ipopt.print_level = 0; 
opts.print_time = 0; 
opts.ipopt.acceptable_tol = 1e-8; 
opts.ipopt.acceptable_obj_change_tol = 1e-6; 

solver = nlpsol('solver','ipopt', nlp_prob, opts); 

args = struct; 
args.lbg = -2;
args.ubg = 2; 

args.lbx(1:2:2*N-1,1) = v_min; args.lbx(2:2:2*N,1)= omega_min; 
args.ubx(1:2:2*N-1,1) = v_max; args.ubx(2:2:2*N,1) = omega_max; 




t0 =0 ; 
x0 = [0; 0 ; 0.0]; 
xs= [1.5; 1.5; 0.0];
xx(:,1) = x0; 
t(1) = t0; 
u0 = zeros(N,2); 
sim_tim = 20; 
mpciter = 0; 
xx1 = [] ;
u_cl = [] ;

main_loop = tic;
while(norm((x0 - xs),2) > 1e-2 && mpciter < sim_tim /T)
    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0 = reshape(u0',2*N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)',2,N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:,1:3,mpciter+1)= full(ff_value)';
    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0] = shift(T, t0, x0, u,f); % get the initialization of the next optimization step
    
    xx(:,mpciter+2) = x0;  
    mpciter
    mpciter = mpciter + 1;
end;
main_loop_time = toc(main_loop)

ss_error = norm((x0-xs),2)
Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam) % a drawing function



function [t0,x0,u0] = shift(T,t0,x0,u,f)

st = x0; 
con = u(1,:)';
f_value = f(st,con);
st = st + (T*f_value);
x0 = full(st); 

t0 = t0 + T; 
u0 = [u(2:size(u,1),:) ; u(size(u,1),:)];


end 


function Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam)


set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 1.5;
fontsize_labels = 14;

%--------------------------------------------------------------------------
%-----------------------Simulate robots -----------------------------------
%--------------------------------------------------------------------------
x_r_1 = [];
y_r_1 = [];



r = rob_diam/2;  % obstacle radius
ang=0:0.005:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);

figure(500)
% Animate the robot motion
%figure;%('Position',[200 200 1280 720]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55 1]);

for k = 1:size(xx,2)
    h_t = 0.14; w_t=0.09; % triangle parameters
    
    x1 = xs(1); y1 = xs(2); th1 = xs(3);
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];
    fill(x1_tri, y1_tri, 'g'); % plot reference state
    hold on;
    x1 = xx(1,k,1); y1 = xx(2,k,1); th1 = xx(3,k,1);
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];

    plot(x_r_1,y_r_1,'-r','linewidth',line_width);hold on % plot exhibited trajectory
    if k < size(xx,2) % plot prediction
        plot(xx1(1:N,1,k),xx1(1:N,2,k),'r--*')
    end
    
    fill(x1_tri, y1_tri, 'r'); % plot robot position
    plot(x1+xp,y1+yp,'--r'); % plot robot circle
    
   
    hold off
    %figure(500)
    ylabel('$y$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    axis([-0.2 1.8 -0.2 1.8])
    pause(0.1)
    box on;
    grid on
    %aviobj = addframe(aviobj,gcf);
    drawnow
    % for video generation
    F(k) = getframe(gcf); % to get the current frame
end
close(gcf)
%viobj = close(aviobj)
%video = VideoWriter('exp.avi','Uncompressed AVI');

% video = VideoWriter('exp.avi','Motion JPEG AVI');
% video.FrameRate = 5;  % (frames per second) this number depends on the sampling time and the number of frames you have
% open(video)
% writeVideo(video,F)
% close (video)

figure
subplot(211)
stairs(t,u_cl(:,1),'k','linewidth',1.5); axis([0 t(end) -0.35 0.75])
ylabel('v (rad/s)')
grid on
subplot(212)
stairs(t,u_cl(:,2),'r','linewidth',1.5); axis([0 t(end) -0.85 0.85])
xlabel('time (seconds)')
ylabel('\omega (rad/s)')
grid on

end

