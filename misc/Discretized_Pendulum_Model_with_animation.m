

%% Pendulum Equations
g = -9.81; 
L = 1; 
b = 1.2;
dt = .1;
t_tot = 100; 
sim_len = t_tot/dt;

t = [0:dt:(t_tot-.1)];

% % Define Initial Conditions 
% x1_0 = pi/4 ; 
% x2_0 = 0; 




x1 = zeros(1,sim_len); 
x2 = zeros(1,sim_len); 



%% Euler Discretization for Pendulum 

% Define Initial Conditions 
x1_0 = pi/4 ; 
x2_0 = 0;

% Store Initial Condition in State Vectors
x1(1) = x1_0; 
x2(1) = x2_0; 

for k=2:sim_len
    
    % State Dynamics
    x1_dot = x2(k-1);
    x2_dot = (g/L)*sin(x1(k-1)) - b*x1_dot;  
        
    % Propagate State Update through Discretized Model 
    x1(k) = x1(k-1) + x1_dot*dt; 
    x2(k) = x2(k-1) + x2_dot*dt;
end 


figure(1);
plot(t, x1);
xlabel('Time [sec]');
ylabel('Angle-Theta [rads]');
title('Theta Vs Time');

figure(2);
plot(t, x2);
xlabel('Time [sec]');
ylabel('Angular-Velocity [rads/sec]');
title('Theta_dot Vs Time');

%% Animation 


% Define Initial Conditions 
x1_0 = pi/4 ; 
x2_0 = 0;

% Store Initial Condition in State Vectors
x1(1) = x1_0; 
x2(1) = x2_0; 


v = VideoWriter('pendulum.avi'); 
v.FrameRate = 30;
open(v); 

for k=2:sim_len
%     
%     % State Dynamics
%     x1_dot = x2(k-1);
%     x2_dot = (g/L)*sin(x1(k-1));  
%         
%     % Propagate State Update through Discretized Model 
%     x1(k) = x1(k-1) + x1_dot*dt; 
%     x2(k) = x2(k-1) + x2_dot*dt;
    
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
    
end 


close(v); 




%% Runge-Kutta 4 Discretization for Pendulum 

% % Define Initial Conditions 
% x1_0 = pi/4 ; 
% x2_0 = 0;
% 
% % Store Initial Condition in State Vectors
% x1(1) = x1_0; 
% x2(1) = x2_0; 
% 
% for k=2:sim_len
%     % Define the TimeStep
%     h = dt; 
%     
%     % State Dynamics
%     
%         % At time 'K'
%         x1_dot_k = ;
%         x2_dot_k = ;  
%         
%         % At time 'K+.5h', 
%         x1_dot_k = x2(k-1 );
%         x2_dot_k = (g/L)*sin(x1(k-1)); 
%     
%     
%     
%     %RK4 Parameter Weights
%     k1_1 = x2(k-1); 
%     k2_1 = (g/L)*sin(x1(k-1)); 
%     k3_1 = ; 
%     k4_1 = ; 
%     
%     k1_2 = ; 
%     k2_2 = ; 
%     k3_2 = ; 
%     k4_2 = ; 
%         
%     % Propagate State Update through Discretized Model 
%     x1(k) = x1(k-1) + h*(1/6)*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1); 
%     x2(k) = x2(k-1) + h*(1/6)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2); 
% end 
% 
% 
% figure(1);
% plot(t, x1);
% xlabel('Time [sec]');
% ylabel('Angle-Theta [rads]');
% title('Theta Vs Time');
% 
% figure(2);
% plot(t, x2);
% xlabel('Time [sec]');
% ylabel('Angular-Velocity [rads/sec]');
% title('Theta_dot Vs Time');


