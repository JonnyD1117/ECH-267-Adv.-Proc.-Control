%% Problem6.2 

[x1, x2] = meshgrid(-.25:.1:1.5, -.25:.1:1.5); 

x1_dot = -x1.^3 + x2; 
x2_dot = x1.^6 + -x2.^3; 


  figure(1)

 hold on
 quiver(x1, x2, x1_dot, x2_dot)
 plot( 0, 0, 'r*')
  plot( 1, 1, 'r*')
  
  
    x1_0_list = [.25, .5 .75]; 
  x2_0_list = [.25 .5 .75];
  
  dt = .1; 
  
  sim_time = 10;
  
  num_step = 10/dt;
  
  num_traj = length(x1_0_list); 
  
  for j = 1:1:num_traj
      
      x1_0 =x1_0_list(j) ; 
      x2_0 = x2_0_list(j); 
        for i = 1:1:num_step
      
      x1 = x1_0; 
      x2 = x2_0; 
      
      x1_dot = -x1.^3 + x2; 
    x2_dot = x1.^6 + -x2.^3;  
      
      x1_traj(i) = x1_0 + dt*x1_dot; 
      x2_traj(i) = x2_0 + dt*x2_dot;
      
      x1_0 =x1_traj(i); 
      x2_0 = x2_traj(i);
        end 
      
        plot(x1_traj, x2_traj, 'k')

  end 
  
  
  title("Problem 6.2")
  
 




%% 6.2 continued

x1 = [0:.1:1];

c1_x2 = x1.^3;
c2_x2 = x1.^2; 



figure()
hold on 

plot(x1, c1_x2)
plot(x1, c2_x2)

x = [x1, fliplr(x1)];
inbetween = [c1_x2, fliplr(c2_x2)];



% fill(x, inbetween, 'g');
fill(x, inbetween, '');
xlabel("X1-Axis")
ylabel("X2-Axis")
title("Problem 6.2:")

legend('$x_2 = x_1^3$','$x_2 = x_1^2$','$\Gamma$','Interpreter','latex', 'Location', 'northwest')


%%

x = 1 : 300;
curve1 = log(x);
curve2 = 2*log(x);
plot(x, curve1, 'r', 'LineWidth', 2);
hold on;
plot(x, curve2, 'b', 'LineWidth', 2);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'g');



%% Problem 10 
[x1, x2] = meshgrid(-1:0.2:1, -1:0.2:1);


x1_dot = x2; 
x2_dot = -x1 + (1/3)*x1^3 -x2; 

  figure(1)

 hold on
 quiver(x1, x2, x1_dot, x2_dot)
 plot( 0, 0, 'r*')
  title(' Problem 10')


  
  x1_0_list = [-1,-.5,.5,1, 1, 1, 1, 1]; 
  x2_0_list = [-1, -1, -1, -1, -1,-.5,.5,1];
  
  dt = .1; 
  
  sim_time = 10;
  
  num_step = 10/dt;
  
  num_traj = length(x1_0_list); 
  
  for j = 1:1:num_traj
      
      x1_0 =x1_0_list(j) ; 
      x2_0 = x2_0_list(j); 
        for i = 1:1:num_step
      
      x1 = x1_0; 
      x2 = x2_0; 
      
      x1_dot = x2; 
      x2_dot = -x1 + (1/3)*x1^3 -x2; 
      
      x1_traj(i) = x1_0 + dt*x1_dot; 
      x2_traj(i) = x2_0 + dt*x2_dot;
      
      x1_0 =x1_traj(i); 
      x2_0 = x2_traj(i);
        end 
      
        plot(x1_traj, x2_traj, 'k')

  end 
    
 V = (3/4).*x1.^2 - (1/12).*x1.^4 + .5.*x1.*x2 + .5.*x2.^2; 
  
  contour(V) 
  %% 10 Contour Plot 
  
  
  [x1, x2] = meshgrid(-2.5:.2:2.5,-2.5:.2:2.5);
  
  
  
  
  V = (3/4).*x1.^2 - (1/12).*x1.^4 + .5.*x1.*x2 + .5.*x2.^2; 

  
  
  figure(2)
  
  hold on
  contour(V, 25)
  
  
  %% 17 
  
  
  syms a b
  
  x1 = sqrt(-a/b);
  
  A = [-3*x1^2 1 ; -a -b]
  
  latex(eig(A))
  
  eig(A)
  
  
  