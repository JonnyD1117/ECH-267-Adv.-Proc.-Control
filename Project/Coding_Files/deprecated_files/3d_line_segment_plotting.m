
t = [0:.1:1]; 

% x_new = 8 + 0*t;
% z_new = 3 + 0*t + .5*(-.5)*t; 
% y_new = zeros(1, length(t));


dt = .1;

p1 = [0 0 4 8]; % List of all X points in Sequence
p2 = [0 0 0 0]; % List of all Y points in Sequence
p3 = [0 2 2 3]; % List of all Z points in Sequence

x0 = 8; 
z0 = 3; 


theta = atan2(z0,x0);

for i = 1:1:length(t)
    
plot3(p1, p2, p3)
xlabel("X Label")
ylabel("Y Label")
zlabel("Z Label")

    pause(0.5);

    
    x_new = x0 - 4.123*sin(theta)*dt;
    z_new = z0 - 4.123*cos(theta)*dt;
    y_new = 0;
    
    x0 = x_new;
    z0 = y_new; 
    
    
    theta = theta + 1*dt;

    p1(1,4)= x_new;
    p2(1,4)= y_new;
    p3(1,4)= z_new;
    
    
    
end 

