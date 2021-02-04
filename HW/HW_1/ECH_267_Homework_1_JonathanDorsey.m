


%% ECH 267 Homework #1 Matlab Work 


%% Problem #6

% 6.1 
 [x1, x2] = meshgrid(-3:0.2:3, -1.5:0.2:1.5);
 
 x1_dot = x2; 
 x2_dot = -x1 +(x1).^3/6 -x2;
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)

 plot( sqrt(6), 0, 'r*')
  plot( -sqrt(6), 0, 'r*')
   plot( 0, 0, 'r*')

%% Problem 7

% Matrices
A_1 = [0 1; -2 -3];
A_2 = [0 -1 ;1 2];
A_3 = [1 1; 0 -1];
A_4 = [1 5 ;-1 -1];
A_5 = [2 -1 ;2 0];

% Eigen Values & Vectors
[v1, d1] = eig(A_1);
[v2, d2] = eig(A_2);
[v3, d3] = eig(A_3);
[v4, d4] = eig(A_4);
[v5, d5] = eig(A_5);

% Modal Matrices

M_1 = v1;
M_1_inv = inv(M_1) ;
T_1 = M_1\A_1*M_1 ;

M_2 = v2 ;
M_2_inv =inv(M_2); 
T_2 = M_2\A_2*M_2;


M_3 = v3;
M_3_inv = inv(M_3) ;
T_3 = M_3\A_3*M_3;


M_4 = v4;
M_4_inv = inv(M_4);
T_4 = M_4\A_4*M_4;


M_5 = v5;
M_5_inv = inv(M_5);
T_5 = M_5\A_5*M_5;


% Modal Matrice Check


% Generate Phase Portraits in both Modal & original Coordinates

% Problem 7 C Plot Generation. 

%Part 1 
[x1, x2] = meshgrid(-1.5:0.2:1.5, -1:0.2:1);
[z1, z2] = meshgrid(-1:0.2:1, -1:0.2:1);
 
% Normal Coordinates 
x1_dot = A_1(1,1)*x1 + A_1(1,2)*x2; 
x2_dot = A_1(2,1)*x1 + A_1(2,2)*x2; 
% Modal Coordinates
z1_dot = T_1(1,1)*z1 + T_1(1,2)*z2; 
z2_dot = T_1(2,1)*z1 + T_1(2,2)*z2; 

 
%  close all
%  hold on
%  figure(1)
%  quiver(x1, x2, x1_dot, x2_dot)
%  title(' Problem 7.C.i Phase Portrait: Original Coordinates ')
%  plot( 0, 0, 'r*')
 
 
 
%  figure(2)
%  hold on 
%  quiver(z1, z2, z1_dot, z2_dot)
%  title(' Problem 7.C.i Phase Portrait: Modal Coordinates ')
%  plot( 0, 0, 'r*')
%  
 

%Part 2 

[x1, x2] = meshgrid(-1.5:0.2:1.5, -1:0.2:1);
[z1, z2] = meshgrid(-1:0.2:1, -1:0.2:1);
 
% Normal Coordinates 
x1_dot = A_2(1,1)*x1 + A_2(1,2)*x2; 
x2_dot = A_2(2,1)*x1 + A_2(2,2)*x2; 
% Modal Coordinates
z1_dot = T_2(1,1)*z1 + T_2(1,2)*z2; 
z2_dot = T_2(2,1)*z1 + T_2(2,2)*z2; 

  figure(3)

 hold on
 quiver(x1, x2, x1_dot, x2_dot)
 plot( 0, 0, 'r*')
  title(' Problem 7.C.ii Phase Portrait: Original Coordinates ')

   figure(4)

 hold on 

 quiver(z1, z2, z1_dot, z2_dot)
 plot( 0, 0, 'r*')
 title(' Problem 7.C.ii Phase Portrait: Modal Coordinates ')

%Part 3 

[x1, x2] = meshgrid(-1.5:0.2:1.5, -1:0.2:1);
[z1, z2] = meshgrid(-1:0.2:1, -1:0.2:1);
 
% Normal Coordinates 
x1_dot = A_3(1,1)*x1 + A_3(1,2)*x2; 
x2_dot = A_3(2,1)*x1 + A_3(2,2)*x2; 
% Modal Coordinates
z1_dot = T_3(1,1)*z1 + T_3(1,2)*z2; 
z2_dot = T_3(2,1)*z1 + T_3(2,2)*z2; 

  figure(5)

 hold on
 quiver(x1, x2, x1_dot, x2_dot)
 plot( 0, 0, 'r*')
  title(' Problem 7.C.iii Phase Portrait: Original Coordinates ')

   figure(6)

 hold on 

 quiver(z1, z2, z1_dot, z2_dot)
 plot( 0, 0, 'r*')
  title(' Problem 7.C.iii Phase Portrait: Modal Coordinates ')

 
%Part 4 

[x1, x2] = meshgrid(-1.5:0.2:1.5, -1:0.2:1);
[z1, z2] = meshgrid(-1:0.2:1, -1:0.2:1);
 
% Normal Coordinates 
x1_dot = A_4(1,1)*x1 + A_4(1,2)*x2; 
x2_dot = A_4(2,1)*x1 + A_4(2,2)*x2; 
% Modal Coordinates
z1_dot = T_4(1,1)*z1 + T_4(1,2)*z2; 
z2_dot = T_4(2,1)*z1 + T_4(2,2)*z2; 

  figure(7)

 hold on
 quiver(x1, x2, x1_dot, x2_dot)
 plot( 0, 0, 'r*')
  title(' Problem 7.C.iv Phase Portrait: Original Coordinates ')

   figure(8)

 hold on 

 quiver(z1, z2, z1_dot, z2_dot)
 plot( 0, 0, 'r*')
  title(' Problem 7.C.iv Phase Portrait: Modal Coordinates ')

 
%Part 5 
[x1, x2] = meshgrid(-1.5:0.2:1.5, -1:0.2:1);
[z1, z2] = meshgrid(-1:0.2:1, -1:0.2:1);
 
% Normal Coordinates 
x1_dot = A_5(1,1)*x1 + A_5(1,2)*x2; 
x2_dot = A_5(2,1)*x1 + A_5(2,2)*x2; 
% Modal Coordinates
z1_dot = T_5(1,1)*z1 + T_5(1,2)*z2; 
z2_dot = T_5(2,1)*z1 + T_5(2,2)*z2; 

  figure(9)

 hold on
 quiver(x1, x2, x1_dot, x2_dot)
 plot( 0, 0, 'r*')
  title(' Problem 7.C.v Phase Portrait: Original Coordinates ')

 figure(10)
 hold on 
 quiver(z1, z2, z1_dot, z2_dot)
 plot( 0, 0, 'r*')
 title(' Problem 7.C.v Phase Portrait: Modal Coordinates ')


%% Problem #8


% [x1, x2] = meshgrid(-2.:0.2:2., -2:0.2:2);
% [z1, z2] = meshgrid(-2:0.2:2, -2:0.2:2);
syms x1 x2

% Normal Coordinates 
f1 = - x1 -(x2/(log(sqrt(x1^2 + x2^2))))  ;
f2 = - x2 + (x1/(log(sqrt(x1^2 + x2^2))))  ;
 
df1_dx1 = diff(f1, x1);
df1_dx2 = diff(f1, x2);

df2_dx1 = diff(f2, x1);
df2_dx2 = diff(f2, x2);

latex_1 = latex(df1_dx1);
latex_2 = latex(df1_dx2);
latex_3 = latex(df2_dx1);
latex_4 = latex(df2_dx2);


x1_eq = 0.1; x2_eq =0.1; 

df1_dx1_out = subs(df1_dx1, [x1, x2], [x1_eq, x2_eq]);
df1_dx2_out = subs(df1_dx2, [x1, x2], [x1_eq, x2_eq]);
df2_dx1_out = subs(df2_dx1, [x1, x2], [x1_eq, x2_eq]);
df2_dx2_out = subs(df2_dx2, [x1, x2], [x1_eq, x2_eq]);



df1_dx1_val = eval(df1_dx1_out);
df1_dx2_val = eval(df1_dx2_out);
df2_dx1_val = eval(df2_dx1_out);
df2_dx2_val = eval(df2_dx2_out);

A_j = [df1_dx1_val df1_dx2_val; df2_dx1_val df2_dx2_val];


[vect, val] = eig(A_j);


[x1, x2] = meshgrid(-1:0.2:1, -1:0.2:1);
 
% Normal Coordinates 
x1_dot = A_j(1,1)*x1 + A_j(1,2)*x2; 
x2_dot = A_j(2,1)*x1 + A_j(2,2)*x2; 


figure(9)

hold on
quiver(x1, x2, x1_dot, x2_dot)
plot( 0, 0, 'r*')
title(' Problem 8 Linearized Phase Portrait')




clear x1 x2 
[x1, x2] = meshgrid(.1:0.2:2, .1:0.2:2);
 
x1_dot = - x1 -(x2/(log(sqrt(x1^2 + x2^2))))  ;
x2_dot = - x2 + (x1/(log(sqrt(x1^2 + x2^2))))  ;

figure(10)

hold on
quiver(x1, x2, x1_dot, x2_dot)
plot( 0, 0, 'r*')
title(' Problem 8 Nonlinear Phase Portrait')

%%
 

r_list = [.1, .25, .5, .75, 1, 1.5]; 
theta_list= [ .1, .25*pi, .5*pi, pi,(3/2)*pi, 2*pi ];



[r, theta] = meshgrid(r_list, theta_list);

x1 = r.*cos(theta); 
x2 = r.*sin(theta); 

d_x1_list = - x1 -(x2)./(log(sqrt(x1.^2 + x2.^2))) ;
d_x2_list = - x2 +(x1)./(log(sqrt(x1.^2 + x2.^2))) ;

figure(11)
hold on
quiver(x1, x2, d_x1_list, d_x2_list)
plot( 0, 0, 'r*')
title(' Problem 8 Nonlinear Phase Portrait: Polar Form')








%% Problem 9.1

[x1, x2] = meshgrid(-2:0.2:2, -2:0.2:2);
 
 x1_dot = x2; 
 x2_dot = x1-2*atan(x1+x2);
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)


   plot( 0, 0, 'r*')
   title("Problem 9.1 Phase Portrait")

%% Problem 9.2

[x1, x2] = meshgrid(-2.5:0.2:2.5, -2.5:0.2:2.5);
 
 x1_dot = x2; 
 x2_dot = -x1 + x2.*(1 - 3*x1.^2 -2*x2.^2 );
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)

 plot( 0, 0, 'r*')

      title("Problem 9.2 Phase Portrait")

   
   %% Problem 9.3
[x1, x2] = meshgrid(-2.5:0.2:2.5, -2.5:0.2:2.5);
 
 x1_dot = x1 - x1*x2; 
 x2_dot = 2.*x1.^2 - x2;
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)

 plot( 0, 0, 'r*')

      title("Problem 9.3 Phase Portrait")
   %% Problem 9.4

[x1, x2] = meshgrid(-2.5:0.2:2.5, -2.5:0.2:2.5);
 
 x1_dot = x1 + x2 - x1*(abs(x1) + abs(x2));
 x2_dot = -2.*x1 + x2 - x2*(abs(x1) + abs(x2));
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)

 plot( 0, 0, 'r*')

      title("Problem 9.4 Phase Portrait")
%% Problem #3 


k1 = 2;
k2 =1; 

A = [-k1 0 0; k1 -k2 0; 0 k2 0];
B = [0;0;0];
% C = [1 0 0];
C = eye(3);

D = [0];

ss_system = ss(A, B,C, D);

x0 = [1;0;0]; 


initial(ss_system,x0);


