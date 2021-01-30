


%% ECH 267 Homework #1 Matlab Work 


%% Problem #1
%% Problem #2
%% Problem #3
%% Problem #4
%% Problem #5
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
M_1_inv = inv(M_1); 

M_2 = v2 ;
M_2_inv =inv(M_2) ; 

M_3 = v3;
M_3_inv = inv(M_3); 

M_4 = v4 ;
M_4_inv = inv(M_4); 

M_5 = v5;
M_5_inv = inv(M_5); 

% Modal Matrice Check


check1 = M_1_inv*A_1*M_1
check2 = M_2_inv*A_2*M_2
check3 = M_3_inv*A_3*M_3
check4 = M_4_inv*A_4*M_4
check5 = M_5_inv*A_5*M_5


% Generate Phase Portraits in both Modal & original Coordinates

%% Problem #8


%% Problem 9.1

[x1, x2] = meshgrid(-2:0.2:2, -2:0.2:2);
 
 x1_dot = x2; 
 x2_dot = x1-2*atan(x1+x2);
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)


   plot( 0, 0, 'r*')

%% Problem 9.2

[x1, x2] = meshgrid(-8:0.2:10, -6:0.2:7);
 
 x1_dot = x2; 
 x2_dot = -10*sin(x1) - x2;
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)

 plot( -6, 0, 'r*')
  plot( 6, 0, 'r*')
   plot( 0, 0, 'r*')
   
   %% Problem 9.3

[x1, x2] = meshgrid(-8:0.2:10, -6:0.2:7);
 
 x1_dot = x2; 
 x2_dot = -10*sin(x1) - x2;
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)

 plot( -6, 0, 'r*')
  plot( 6, 0, 'r*')
   plot( 0, 0, 'r*')
   %% Problem 9.4

[x1, x2] = meshgrid(-8:0.2:10, -6:0.2:7);
 
 x1_dot = x2; 
 x2_dot = -10*sin(x1) - x2;
 
 close all
 hold on
 quiver(x1, x2, x1_dot, x2_dot)

 plot( -6, 0, 'r*')
  plot( 6, 0, 'r*')
   plot( 0, 0, 'r*')
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


