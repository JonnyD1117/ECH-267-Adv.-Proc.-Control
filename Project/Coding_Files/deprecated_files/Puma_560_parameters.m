%% PUMA 560 Robot Paramters 


%% DH Parameters 
alpha = [0, -pi/2, 0];
A = [0, 0, .4318]; 
D = [0, .2435, -.0934]; 


%% Link Masses (Kg)
% m1; 
m2 = 17.40; 
m3 = 4.80; 


%% Link Moments of Inertia about CG (Kg m^2)
Izz_1 = 0.350; 

Ixx_2 = .130;
Iyy_2 = .524;
Izz_2 = .539;

Ixx_3 = 66.0e-3;
Iyy_3 = 12.5e-3;
Izz_3 = 86e-3;

S2 = [0.068;.006; -.016] ; 
S3 = [0; -.070; .014]; 