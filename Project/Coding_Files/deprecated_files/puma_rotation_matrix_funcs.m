%% PUMA Robot Homogeneous & Rotation Matrice Transforms

syms q_d1 q_d2 q_d3




% 1-0 Reference Transform 
T_01  = T_h(0,0, 0, q_d1);
R_01 = T_01(1:3,1:3);
R_01_simp = simplify(R_01);


% 2-1 Reference Transform 
T_12  = T_h(-pi/2, 0, .2435, q_d2);
R_12 = T_12(1:3,1:3);


% 3-2 Reference Transform 
T_23  = T_h(0,.4318, -.0934, q_d3);
R_23 = T_23(1:3, 1:3);


% 2-0 Reference Transform 

T_02 = T_01*T_12;
R_02 = T_02(1:3, 1:3)
R_02_simp = simplify(R_02);

% 3-0 Reference Transform 
T_03 = T_01*T_12*T_23;
R_03 = T_03(1:3, 1:3);
R_03_simp = simplify(R_03);

%% Test Rotation Functions 


rot_1 = R_01_func(0);
rot_2 = R_02_func(0,0);
rot_3 = R_03_func(0,0,90);



%% Rotation Functions 

function out = R_01_func(q_1)

out =  [cos(q_1), -sin(q_1), 0; ...
        sin(q_1),  cos(q_1), 0; ...
           0,       0,       1]; 

end 

function out = R_02_func(q_1, q_2)


out = [ cos(q_1)*cos(q_2) - (4967757600021511*sin(q_1)*sin(q_2))/81129638414606681695789005144064, - cos(q_1)*sin(q_2) - (4967757600021511*cos(q_2)*sin(q_1))/81129638414606681695789005144064,                                        -sin(q_1); ...
       (4967757600021511*cos(q_1)*sin(q_2))/81129638414606681695789005144064 + cos(q_2)*sin(q_1),   (4967757600021511*cos(q_1)*cos(q_2))/81129638414606681695789005144064 - sin(q_1)*sin(q_2),                                         cos(q_1); ...
                                            -sin(q_2),                                                              -cos(q_2),           4967757600021511/81129638414606681695789005144064];
 


end 

function out = R_03_func(q_1, q_2, q_3)

out = [ cos(q_3)*(cos(q_1)*cos(q_2) - (4967757600021511*sin(q_1)*sin(q_2))/81129638414606681695789005144064) - sin(q_3)*(cos(q_1)*sin(q_2) + (4967757600021511*cos(q_2)*sin(q_1))/81129638414606681695789005144064), - sin(q_3)*(cos(q_1)*cos(q_2) - (4967757600021511*sin(q_1)*sin(q_2))/81129638414606681695789005144064) - cos(q_3)*(cos(q_1)*sin(q_2) + (4967757600021511*cos(q_2)*sin(q_1))/81129638414606681695789005144064),                                        -sin(q_1); ...
 sin(q_3)*((4967757600021511*cos(q_1)*cos(q_2))/81129638414606681695789005144064 - sin(q_1)*sin(q_2)) + cos(q_3)*((4967757600021511*cos(q_1)*sin(q_2))/81129638414606681695789005144064 + cos(q_2)*sin(q_1)),   cos(q_3)*((4967757600021511*cos(q_1)*cos(q_2))/81129638414606681695789005144064 - sin(q_1)*sin(q_2)) - sin(q_3)*((4967757600021511*cos(q_1)*sin(q_2))/81129638414606681695789005144064 + cos(q_2)*sin(q_1)),                                         cos(q_1); ...
                                                                                                                                                                         - cos(q_2)*sin(q_3) - cos(q_3)*sin(q_2),                                                                                                                                                                               sin(q_2)*sin(q_3) - cos(q_2)*cos(q_3), 4967757600021511/81129638414606681695789005144064];
 


end 
%% Functions 

function output = T_Rz(theta)

R_z = [cos(theta) -sin(theta) 0; 
       sin(theta) cos(theta) 0; 
           0          0      1 ]; 
       
    output = [R_z [0;0;0]; 
              [0 0 0] 1];  
end 

function output = T_z(dist)

output = [1 0 0 0; 
          0 1 0 0; 
          0 0 1 dist;
          0 0 0 1];
end 

function output = T_Rx(theta)

R_x = [ 1 0 0 ;
        0 cos(theta) -sin(theta); 
       0 sin(theta) cos(theta) ]; 
       
    output = [R_x [0;0;0]; 
              [0 0 0] 1];  
end 

function output = T_x(dist)

output = [1 0 0 dist; 
          0 1 0 0; 
          0 0 1 0;
          0 0 0 1];
end 

function output = T_h(alpha, a, d, theta)


output = T_Rx(alpha)*T_x(a)*T_Rz(theta)*T_z(d);

end 