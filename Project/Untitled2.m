clear all 

close all 


syms m1 m2 m3 L1 L2 L3 q1 q2 q3 q1_dot q2_dot g t1 t2

tau = [t1;t2];

eps = .1 ;

m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(q2) + (L2)^2) + eps, m2*(L1*L2*cos(q2)+ (L2)^2); ...
     m2*(L1*L2*cos(q2) + (L2)^2), m2*(L2)^2+eps];

v = [-m2*L1*L2*sin(q2)*(2*q1_dot*q2_dot + (q2_dot)^2);...
     m2*L1*L2*(q1_dot)^2*sin(q2)] ; 

g = [(m1 + m2)*L1*g*cos(q1)+ m2*g*L2*cos(q1+ q2); ...
     m2*g*L2*cos(q1+q2)];
 

f = [q1_dot*.25;q2_dot*.25];
 
 q_ddot = m\(tau - v -g - f);
 
 
 f1 = q1_dot; 
 f2 = q2_dot;
 f3 = q_ddot(1);
 f4 = q_ddot(2);
 
 d_f1_dx1 = diff(f1,q1);
 d_f1_dx2 = diff(f1,q2);
 d_f1_dx3 = diff(f1,q1_dot);
 d_f1_dx4 = diff(f1,q2_dot);
 

 d_f2_dx1 = diff(f2,q1);
 d_f2_dx2 = diff(f2,q2);
 d_f2_dx3 = diff(f2,q1_dot);
 d_f2_dx4 = diff(f2,q2_dot);
 
 
 d_f3_dx1 = diff(f3,q1);
 d_f3_dx2 = diff(f3,q2);
 d_f3_dx3 = diff(f3,q1_dot);
 d_f3_dx4 = diff(f3,q2_dot);


 d_f4_dx1 = diff(f4,q1);
 d_f4_dx2 = diff(f4,q2);
 d_f4_dx3 = diff(f4,q1_dot);
 d_f4_dx4 = diff(f4,q2_dot);
 
 
 
 A = [d_f1_dx1, d_f1_dx2, d_f1_dx3, d_f1_dx4;...
      d_f2_dx1, d_f2_dx2, d_f2_dx3, d_f2_dx4;...
      d_f3_dx1, d_f3_dx2, d_f3_dx3, d_f3_dx4; ...
      d_f4_dx1, d_f4_dx2, d_f4_dx3, d_f4_dx4];
  
  A = simplify(A); 
  
  
  
 d_f1_du1 = diff(f1,t1);
 d_f1_du2 = diff(f1,t2);

 d_f2_du1 = diff(f2,t1);
 d_f2_du2 = diff(f2,t2);
 
 d_f3_du1 = diff(f3,t1);
 d_f3_du2 = diff(f3,t2);

 d_f4_du1 = diff(f4,t1);
 d_f4_du2 = diff(f4,t2);
 
 B = [d_f1_du1,d_f1_du2;...
      d_f2_du1,d_f2_du2;...
      d_f3_du1,d_f3_du2;...
      d_f4_du1,d_f4_du2];
  
  C = eye(4); 
  
   dt = .1; 
  
%   sys = ss(A,B,C, [0;0;0;0] )
 
 
 