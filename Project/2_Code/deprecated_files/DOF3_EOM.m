




q1= 0;
q2= 0;
q3= 0;


q1_dot = 0;
q2_dot = 0;
q3_dot = 0;

q1_dot_dot= 0; 
q2_dot_dot= 0; 
q3_dot_dot= 0; 


Tau1 = I1*q1_dot_dot + 0.25*L2^2*m2*q1_dot_dot + 1.0*L2^2*m3*q1_dot_dot + 0.5*L2*L3*m3*sin(q1 - q2)*q1_dot^2 + 0.5*L2*L3*m3*cos(q1 - q2)*q1_dot_dot;
Tau2 = I2*q2_dot_dot - 0.5*m3*(-1.0*L3*(-L2*sin(q1)*q1_dot - 0.5*L3*sin(q2)*q2_dot)*cos(q2)*q2_dot - 1.0*L3*(L2*cos(q1)*q1_dot + 0.5*L3*cos(q2)*q2_dot)*sin(q2)*q2_dot)...
    + 0.5*m3*(-1.0*L3*(-L2*sin(q1)*q1_dot - 0.5*L3*sin(q2)*q2_dot)*cos(q2)*q2_dot - 1.0*L3*(L2*cos(q1)*q1_dot + 0.5*L3*cos(q2)*q2_dot)*sin(q2)*q2_dot + 1.0*L3*(-L2*sin(q1)*q1_dot^2 ...
    + L2*cos(q1)*q1_dot_dot - 0.5*L3*sin(q2)*q2_dot^2 + 0.5*L3*cos(q2)*q2_dot_dot)*cos(q2) - 1.0*L3*(-L2*sin(q1)*q1_dot_dot - L2*cos(q1)*q1_dot^2 - 0.5*L3*sin(q2)*q2_dot_dot - ...
    0.5*L3*cos(q2)*q2_dot^2)*sin(q2));
Tau3 = I3*q3_dot_dot; 


Tau1 
Tau2 
Tau3