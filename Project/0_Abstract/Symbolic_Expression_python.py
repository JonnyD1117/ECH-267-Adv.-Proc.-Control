from sympy import symbols
from sympy import *
from sympy.physics.vector import dynamicsymbols


init_printing()

q1, q2, q3, m1, m2, m3, L1, L2, L3, g, I1, I2, I3,t, Tau1, Tau2, Tau3 = symbols('q1 q2 q3 m1 m2 m3 L1 L2 L3 g I1 I2 I3 t Tau1 Tau2 Tau3')

q1 = dynamicsymbols('q1')
q2 = dynamicsymbols('q2')
q3 = dynamicsymbols('q3')

q1_dot = diff(q1, Symbol('t'))
q2_dot = diff(q2, Symbol('t'))
q3_dot = diff(q3, Symbol('t'))


x1 = 0
x2  =.5*L2*cos(q1)
x3 = L2*cos(q1) + .5*L3*cos(q2)

y1 = L1
y2 = L1 + .5*L2*sin(q1)
y3 = L1 + L2*sin(q1) + .5*L3*sin(q2)

x1_dot = diff(x1, Symbol('t'))
x2_dot = diff(x2, Symbol('t'))
x3_dot = diff(x3, Symbol('t'))
y1_dot = diff(y1, Symbol('t'))
y2_dot = diff(y2, Symbol('t'))
y3_dot = diff(y3, Symbol('t'))


v1_sqr = x1_dot**2 + y1_dot**2
v2_sqr = x2_dot**2 + y2_dot**2
v3_sqr = x3_dot**2 + y3_dot**2

# Potential Energy
U1 = m1*g*y1
U2 = m2*g*y2
U3 = m3*g*y3

U = U1 + U1 + U1

# Kinetic Energy
K1 = .5*m1*v1_sqr + .5*I1*(q1_dot)**2
K2 = .5*m2*v2_sqr + .5*I2*(q2_dot)**2
K3 = .5*m3*v3_sqr + .5*I3*(q3_dot)**2

T = K1 + K2 + K3

dU_dq1 = diff(U, q1)
dU_dq2 = diff(U, q2)
dU_dq3 = diff(U, q3)

dT_dq1 = diff(T, q1)
dT_dq2 = diff(T, q2)
dT_dq3 = diff(T, q3)

dT_d_q1_dot = diff(T, q1_dot)
dT_d_q2_dot = diff(T, q2_dot)
dT_d_q3_dot = diff(T, q3_dot)

dt_T_q1_dot = diff(dT_d_q1_dot, t)
dt_T_q2_dot = diff(dT_d_q2_dot, t)
dt_T_q3_dot = diff(dT_d_q3_dot, t)



eqn1 = -dt_T_q1_dot + dT_dq1 - dU_dq1 + Tau1
eqn2 = -dt_T_q2_dot + dT_dq2 - dU_dq2 + Tau2
eqn3 = -dt_T_q3_dot + dT_dq3 - dU_dq3 + Tau3



eqn1 = simplify(Eq(dt_T_q1_dot - dT_dq1 + dU_dq1, Tau1))
eqn1


eqn2 = Eq(dt_T_q2_dot - dT_dq2 + dU_dq2, Tau2)
eqn2


eqn3 = Eq(dt_T_q3_dot - dT_dq3 + dU_dq3, Tau3)
eqn3


Tau1 = I1*q1_dot_dot + 0.25*L2**2*m2*q1_dot_dot + 1.0*L2**2*m3*q1_dot_dot + 0.5*L2*L3*m3*sin(q1(t) - q2(t))*q1_dot**2 + 0.5*L2*L3*m3*cos(q1(t) - q2(t))*q1_dot_dot)
Tau2 = I2*q2_dot_dot - 0.5*m3*(-1.0*L3*(-L2*sin(q1(t))*q1_dot - 0.5*L3*sin(q2(t))*q2_dot)*cos(q2(t))*q2_dot - 1.0*L3*(L2*cos(q1(t))*q1_dot + 0.5*L3*cos(q2(t))*q2_dot)*sin(q2(t))*q2_dot) + 0.5*m3*(-1.0*L3*(-L2*sin(q1(t))*q1_dot - 0.5*L3*sin(q2(t))*q2_dot)*cos(q2(t))*q2_dot - 1.0*L3*(L2*cos(q1(t))*q1_dot + 0.5*L3*cos(q2(t))*q2_dot)*sin(q2(t))*q2_dot + 1.0*L3*(-L2*sin(q1(t))*q1_dot**2 + L2*cos(q1(t))*q1_dot_dot - 0.5*L3*sin(q2(t))*q2_dot**2 + 0.5*L3*cos(q2(t))*q2_dot_dot)*cos(q2(t)) - 1.0*L3*(-L2*sin(q1(t))*q1_dot_dot - L2*cos(q1(t))*q1_dot**2 - 0.5*L3*sin(q2(t))*q2_dot_dot - 0.5*L3*cos(q2(t))*q2_dot**2)*sin(q2(t)))
Tau3 = I3*q3_dot_dot