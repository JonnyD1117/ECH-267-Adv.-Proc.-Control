% clear all 
clearvars -except x_list
% close all 


syms m1 m2 m3 L1 L2 L3 q1 q2 q3 q1_dot q2_dot g_val t1 t2

tau = [t1;t2];

eps = .1 ;

m = [m1*(L1)^2 + m2*(L1^2 + 2*L1*L2*cos(q2) + (L2)^2) + eps, m2*(L1*L2*cos(q2)+ (L2)^2); ...
     m2*(L1*L2*cos(q2) + (L2)^2), m2*(L2)^2+eps];

v = [-m2*L1*L2*sin(q2)*(2*q1_dot*q2_dot + (q2_dot)^2);...
     m2*L1*L2*(q1_dot)^2*sin(q2)] ; 

g = [(m1 + m2)*L1*g_val*cos(q1)+ m2*g_val*L2*cos(q1+ q2); ...
     m2*g_val*L2*cos(q1+q2)];
 

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
%% PID Controller 

x0 = [0,0,0,0];         % Equilibrium Point
states_sym = [q1, q2, q1_dot, q2_dot]; % State Symbols
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
params_val = [.25,.25, .5, .25, 9.84];


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);

% Linear State Space System
Gp_ss = ss(A_mat, B_mat, C, zeros(4,2));

% Transfer function System 
Gp_tf = tf(Gp_ss);

% Plant TF: wrt Input #1
Gp_x1_u1 = Gp_tf(1,1);
Gp_x2_u1 = Gp_tf(2,1);
Gp_x3_u1 = Gp_tf(3,1);
Gp_x4_u1 = Gp_tf(4,1);

% Plant TF: wrt Input #2
Gp_x1_u2 = Gp_tf(1,2);
Gp_x2_u2 = Gp_tf(2,2);
Gp_x3_u2 = Gp_tf(3,2);
Gp_x4_u2 = Gp_tf(4,2);

% PID Gains and Controller TF
k_i = 5;
k_p = 20; 
k_d = 5; 

G_c = tf([k_d, k_p, k_i],[1, 0]);

% Feedforward Gains: Input #1
FF_x1_u1 =G_c*Gp_x1_u1; 
FF_x2_u1 =G_c*Gp_x2_u1; 
FF_x3_u1 =G_c*Gp_x3_u1; 
FF_x4_u1 =G_c*Gp_x4_u1; 

% Feedforward Gains: Input #2
FF_x1_u2 =G_c*Gp_x1_u2; 
FF_x2_u2 =G_c*Gp_x2_u2; 
FF_x3_u2 =G_c*Gp_x3_u2; 
FF_x4_u2 =G_c*Gp_x4_u2; 

% Closed Loop Gains: Input #1
G_cl_x1_u1 = (FF_x1_u1)/(1+FF_x1_u1);
G_cl_x2_u1 = (FF_x2_u1)/(1+FF_x2_u1);
G_cl_x3_u1 = (FF_x3_u1)/(1+FF_x3_u1);
G_cl_x4_u1 = (FF_x4_u1)/(1+FF_x4_u1);

% Closed Loop Gains: Input #2
G_cl_x1_u2 = (FF_x1_u2)/(1+FF_x1_u2);
G_cl_x2_u2 = (FF_x2_u2)/(1+FF_x2_u2);
G_cl_x3_u2 = (FF_x3_u2)/(1+FF_x3_u2);
G_cl_x4_u2 = (FF_x4_u2)/(1+FF_x4_u2);



% % Bode Plots
% figure()
% bode(G_cl_x1_u1)
% figure()
% nyquist(G_cl_x1_u1)

% figure()
% hold on 
% step(G_cl_x1_u1)
% step(G_cl_x2_u1)
% step(G_cl_x3_u1)
% step(G_cl_x4_u1)
% 
% 
% legend('q1','q2','q1_dot','q2_dot')

%% 
close all
x_list  =[0.436332312998582,0.436332312998582,0.820740251679684,1.17012419486637,1.46715299960859,1.71874300527458,1.93206754943125,2.11315688994308,2.26701502505225,2.39781203288825,2.50904325001313,2.60365317586080,2.68413174918613,2.75259017252834,2.81082169740738,2.86035118994227,2.90247612295057,2.93830084181358,2.96876542510260,2.99467011597050,3.01669607083396,3.03542301467511,3.05134428035439,3.06487962620156,3.07638616193970,3.08616766172851,3.09448250113355,3.10155041983562,3.10755828237970,3.11266498418282,3.11700562861715,3.12069508267713,3.12383100306408,3.12649641109811,3.12876188337353,3.13068741523860,3.13232400576934,3.13371500571801,3.13489726377689,3.13590210125686,3.13675614080927,3.13748201100644,3.13809894534730,3.13862329148560,3.13906894412051,3.13944771298175,3.13976963563260,3.14004324335959,3.14027578718050,3.14047342994948,3.14064140964262,3.14078417814597,3.14090551921991,3.14100864876358,3.14109630003440,3.14117079607983,3.14123411129984,3.14128792377103,3.14133365971847,3.14137253131364,3.14140556879986,3.14143364779641,3.14145751250502,3.14147779543343,3.14149503415893,3.14150968557594,3.14152213800532,3.14153272148623,3.14154171652339,3.14154936152152;
          0.785398163397448,0.785398163397448,0.552454697426174,0.352560752853411,0.182510181855351,0.0377353948107595,-0.0854796071542897,-0.190304300941552,-0.279458097821766,-0.355269758222179,-0.419728712688828,-0.474531389533293,-0.521122588216915,-0.560731854309418,-0.594405067108329,-0.623031639656576,-0.647367815332755,-0.668056556292761,-0.685644489984407,-0.700596333301256,-0.713307162538916,-0.724112847618225,-0.733298923833359,-0.741108134543664,-0.747746843712304,-0.753390487566873,-0.758188209350971,-0.762266799576903,-0.765734045847164,-0.768681580715570,-0.771187302800086,-0.773317435088686,-0.775128274798811,-0.776667681006375,-0.777976339336280,-0.779088837119830,-0.780034577419780,-0.780838556068768,-0.781522022249275,-0.782103040067524,-0.782596965958766,-0.783016854538094,-0.783373803620763,-0.783677247528949,-0.783935206435640,-0.784154498334812,-0.784340919239527,-0.784499396370095,-0.784634118380687,-0.784748646066050,-0.784846006474121,-0.784928772911826,-0.784999132958541,-0.785058946284749,-0.785109793804031,-0.785153019457440,-0.785189765734628,-0.785221003870526,-0.785247559515702,-0.785270134558851,-0.785289325678191,-0.785305640112092,-0.785319509065754,-0.785331299108286,-0.785341321861405,-0.785349842235842,-0.785357085433144,-0.785363242897933,-0.785368477377947,-0.785372927225595;
          0,3.84407938681102,3.49383943186684,2.97028804742218,2.51590005665992,2.13324544156675,1.81089340511828,1.53858135109169,1.30797007836004,1.11231217124876,0.946099258476672,0.804785733253333,0.684584233422141,0.582315248790312,0.495294925348982,0.421249330082945,0.358247188630143,0.304645832890189,0.259046908678966,0.220259548634600,0.187269438411493,0.159212656792827,0.135353458471720,0.115065357381375,0.0978149978881117,0.0831483940503659,0.0706791870207202,0.0600786254408117,0.0510670180311530,0.0434064443433517,0.0368945405997629,0.0313592038694784,0.0266540803403916,0.0226547227541931,0.0192553186506526,0.0163659053074096,0.0139099994867408,0.0118225805887117,0.0100483747997542,0.00854039552404483,0.00725870197170560,0.00616934340865937,0.00524346138302436,0.00445652634901646,0.00378768861239831,0.00321922650853071,0.00273607726993446,0.00232543820912495,0.00197642768972008,0.00167979693146369,0.00142768503350101,0.00121341073939470,0.00103129543667344,0.000876512708233365,0.000744960454298925,0.000633152200119333,0.000538124711851318,0.000457359474398216,0.000388715951748604,0.000330374862141969,0.000280789965516890,0.000238647086041280,0.000202829284103676,0.000172387255024682,0.000146514170163033,0.000124524293785107,0.000105834809062383,8.99503716209979e-05,7.64499813015662e-05,6.49758242260821e-05;
          0,-2.32943465971274,-1.99893944572763,-1.70050570998060,-1.44774787044592,-1.23215001965049,-1.04824693787263,-0.891537968802136,-0.758116604004130,-0.644589544666489,-0.548026768444649,-0.465911986836222,-0.396092660925028,-0.336732127989116,-0.286265725482466,-0.243361756761796,-0.206887409600060,-0.175879336916454,-0.149518433168486,-0.127108292376606,-0.108056850793091,-0.0918607621513327,-0.0780921071030555,-0.0663870916864002,-0.0564364385456918,-0.0479772178409722,-0.0407859022593234,-0.0346724627026106,-0.0294753486840607,-0.0250572208451589,-0.0213013228860031,-0.0181083971012494,-0.0153940620756413,-0.0130865832990471,-0.0111249778354967,-0.00945740299949986,-0.00803978648988660,-0.00683466180506548,-0.00581017818249237,-0.00493925891241902,-0.00419888579327957,-0.00356949082668881,-0.00303443908185787,-0.00257958906691529,-0.00219291899171759,-0.00186420904714938,-0.00158477130567648,-0.00134722010592017,-0.00114527685362927,-0.000973604080710030,-0.000827664377057545,-0.000703600467144834,-0.000598133262078738,-0.000508475192818976,-0.000432256534099527,-0.000367462771875419,-0.000312381358977140,-0.000265556451765155,-0.000225750431492365,-0.000191911193399951,-0.000163144339005885,-0.000138689536622746,-0.000117900425320552,-0.000100227531188105,-8.52037443678623e-05,-7.24319730189044e-05,-6.15746478935270e-05,-5.23448001365405e-05,-4.44984764801775e-05,-3.78282913748796e-05];

dt = .1;
t = 0:dt:(7-dt);

u1 = x_list(1,:);
u2 = x_list(2,:);


figure() 
hold on
y1 = lsim(G_cl_x1_u1,u1,t);
y2 = lsim(G_cl_x2_u2,u2,t);
% lsim(G_cl_x3_u1,u1,t)
% lsim(G_cl_x4_u2,u2,t)
legend('q1','q2')


figure() 
hold on 

plot(y1, 'r')
plot(y2, 'b')
plot(x_list(1,:),'-k')
plot(x_list(2,:),'-')

legend('Theta_1','Theta_2','Theta_{1ref}','Theta_{2ref}')


hold off 

%%
p = struct ;

p.m1 = .5;
p.m2 = .25; 
p.L1 = .25;
p.L2 = .25; 
p.g = 9.81;




for i=2:1:length(y1)
    
    
    theta = [y1(i); y2(i)]; 
    theta_dot = x_list(3:4,i); 
    pos = forward_kinematics(p,theta, theta_dot);
    
    x_vec = pos.x;
    y_vec = pos.y;
    
    plot_robot(x_vec, y_vec);
    pause(.1)
end 




%% MIMO Simulation 



dt = .1;
t = [0:dt:(7-dt)]';

u1 = x_list(1,:);
u2 = x_list(2,:);

sys_1 = [G_cl_x1_u1 G_cl_x2_u2];
inputs = [u1', u2'];

figure() 
lsim(sys_1,inputs, t)


%%


x_targ = [pi; deg2rad(-45);0;0];

x_init = [0;-pi/2; 0; 0] ;

x_new = [];
u_new = []; 
dt = .1;

for k = 1:1:100
    u = -K*(x_init-x_targ);
    x_dot = A_mat*x_init + B_mat*(u);
    
    x_new(:,k) = x_init + dt.*x_dot;
    u_new(:,k) = u; 
    x_init = x_new(:,k);
end 

figure()
subplot(4,1,1)
plot(x_new(1,:))
ylabel('Q1')
title('Full State Feedback Stabilizing Controller X_0=[0,-.5*pi,0,0], X_{ref}=[pi, -.25*pi,0,0]')

subplot(4,1,2)
plot(x_new(2,:))
ylabel('Q2')

subplot(4,1,3)
plot(x_new(3,:))
ylabel('Q1_{dot}')

subplot(4,1,4)
plot(x_new(4,:))
ylabel('Q2_{dot}')
xlabel('Time')



figure()
subplot(2,1,1)
plot(u_new(1,:))
ylabel('Tau_1')
title('Full State Feedback Stabilizing Controller Torque Inputs')

subplot(2,1,2)
plot(u_new(2,:))
ylabel('Tau_2')
 


  
  
  
  %% Full-State Feedback Regulator

x0 = [0,0,0,0];         % Equilibrium Point
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
states_sym = [q1, q2, q1_dot, q2_dot, 9.82]; % State Symbols

params_val = [.25,.25, .5, .25]; 
% states_val = []


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);


G_ss = ss(A_mat, B_mat, C, zeros(4,2));


K = place(A_mat, B_mat, [-5, -4.5, -1, -1]);

% X_stable = (A_mat-B_mat*K);
% 
% eig(X_stable)

r = [-.5;0];




x_init = [pi;-pi/4; 0; 0] ;
X_stable = ss((A_mat-B_mat*K) + B_mat*r, eye(4,2) , C, zeros(4,2));

figure(5)
initial(X_stable, x_init)
title('Full State Feedback Stabilizing Controller X_0=[pi,-.25*pi ,0,0], X_ref=[0,0,0,0]')




%% Full-State Feedback Steady State Reference 
x0 = [0,0,0,0];         % Equilibrium Point
states_sym = [q1, q2, q1_dot, q2_dot]; % State Symbols
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
params_val = [.25,.25, .5, .25, 9.84];


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);


G_ss = ss(A_mat, B_mat, C, zeros(4,2));


% K = place(A_mat, B_mat, [-5, -4.5, -1, -1]);
K = place(A_mat, B_mat, [-7, -7, -2, -2]);


x_targ = [pi; deg2rad(-45);0;0];

x_init = [0;-pi/2; 0; 0] ;

x_new = [];
u_new = []; 
dt = .1;

for k = 1:1:100
    u = -K*(x_init-x_targ);
    x_dot = A_mat*x_init + B_mat*(u);
    
    x_new(:,k) = x_init + dt.*x_dot;
    u_new(:,k) = u; 
    x_init = x_new(:,k);
end 

figure()
subplot(4,1,1)
plot(x_new(1,:))
ylabel('Q1')
title('Full State Feedback Stabilizing Controller X_0=[0,-.5*pi,0,0], X_{ref}=[pi, -.25*pi,0,0]')

subplot(4,1,2)
plot(x_new(2,:))
ylabel('Q2')

subplot(4,1,3)
plot(x_new(3,:))
ylabel('Q1_{dot}')

subplot(4,1,4)
plot(x_new(4,:))
ylabel('Q2_{dot}')
xlabel('Time')



figure()
subplot(2,1,1)
plot(u_new(1,:))
ylabel('Tau_1')
title('Full State Feedback Stabilizing Controller Torque Inputs')

subplot(2,1,2)
plot(u_new(2,:))
ylabel('Tau_2')
 



%% Linear Quadratic Regulator

x0 = [0,0,0,0];         % Equilibrium Point
states_sym = [q1, q2, q1_dot, q2_dot]; % State Symbols
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
params_val = [.25,.25, .5, .25, 9.84];


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);


G_ss = ss(A_mat, B_mat, C, zeros(4,2));


Q = .5.*eye(4); 
R = 1.*eye(2); 

[K,S,CLP] = lqr(G_ss,Q,R); 

x_targ = [0; 0;0;0];

x_init = [pi;-pi/2; 0; 0] ;

x_new = [];
u_new = []; 
dt = .1;

for k = 1:1:100
    u = -K*(x_init-x_targ);
    x_dot = A_mat*x_init + B_mat*(u);
    
    x_new(:,k) = x_init + dt.*x_dot;
    u_new(:,k) = u; 
    x_init = x_new(:,k);
end 

figure()
subplot(4,1,1)
plot(x_new(1,:))
ylabel('Q1')
title('LQR X_0=[0,-.5*pi,0,0]')

subplot(4,1,2)
plot(x_new(2,:))
ylabel('Q2')

subplot(4,1,3)
plot(x_new(3,:))
ylabel('Q1_{dot}')

subplot(4,1,4)
plot(x_new(4,:))
ylabel('Q2_{dot}')
xlabel('Time')

figure()
subplot(2,1,1)
plot(u_new(1,:))
ylabel('Tau_1')
title('LQR Torque Inputs')

subplot(2,1,2)
plot(u_new(2,:))
ylabel('Tau_2')
 
 


%% LQR State Following 

x0 = [0,0,0,0];         % Equilibrium Point
states_sym = [q1, q2, q1_dot, q2_dot]; % State Symbols
params_sym = [L1, L2, m1, m2, g_val];  % Parameter Symbols
params_val = [.25,.25, .5, .25, 9.84];


A_mat = vpa(subs(A, [states_sym, params_sym],[x0, params_val]),4);
B_mat = vpa(subs(B, [states_sym, params_sym],[x0, params_val]),4);
A_mat = double(A_mat);
B_mat = double(B_mat);


G_ss = ss(A_mat, B_mat, C, zeros(4,2));


Q = .5.*eye(4); 
R = 1.*eye(2); 

[K,S,CLP] = lqr(G_ss,Q,R); 

x_targ = [pi; deg2rad(-45);0;0];

x_init = [0;-pi/2; 0; 0] ;

x_new = [];
u_new = []; 
dt = .1;

for k = 1:1:100
    u = -K*(x_init-x_targ);
    x_dot = A_mat*x_init + B_mat*(u);
    
    x_new(:,k) = x_init + dt.*x_dot;
    u_new(:,k) = u; 
    x_init = x_new(:,k);
end 

figure()
subplot(4,1,1)
plot(x_new(1,:))
ylabel('Q1')
title('LQR X_0=[0,-.5*pi,0,0], X_{ref}=[pi, -.25*pi,0,0]')

subplot(4,1,2)
plot(x_new(2,:))
ylabel('Q2')

subplot(4,1,3)
plot(x_new(3,:))
ylabel('Q1_{dot}')

subplot(4,1,4)
plot(x_new(4,:))
ylabel('Q2_{dot}')
xlabel('Time')



figure()
subplot(2,1,1)
plot(u_new(1,:))
ylabel('Tau_1')
title('LQR Torque Inputs')

subplot(2,1,2)
plot(u_new(2,:))
ylabel('Tau_2')
 

%% Functions 


function plot_robot(x_vec, y_vec)

x1 = x_vec(1); 
x2 = x_vec(2); 
y1 = y_vec(1); 
y2 = y_vec(2); 



hold on
fill([-1-0.2*1 1+0.2*1 1+0.2*1 -1-0.2*1], [-1-0.2*1 -1-0.2*1 1+0.2*1 1+0.2*1], 'w'); % Clears Background

plot([ 0, x1], [0, y1], "b");
plot([ x1, x2], [y1, y2], "r");
hold off

title('Planar 2D Robot in Workspace')
xlim([-.6,.6])
ylim([-.6,.6])

xlabel("X-Axis")
ylabel("Y-Axis") 



end

function pos = forward_kinematics(p,theta, theta_dot)
 
L1 = p.L1; 
L2 = p.L2; 
g = p.g; 

q1 = theta(1);
q2 = theta(2);

q1_dot = theta_dot(1); 
q2_dot = theta_dot(2); 


x = [L1*cos(q1) ;L1*cos(q1) + L2*cos(q1 + q2)];
y = [L1*sin(q1) ;L1*sin(q1) + L2*sin(q1 + q2)];

pos = struct ;

pos.x = x; 
pos.y = y; 

end 

