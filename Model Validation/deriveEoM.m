clear all; clc;

%DEFINE VARIABLES
syms L1 L2 L3
syms t1(t) t2(t) t3(t)  % theta1 theta2 theta2
syms m1 mA m2 mB m3
syms g                  % gravity time
syms b
syms tau1(t) tau2(t) tau3(t)

%GENERATE POSITION VECTORS
P_1 = .5*L1*[cos(t1);sin(t1)];
P_A = L1*[cos(t1);sin(t1)];
P_2 = P_A + .5*L2*[cos(t1+t2);sin(t1+t2)];
P_B = P_A + L2*[cos(t1+t2);sin(t1+t2)];
P_3 = P_B + .5*L3*[cos(t1+t2+t3);sin(t1+t2+t3)];

%GENERATE VELOCITY VECTORS
dP_1 = diff(P_1,t);
dP_A = diff(P_A,t);
dP_2 = diff(P_2,t);
dP_B = diff(P_B,t);
dP_3 = diff(P_3,t);

%CALCULATE POTENTIAL ENERGY
y1 = P_1(t);
yA = P_A(t);
y2 = P_2(t);
yB = P_B(t);
y3 = P_3(t);

U1 = m1*g*y1(2);
UA = mA*g*yA(2);
U2 = m2*g*y2(2);
UB = mB*g*yB(2);
U3 = m3*g*y3(2);
U = U1 + UA + U2 + UB + U3;

%CACULATE KINETIC ENERGY
T1 = .5*m1*(dP_1.'*dP_1);
TA = .5*mA*(dP_A.'*dP_A);
T2 = .5*m2*(dP_2.'*dP_2);
TB = .5*mB*(dP_B.'*dP_B);
T3 = .5*m3*(dP_3.'*dP_3);
T = T1 + TA + T2 + TB + T3;

%DEFINE LAGRANGIAN
L = T - U;

%CALCULATE (d/dt)(dL/d(dt))
L_dt1 = simplify(diff(L,diff(t1(t),t)));
ddt_L_dt1 = diff(L_dt1,t);

L_dt2 = simplify(diff(L,diff(t2(t),t)));
ddt_L_dt2 = diff(L_dt2,t);

L_dt3 = simplify(diff(L,diff(t3(t),t)));
ddt_L_dt3 = diff(L_dt3,t);

%CALCULATE dL/Dt
L_t1 = simplify(diff(L,t1(t)));
L_t2 = simplify(diff(L,t2(t)));
L_t3 = simplify(diff(L,t3(t)));

alleqn(1) = ddt_L_dt1 - L_t1 + b*diff(t1, t);
alleqn(2) = ddt_L_dt2 - L_t2 + b*diff(t2, t);
alleqn(3) = ddt_L_dt3 - L_t3 + b*diff(t3, t);
alleqn = alleqn.';

%ISOLATE MASS MATRIX TERMS
MMrow1 = collect(collect(collect(ddt_L_dt1,diff(t1(t), t, t)),diff(t2(t), t, t)),diff(t3(t), t, t));
massMatrix(1,1) = diff(MMrow1,diff(t1(t), t, t));
massMatrix(1,2) = diff(MMrow1,diff(t2(t), t, t));
massMatrix(1,3) = diff(MMrow1,diff(t3(t), t, t));

MMrow2 = collect(collect(collect(ddt_L_dt2,diff(t1(t), t, t)),diff(t2(t), t, t)),diff(t3(t), t, t));
massMatrix(2,1) = diff(MMrow2,diff(t1(t), t, t));
massMatrix(2,2) = diff(MMrow2,diff(t2(t), t, t));
massMatrix(2,3) = diff(MMrow2,diff(t3(t), t, t));

MMrow3 = collect(collect(collect(ddt_L_dt3,diff(t1(t), t, t)),diff(t2(t), t, t)),diff(t3(t), t, t));
massMatrix(3,1) = diff(MMrow3,diff(t1(t), t, t));
massMatrix(3,2) = diff(MMrow3,diff(t2(t), t, t));
massMatrix(3,3) = diff(MMrow3,diff(t3(t), t, t));

ddtheta = [diff(t1(t), t, t);diff(t2(t), t, t);diff(t3(t), t, t)];

E = simplify(-alleqn + massMatrix*ddtheta);

input = [tau1;tau2;tau3];
state = [t1;t2;t3;diff(t1, t);diff(t2, t);diff(t3, t)];
F = [diff(t1, t); diff(t2, t); diff(t3, t); inv(massMatrix)*(E + input)];
F = simplify(F);

A = jacobian(F,state);
A = simplify(A);
B = jacobian(F,input);
B = simplify(B);

%SIMPLIFYING EXPRESSION FOR COPYPASTE
syms dt1 dt2 dt3
E = subs(E, diff(t1),dt1);
E = subs(E, diff(t2),dt2);
E = subs(E, diff(t3),dt3);
E = simplify(E);