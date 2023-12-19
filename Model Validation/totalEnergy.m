function [U,T] = totalEnergy(t1,t2,t3,dt1,dt2,dt3,param)
    L1 = param.L1;
    L2 = param.L2;
    L3 = param.L3;
    m1 = param.m1;
    mA = param.mA;
    m2 = param.m2;
    mB = param.mB;
    m3 = param.m3;
    g  = param.g;

    %CALCULATE ENERGY IN THE SYSTEM TO CONFIRM MODEL ACCURACY
   
    %CALCULATE POSITION VECTORS
    P_1 = .5*L1*[cos(t1);sin(t1)];
    P_A = L1*[cos(t1);sin(t1)];
    P_2 = P_A + .5*L2*[cos(t1+t2);sin(t1+t2)];
    P_B = P_A + L2*[cos(t1+t2);sin(t1+t2)];
    P_3 = P_B + .5*L3*[cos(t1+t2+t3);sin(t1+t2+t3)];
    
    %CALCULATE POTENTIAL ENERGY
    U1 = m1*g*P_1(2);
    UA = mA*g*P_A(2);
    U2 = m2*g*P_2(2);
    UB = mB*g*P_B(2);
    U3 = m3*g*P_3(2);
    U = U1+UA+U2+UB+U3;

    %CALCULATE VELOCITY VECTORS
    dP_1 = .5*L1*dt1*[-sin(t1);cos(t1)];
    dP_A = L1*dt1*[-sin(t1);cos(t1)];
    dP_2 = dP_A + .5*L2*(dt1+dt2)*[-sin(t1+t2);cos(t1+t2)];
    dP_B = dP_A + L2*(dt1+dt2)*[-sin(t1+t2);cos(t1+t2)];
    dP_3 = dP_B + .5*L3*(dt1+dt2+dt3)*[-sin(t1+t2+t3);cos(t1+t2+t3)];

    %CALCULATE TRANSLATIONAL KINETIC ENERGY
    T1 = .5*m1*(dP_1.'*dP_1);
    TA = .5*mA*(dP_A.'*dP_A);
    T2 = .5*m2*(dP_2.'*dP_2);
    TB = .5*mB*(dP_B.'*dP_B);
    T3 = .5*m3*(dP_3.'*dP_3);
    T = T1+TA+T2+TB+T3;
end