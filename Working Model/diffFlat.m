function [state,input] = diffFlat(initialState,finalState,timeList,param)
%NATHAN WOLF-SONKIN & JIAH JIN DIFF FLATTNESS CALCULATION

theta1_0 = initialState(1);
theta2_0 = initialState(2);
theta3_0 = initialState(3);
dtheta1_0 = initialState(4);
dtheta2_0 = initialState(5);
dtheta3_0 = initialState(6);

theta1_f = finalState(1);
theta2_f = finalState(2);
theta3_f = finalState(3);
dtheta1_f = finalState(4);
dtheta2_f = finalState(5);
dtheta3_f = finalState(6);

L1 = param.L1;
L2 = param.L2;
L3 = param.L3;
m1 = param.m1;
mA = param.mA;
m2 = param.m2;
mB = param.mB;
m3 = param.m3;
b = param.b;
g  = param.g;

finalTime = timeList(length(timeList));

t1Traj = trajGen(theta1_0,theta1_f,dtheta1_0,dtheta1_f,finalTime);
dt1Traj = polyder(t1Traj)';
ddt1Traj = polyder(dt1Traj)';

t2Traj = trajGen(theta2_0,theta2_f,dtheta2_0,dtheta2_f,finalTime);
dt2Traj = polyder(t2Traj)';
ddt2Traj = polyder(dt2Traj)';

t3Traj = trajGen(theta3_0,theta3_f,dtheta3_0,dtheta3_f,finalTime);
dt3Traj = polyder(t3Traj)';
ddt3Traj = polyder(dt3Traj)';

% tau1 = @(t) (L1^2*m1*polyval(ddt1,t))/4 + L1^2*m2*polyval(ddt1,t) + L1^2*m3*polyval(ddt1,t) + (L2^2*m2*polyval(ddt1,t))/4 + (L2^2*m2*polyval(ddt2,t))/4 + L2^2*m3*polyval(ddt1,t) + L2^2*m3*polyval(ddt2,t) + (L3^2*m3*polyval(ddt1,t))/4 + (L3^2*m3*polyval(ddt2,t))/4 + (L3^2*m3*polyval(ddt3,t))/4 + L1^2*mA*polyval(ddt1,t) + L1^2*mB*polyval(ddt1,t) + L2^2*mB*polyval(ddt1,t) + L2^2*mB*polyval(ddt2,t) + (L2*g*m2*cos(polyval(t1,t) + polyval(t2,t)))/2 + L2*g*m3*cos(polyval(t1,t) + polyval(t2,t)) + L2*g*mB*cos(polyval(t1,t) + polyval(t2,t)) + (L1*g*m1*cos(polyval(t1,t)))/2 + L1*g*m2*cos(polyval(t1,t)) + L1*g*m3*cos(polyval(t1,t)) + L1*g*mA*cos(polyval(t1,t)) + L1*g*mB*cos(polyval(t1,t)) + (L3*g*m3*cos(polyval(t1,t) + polyval(t2,t) + polyval(t3,t)))/2 + L1*L3*m3*cos(polyval(t2,t) + polyval(t3,t))*polyval(ddt1,t) + (L1*L3*m3*cos(polyval(t2,t) + polyval(t3,t))*polyval(ddt2,t))/2 + (L1*L3*m3*cos(polyval(t2,t) + polyval(t3,t))*polyval(ddt3,t))/2 - (L1*L2*m2*sin(polyval(t2,t))*polyval(dt2,t)^2)/2 - L1*L2*m3*sin(polyval(t2,t))*polyval(dt2,t)^2 - (L2*L3*m3*sin(polyval(t3,t))*polyval(dt3,t)^2)/2 - L1*L2*mB*sin(polyval(t2,t))*polyval(dt2,t)^2 + L1*L2*m2*cos(polyval(t2,t))*polyval(ddt1,t) + (L1*L2*m2*cos(polyval(t2,t))*polyval(ddt2,t))/2 + 2*L1*L2*m3*cos(polyval(t2,t))*polyval(ddt1,t) + L1*L2*m3*cos(polyval(t2,t))*polyval(ddt2,t) + L2*L3*m3*cos(polyval(t3,t))*polyval(ddt1,t) + L2*L3*m3*cos(polyval(t3,t))*polyval(ddt2,t) + (L2*L3*m3*cos(polyval(t3,t))*polyval(ddt3,t))/2 + 2*L1*L2*mB*cos(polyval(t2,t))*polyval(ddt1,t) + L1*L2*mB*cos(polyval(t2,t))*polyval(ddt2,t) - L1*L2*m2*sin(polyval(t2,t))*polyval(dt1,t)*polyval(dt2,t) - 2*L1*L2*m3*sin(polyval(t2,t))*polyval(dt1,t)*polyval(dt2,t) - L2*L3*m3*sin(polyval(t3,t))*polyval(dt1,t)*polyval(dt3,t) - L2*L3*m3*sin(polyval(t3,t))*polyval(dt2,t)*polyval(dt3,t) - 2*L1*L2*mB*sin(polyval(t2,t))*polyval(dt1,t)*polyval(dt2,t) - L1*L3*m3*sin(polyval(t2,t) + polyval(t3,t))*polyval(dt1,t)*(polyval(dt2,t) + polyval(dt3,t)) - (L1*L3*m3*sin(polyval(t2,t) + polyval(t3,t))*polyval(dt2,t)*(polyval(dt2,t) + polyval(dt3,t)))/2 - (L1*L3*m3*sin(polyval(t2,t) + polyval(t3,t))*polyval(dt3,t)*(polyval(dt2,t) + polyval(dt3,t)))/2;
% tau2 = @(t) (L2^2*m2*polyval(ddt1,t))/4 + (L2^2*m2*polyval(ddt2,t))/4 + L2^2*m3*polyval(ddt1,t) + L2^2*m3*polyval(ddt2,t) + (L3^2*m3*polyval(ddt1,t))/4 + (L3^2*m3*polyval(ddt2,t))/4 + (L3^2*m3*polyval(ddt3,t))/4 + L2^2*mB*polyval(ddt1,t) + L2^2*mB*polyval(ddt2,t) + (L2*g*m2*cos(polyval(t1,t) + polyval(t2,t)))/2 + L2*g*m3*cos(polyval(t1,t) + polyval(t2,t)) + L2*g*mB*cos(polyval(t1,t) + polyval(t2,t)) + (L3*g*m3*cos(polyval(t1,t) + polyval(t2,t) + polyval(t3,t)))/2 + (L1*L3*m3*sin(polyval(t2,t) + polyval(t3,t))*polyval(dt1,t)^2)/2 + (L1*L3*m3*cos(polyval(t2,t) + polyval(t3,t))*polyval(ddt1,t))/2 + (L1*L2*m2*sin(polyval(t2,t))*polyval(dt1,t)^2)/2 + L1*L2*m3*sin(polyval(t2,t))*polyval(dt1,t)^2 - (L2*L3*m3*sin(polyval(t3,t))*polyval(dt3,t)^2)/2 + L1*L2*mB*sin(polyval(t2,t))*polyval(dt1,t)^2 + (L1*L2*m2*cos(polyval(t2,t))*polyval(ddt1,t))/2 + L1*L2*m3*cos(polyval(t2,t))*polyval(ddt1,t) + L2*L3*m3*cos(polyval(t3,t))*polyval(ddt1,t) + L2*L3*m3*cos(polyval(t3,t))*polyval(ddt2,t) + (L2*L3*m3*cos(polyval(t3,t))*polyval(ddt3,t))/2 + L1*L2*mB*cos(polyval(t2,t))*polyval(ddt1,t) - L2*L3*m3*sin(polyval(t3,t))*polyval(dt1,t)*polyval(dt3,t) - L2*L3*m3*sin(polyval(t3,t))*polyval(dt2,t)*polyval(dt3,t);
% tau3 = @(t) (L3^2*m3*polyval(ddt1,t))/4 + (L3^2*m3*polyval(ddt2,t))/4 + (L3^2*m3*polyval(ddt3,t))/4 + (L3*g*m3*cos(polyval(t1,t) + polyval(t2,t) + polyval(t3,t)))/2 + (L1*L3*m3*sin(polyval(t2,t) + polyval(t3,t))*polyval(dt1,t)^2)/2 + (L1*L3*m3*cos(polyval(t2,t) + polyval(t3,t))*polyval(ddt1,t))/2 + (L2*L3*m3*sin(polyval(t3,t))*polyval(dt1,t)^2)/2 + (L2*L3*m3*sin(polyval(t3,t))*polyval(dt2,t)^2)/2 + (L2*L3*m3*cos(polyval(t3,t))*polyval(ddt1,t))/2 + (L2*L3*m3*cos(polyval(t3,t))*polyval(ddt2,t))/2 + L2*L3*m3*sin(polyval(t3,t))*polyval(dt1,t)*polyval(dt2,t);
% 
% theta1d = @(t) polyval(t1,t);
% theta2d = @(t) polyval(t2,t);
% theta3d = @(t) polyval(t3,t);

state = cell(length(timeList),1);
input = cell(length(timeList),1);
index = 0;
for i = timeList
    index = index + 1;
    t1 = polyval(t1Traj,i);
    t2 = polyval(t2Traj,i);
    t3 = polyval(t3Traj,i);
    dt1 = polyval(dt1Traj,i);
    dt2 = polyval(dt2Traj,i);
    dt3 = polyval(dt3Traj,i);
    ddt1 = polyval(ddt1Traj,i);
    ddt2 = polyval(ddt2Traj,i);
    ddt3 = polyval(ddt3Traj,i);

    tau = [b*dt1 + (L1^2*m1*ddt1)/4 + L1^2*m2*ddt1 + L1^2*m3*ddt1 + (L2^2*m2*ddt1)/4 + (L2^2*m2*ddt2)/4 + L2^2*m3*ddt1 + L2^2*m3*ddt2 + (L3^2*m3*ddt1)/4 + (L3^2*m3*ddt2)/4 + (L3^2*m3*ddt3)/4 + L1^2*mA*ddt1 + L1^2*mB*ddt1 + L2^2*mB*ddt1 + L2^2*mB*ddt2 + (L2*g*m2*cos(t1 + t2))/2 + L2*g*m3*cos(t1 + t2) + L2*g*mB*cos(t1 + t2) + (L1*g*m1*cos(t1))/2 + L1*g*m2*cos(t1) + L1*g*m3*cos(t1) + L1*g*mA*cos(t1) + L1*g*mB*cos(t1) + (L3*g*m3*cos(t1 + t2 + t3))/2 + L1*L3*m3*cos(t2 + t3)*ddt1 + (L1*L3*m3*cos(t2 + t3)*ddt2)/2 + (L1*L3*m3*cos(t2 + t3)*ddt3)/2 - (L1*L2*m2*sin(t2)*dt2^2)/2 - L1*L2*m3*sin(t2)*dt2^2 - (L2*L3*m3*sin(t3)*dt3^2)/2 - L1*L2*mB*sin(t2)*dt2^2 + L1*L2*m2*cos(t2)*ddt1 + (L1*L2*m2*cos(t2)*ddt2)/2 + 2*L1*L2*m3*cos(t2)*ddt1 + L1*L2*m3*cos(t2)*ddt2 + L2*L3*m3*cos(t3)*ddt1 + L2*L3*m3*cos(t3)*ddt2 + (L2*L3*m3*cos(t3)*ddt3)/2 + 2*L1*L2*mB*cos(t2)*ddt1 + L1*L2*mB*cos(t2)*ddt2 - L1*L2*m2*sin(t2)*dt1*dt2 - 2*L1*L2*m3*sin(t2)*dt1*dt2 - L2*L3*m3*sin(t3)*dt1*dt3 - L2*L3*m3*sin(t3)*dt2*dt3 - 2*L1*L2*mB*sin(t2)*dt1*dt2 - L1*L3*m3*sin(t2 + t3)*dt1*(dt2 + dt3) - (L1*L3*m3*sin(t2 + t3)*dt2*(dt2 + dt3))/2 - (L1*L3*m3*sin(t2 + t3)*dt3*(dt2 + dt3))/2;
           b*dt2 + (L2^2*m2*ddt1)/4 + (L2^2*m2*ddt2)/4 + L2^2*m3*ddt1 + L2^2*m3*ddt2 + (L3^2*m3*ddt1)/4 + (L3^2*m3*ddt2)/4 + (L3^2*m3*ddt3)/4 + L2^2*mB*ddt1 + L2^2*mB*ddt2 + (L2*g*m2*cos(t1 + t2))/2 + L2*g*m3*cos(t1 + t2) + L2*g*mB*cos(t1 + t2) + (L3*g*m3*cos(t1 + t2 + t3))/2 + (L1*L3*m3*sin(t2 + t3)*dt1^2)/2 + (L1*L3*m3*cos(t2 + t3)*ddt1)/2 + (L1*L2*m2*sin(t2)*dt1^2)/2 + L1*L2*m3*sin(t2)*dt1^2 - (L2*L3*m3*sin(t3)*dt3^2)/2 + L1*L2*mB*sin(t2)*dt1^2 + (L1*L2*m2*cos(t2)*ddt1)/2 + L1*L2*m3*cos(t2)*ddt1 + L2*L3*m3*cos(t3)*ddt1 + L2*L3*m3*cos(t3)*ddt2 + (L2*L3*m3*cos(t3)*ddt3)/2 + L1*L2*mB*cos(t2)*ddt1 - L2*L3*m3*sin(t3)*dt1*dt3 - L2*L3*m3*sin(t3)*dt2*dt3 + (L1*L3*m3*sin(t2 + t3)*dt1*dt2)/2 + (L1*L3*m3*sin(t2 + t3)*dt1*dt3)/2 - (L1*L3*m3*sin(t2 + t3)*dt1*(dt2 + dt3))/2;
           (m3*((L3^2*ddt1)/2 + (L3^2*ddt2)/2 + (L3^2*ddt3)/2 + L2*L3*cos(t3)*ddt1 + L2*L3*cos(t3)*ddt2 + L1*L3*cos(t2 + t3)*ddt1 - L1*L3*sin(t2 + t3)*dt1*(dt2 + dt3) - L2*L3*sin(t3)*dt1*dt3 - L2*L3*sin(t3)*dt2*dt3))/2 + b*dt3 + (L3*g*m3*cos(t1 + t2 + t3))/2 + (L1*L3*m3*sin(t2 + t3)*dt1^2)/2 + (L2*L3*m3*sin(t3)*dt1^2)/2 + (L2*L3*m3*sin(t3)*dt2^2)/2 + L2*L3*m3*sin(t3)*dt1*dt2 + (L2*L3*m3*sin(t3)*dt1*dt3)/2 + (L2*L3*m3*sin(t3)*dt2*dt3)/2 + (L1*L3*m3*sin(t2 + t3)*dt1*dt2)/2 + (L1*L3*m3*sin(t2 + t3)*dt1*dt3)/2];

    %%%%%%%%%%%%%%%%%%%%%
    currentState = [t1; t2; t3; dt1; dt2; dt3];
    currentInput = [tau(1);tau(2);tau(3)];

    state{index} = currentState;
    input{index} = currentInput;
end

% figure;
% hold on;
% h=gca;h.LineWidth=2;h.FontSize=18;
% plot(time,tau1,'LineWidth',2)
% plot(time,tau2,'LineWidth',2)
% plot(time,tau3,'LineWidth',2)
% xlabel('Time (s)')
% ylabel('Torque (N*m)')
% legend('Joint 1','Joint 2','Joint 3','Location','best')