clear all; clc;
%SET PARAMETERS
param.L1 = 1;               % [m]
param.L2 = 1;               % [m]
param.L3 = 1;               % [m]
param.m1 = 1;               % [kg]
param.mA = 1;               % [kg]
param.m2 = 1;               % [kg]
param.mB = 1;               % [kg]
param.m3 = 1;               % [kg]
param.b = 50;               % [N*m/(rad/s)]
param.g = 9.81;             % [m/s^2]

%SET SATURATION BOUNDS
param.highSat = 100;        % [N*m]
param.lowSat = -100;        % [N*m]

%SET FINAL TIME
finalTime = 30;     % [s]

%SET BOUNDRY CONDITIONS
theta1_0 = -pi/4;       % [rad]
theta1_f = pi/3;       % [rad]
dtheta1_0 = 0;      % [rad/sec]
dtheta1_f = 0;      % [rad/sec]
theta2_0 = pi/2;       % [rad]
theta2_f = pi/4;       % [rad]
dtheta2_0 = 0;      % [rad/sec]
dtheta2_f = 0;      % [rad/sec]
theta3_0 = -5*pi/6;       % [rad]
theta3_f = 0;       % [rad]
dtheta3_0 = 0;      % [rad/sec]
dtheta3_f = 0;      % [rad/sec]

initialState = [theta1_0;theta2_0;theta3_0;dtheta1_0;dtheta2_0;dtheta3_0];
finalState = [theta1_f;theta2_f;theta3_f;dtheta1_f;dtheta2_f;dtheta3_f];

%CALCULATE TRAJECTORY FOR INPUT TORQUES/DESIRED ANGLES
timeList = linspace(0,finalTime,1000);
[state_des,input_des] = diffFlat(initialState,finalState,timeList,param);

%LINEARIZE SYSTEM ABOUT ERROR FOR ALL TIMES
K = cell(length(timeList),1);
index = 0;
for i = timeList
    index = index + 1;
    [A,B] = linearize(state_des{index},input_des{index},param);
    Q = diag([1 1 1 1 1 1]);
    R = diag([1 1 1]);
    K{index} = lqr(A,B,Q,R);
    %x'Qx + u'Ru
end 

%SET ODE INITIAL CONDITIONS
disturb(1) = -5 * pi/180;
disturb(2) = -5 * pi/180;
disturb(3) = -5 * pi/180;
disturb(4) = 0;
disturb(5) = 0;
disturb(6) = 0;
disturb = disturb';
initialState = [theta1_0;theta2_0;theta3_0;dtheta1_0;dtheta2_0;dtheta3_0] + disturb; %USE LATER FOR PERTURBATIONS

%RUN OPTIMAL CONTROL SIMULATION
[~,optimalSim] = ode45(@(t,output) threeArmSim(t,output,state_des,input_des,K,timeList,param),timeList,initialState);

%RUN FF PROPORTIONAL SIMULATION
Kp = 50;
for i = 1:length(K)
    ff.K{i} = zeros(3,6);
    ff.K{i}(1,1) = Kp;
    ff.K{i}(2,2) = Kp;
    ff.K{i}(3,3) = Kp;
end
[~,forwardSim] = ode45(@(t,output) threeArmSim(t,output,state_des,input_des,ff.K,timeList,param),timeList,initialState);

%RUN PROPORTIONAL ONLY SIMULATION
Kp = 400;
for i = 1:length(K)
    fb.K{i} = zeros(3,6);
    fb.K{i}(1,1) = Kp;
    fb.K{i}(2,2) = Kp;
    fb.K{i}(3,3) = Kp;
    fb.input_des{i} = [0;0;0];
end
[t,backSim] = ode45(@(t,output) threeArmSim(t,output,state_des,fb.input_des,fb.K,timeList,param),timeList,initialState);

%RETRIEVE ERROR DATA FOR EACH SIM
opt.stateError = cell(length(t),1);
opt.inputEffort = cell(length(t),1);
ff.stateError = cell(length(t),1);
ff.inputEffort = cell(length(t),1);
fb.stateError = cell(length(t),1);
fb.inputEffort = cell(length(t),1);
for i = 1:length(t)
    [~,opt.stateError{i},opt.inputEffort{i}] = threeArmSim(t(i),optimalSim(i,:).',state_des,input_des,K,timeList,param);
    [~,ff.stateError{i},ff.inputEffort{i}] = threeArmSim(t(i),forwardSim(i,:).',state_des,input_des,ff.K,timeList,param);
    [~,fb.stateError{i},fb.inputEffort{i}] = threeArmSim(t(i),backSim(i,:).',state_des,fb.input_des,fb.K,timeList,param);
end

%FORMAT DATA FOR ANIMATION
opt.theta1 = optimalSim(:,1);
opt.theta2 = optimalSim(:,2);
opt.theta3 = optimalSim(:,3);
ff.theta1 = forwardSim(:,1);
ff.theta2 = forwardSim(:,2);
ff.theta3 = forwardSim(:,3);
fb.theta1 = backSim(:,1);
fb.theta2 = backSim(:,2);
fb.theta3 = backSim(:,3);

tau1 = zeros(length(timeList),1);
tau2 = zeros(length(timeList),1);
tau3 = zeros(length(timeList),1);
theta1d = zeros(length(timeList),1);
theta2d = zeros(length(timeList),1);
theta3d = zeros(length(timeList),1);
for i = 1:length(timeList)
    tau1(i) = input_des{i}(1,1);
    tau2(i) = input_des{i}(2,1);
    tau3(i) = input_des{i}(3,1);
    theta1d(i) = state_des{i}(1,1);
    theta2d(i) = state_des{i}(2,1);
    theta3d(i) = state_des{i}(3,1);
end

% PLOT FEEDFORWARD TORQUE TRAJECTORY
figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
plot(timeList,tau1,'LineWidth',2)
plot(timeList,tau2,'LineWidth',2)
plot(timeList,tau3,'LineWidth',2)
xlabel('Time (s)')
ylabel('Torque (N*m)')
legend('Joint 1','Joint 2','Joint 3','Location','best')

for i = 1:length(t)
    %OPTIMAL CONTROLLER
    opt.stateError1(i) = opt.stateError{i}(1) * 180/pi;
    opt.stateError2(i) = opt.stateError{i}(2) * 180/pi;
    opt.stateError3(i) = opt.stateError{i}(3) * 180/pi;
    opt.inputEffort1(i) = opt.inputEffort{i}(1);
    opt.inputEffort2(i) = opt.inputEffort{i}(2);
    opt.inputEffort3(i) = opt.inputEffort{i}(3);
    
    %FEEDFORWRD/FEEDBACK CONTROLLER
    ff.stateError1(i) = ff.stateError{i}(1) * 180/pi;
    ff.stateError2(i) = ff.stateError{i}(2) * 180/pi;
    ff.stateError3(i) = ff.stateError{i}(3) * 180/pi;
    ff.inputEffort1(i) = ff.inputEffort{i}(1);
    ff.inputEffort2(i) = ff.inputEffort{i}(2);
    ff.inputEffort3(i) = ff.inputEffort{i}(3);

    %FEEDBACK CONTROLLER
    fb.stateError1(i) = fb.stateError{i}(1) * 180/pi;
    fb.stateError2(i) = fb.stateError{i}(2) * 180/pi;
    fb.stateError3(i) = fb.stateError{i}(3) * 180/pi;
    fb.inputEffort1(i) = fb.inputEffort{i}(1);
    fb.inputEffort2(i) = fb.inputEffort{i}(2);
    fb.inputEffort3(i) = fb.inputEffort{i}(3);
end
%%
%DOUBLE PLOT SHOWING P CONTROLLER PERFORMANCE
figure;
subplot(2,1,1)
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint Error (°)')
plot(t,abs(fb.stateError1),'r','LineWidth',2)
plot(t,abs(fb.stateError2),'b','LineWidth',2)
plot(t,abs(fb.stateError3),'g','LineWidth',2)
plot(t,abs(fb.stateError1)+abs(fb.stateError2)+abs(fb.stateError3),'k','LineWidth',2)
legend('Joint 1','Joint 2','Joint 3','Total','Location','best')

subplot(2,1,2)
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Actuation Effort (N*m)')
plot(t,abs(fb.inputEffort1),'r','LineWidth',2)
plot(t,abs(fb.inputEffort2),'b','LineWidth',2)
plot(t,abs(fb.inputEffort3),'g','LineWidth',2)
plot(t,abs(fb.inputEffort1)+abs(fb.inputEffort2)+abs(fb.inputEffort3),'k','LineWidth',2)
legend('Joint 1','Joint 2','Joint 3','Total','Location','best')
%%
%DOUBLE PLOT SHOWING FF+P CONTROLLER PERFORMANCE
figure;
subplot(2,1,1)
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint Error (°)')
plot(t,abs(ff.stateError1),'r','LineWidth',2)
plot(t,abs(ff.stateError2),'b','LineWidth',2)
plot(t,abs(ff.stateError3),'g','LineWidth',2)
plot(t,abs(ff.stateError1)+abs(ff.stateError2)+abs(ff.stateError3),'k','LineWidth',2)
legend('Joint 1','Joint 2','Joint 3','Total','Location','best')

subplot(2,1,2)
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Actuation Effort (N*m)')
plot(t,abs(ff.inputEffort1),'r','LineWidth',2)
plot(t,abs(ff.inputEffort2),'b','LineWidth',2)
plot(t,abs(ff.inputEffort3),'g','LineWidth',2)
plot(t,abs(ff.inputEffort1)+abs(ff.inputEffort2)+abs(ff.inputEffort3),'k','LineWidth',2)
legend('Joint 1','Joint 2','Joint 3','Total','Location','best')
%%

%DOUBLE PLOT SHOWING FF+OPT CONTROLLER PERFORMANCE
figure;
subplot(2,1,1)
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint Error (°)')
plot(t,abs(opt.stateError1),'r','LineWidth',2)
plot(t,abs(opt.stateError2),'b','LineWidth',2)
plot(t,abs(opt.stateError3),'g','LineWidth',2)
plot(t,abs(opt.stateError1)+abs(opt.stateError2)+abs(opt.stateError3),'k','LineWidth',2)
legend('Joint 1','Joint 2','Joint 3','Total','Location','best')

subplot(2,1,2)
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Actuation Effort (N*m)')
plot(t,abs(opt.inputEffort1),'r','LineWidth',2)
plot(t,abs(opt.inputEffort2),'b','LineWidth',2)
plot(t,abs(opt.inputEffort3),'g','LineWidth',2)
plot(t,abs(opt.inputEffort1)+abs(opt.inputEffort2)+abs(opt.inputEffort3),'k','LineWidth',2)
legend('Joint 1','Joint 2','Joint 3','Total','Location','best')
%%
% PLOT JOINT ERROR COMPARISON
figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint 1 Error (°)')
plot(t,opt.stateError1,'r','LineWidth',2)
plot(t,ff.stateError1,'b','LineWidth',2)
plot(t,fb.stateError1,'g','LineWidth',2)
legend('Feedforward + Optimal','Feedforward + PID','PID Only','Location','best')

figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint 2 Error (°)')
plot(t,opt.stateError2,'r','LineWidth',2)
plot(t,ff.stateError2,'b','LineWidth',2)
plot(t,fb.stateError2,'g','LineWidth',2)
legend('Feedforward + Optimal','Feedforward + PID','PID Only','Location','best')

figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint 3 Error (°)')
plot(t,opt.stateError3,'r','LineWidth',2)
plot(t,ff.stateError3,'b','LineWidth',2)
plot(t,fb.stateError3,'g','LineWidth',2)
legend('Feedforward + Optimal','Feedforward + PID','PID Only','Location','best')

% PLOT ACTUATION EFFORT
% PLOTTING EACH TORQUE EFFORT
figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint 1 Applied Torque (N*m)')
plot(t,opt.inputEffort1,'r','LineWidth',2)
plot(t,ff.inputEffort1,'b','LineWidth',2)
plot(t,fb.inputEffort1,'g','LineWidth',2)
legend('Feedforward + Optimal','Feedforward + PID','PID Only','Location','best')

figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint 2 Applied Torque (N*m)')
plot(t,opt.inputEffort2,'r','LineWidth',2)
plot(t,ff.inputEffort2,'b','LineWidth',2)
plot(t,fb.inputEffort2,'g','LineWidth',2)
legend('Feedforward + Optimal','Feedforward + PID','PID Only','Location','best')

figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Joint 3 Applied Torque (N*m)')
plot(t,opt.inputEffort3,'r','LineWidth',2)
plot(t,ff.inputEffort3,'b','LineWidth',2)
plot(t,fb.inputEffort3,'g','LineWidth',2)
legend('Feedforward + Optimal','Feedforward + PID','PID Only','Location','best')
%%
%
%PLOT TOTAL ERROR
figure;
subplot(2,1,1)
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Error (°)')
% title('Total Error')
plot(t,abs(opt.stateError1) + abs(opt.stateError2) + abs(opt.stateError3),'r','LineWidth',2)
plot(t,abs(ff.stateError1) + abs(ff.stateError2) + abs(ff.stateError3),'b','LineWidth',2)
plot(t,abs(fb.stateError1) + abs(fb.stateError2) + abs(fb.stateError3),'g','LineWidth',2)
legend('Feedforward + Optimal','Feedforward + P','P Only','Location','best')

%PLOT TOTAL TORQUE INPUT
subplot(2,1,2)
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
xlabel('Time (s)')
ylabel('Applied Torque (N*m)')
% title('Total Applied Torque')
plot(t,abs(opt.inputEffort1) + abs(opt.inputEffort2) + abs(opt.inputEffort3),'r','LineWidth',2)
plot(t,abs(ff.inputEffort1) + abs(ff.inputEffort2) + abs(ff.inputEffort3),'b','LineWidth',2)
plot(t,abs(fb.inputEffort1) + abs(fb.inputEffort2) + abs(fb.inputEffort3),'g','LineWidth',2)
legend('Feedforward + Optimal','Feedforward + P','P Only','Location','best')
%%
animateBot(opt.theta1,opt.theta2,opt.theta3,theta1d,theta2d,theta3d,param,t,'Optimal Controller',"OptimalAnimation.gif",timeList);
animateBot(ff.theta1,ff.theta2,ff.theta3,theta1d,theta2d,theta3d,param,t,'Feedforward + P Controller',"FeedforwardAnimation.gif",timeList);
animateBot(fb.theta1,fb.theta2,fb.theta3,theta1d,theta2d,theta3d,param,t,'P Controller',"FeedbackAnimation.gif",timeList);

function [dXdt,stateError,commandedInput] = threeArmSim(t,X,X_des,u_des,K,timeList,param)
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
    highSat = param.highSat;
    lowSat = param.lowSat;

    t1 = X(1);
    t2 = X(2);
    t3 = X(3);
    dt1 = X(4);
    dt2 = X(5);
    dt3 = X(6);

    %USE CONTINUOUS TIME TO FIND CLOSEST DISCRETE TIME INTERVAL
    i = 2;
    while t > timeList(i)
        i = i+1;
    end

    %LINEAR INTERPOLATION OF K
    newK = K{i-1} + ((K{i} - K{i-1})./(timeList(i)-timeList(i-1))) .* (t - timeList(i-1));
    newU = u_des{i-1} + ((u_des{i} - u_des{i-1})./(timeList(i)-timeList(i-1))) .* (t - timeList(i-1));
    newX = X_des{i-1} + ((X_des{i} - X_des{i-1})./(timeList(i)-timeList(i-1))) .* (t - timeList(i-1));

    %CONTROL LAW
    u = newU - 1*newK*(X - newX);
    
    stateError = X - newX;
    commandedInput = u;

    %APPLY SATURATION IF REQUIRED
    for i = 1:length(u)
        if u(i) > highSat
            u(i) = highSat;
        elseif u(i) < lowSat
            u(i) = lowSat;
        end
    end

    %DAMPED SYSTEM
    E = [(L1*L2*dt2^2*m2*sin(t2))/2 - (L2*g*m2*cos(t1 + t2))/2 - L2*g*m3*cos(t1 + t2) - L2*g*mB*cos(t1 + t2) - (L1*g*m1*cos(t1))/2 - L1*g*m2*cos(t1) - L1*g*m3*cos(t1) - L1*g*mA*cos(t1) - L1*g*mB*cos(t1) - (L3*g*m3*cos(t1 + t2 + t3))/2 - b*dt1 + L1*L2*dt2^2*m3*sin(t2) + (L2*L3*dt3^2*m3*sin(t3))/2 + L1*L2*dt2^2*mB*sin(t2) + (L1*L3*dt2^2*m3*sin(t2 + t3))/2 + (L1*L3*dt3^2*m3*sin(t2 + t3))/2 + L1*L2*dt1*dt2*m2*sin(t2) + 2*L1*L2*dt1*dt2*m3*sin(t2) + L2*L3*dt1*dt3*m3*sin(t3) + L2*L3*dt2*dt3*m3*sin(t3) + 2*L1*L2*dt1*dt2*mB*sin(t2) + L1*L3*dt1*dt2*m3*sin(t2 + t3) + L1*L3*dt1*dt3*m3*sin(t2 + t3) + L1*L3*dt2*dt3*m3*sin(t2 + t3);
                                                                                                                                                                                                                                                                                                                                                                     (L2*L3*dt3^2*m3*sin(t3))/2 - (L2*g*m2*cos(t1 + t2))/2 - L2*g*m3*cos(t1 + t2) - L2*g*mB*cos(t1 + t2) - (L3*g*m3*cos(t1 + t2 + t3))/2 - (L1*L2*dt1^2*m2*sin(t2))/2 - L1*L2*dt1^2*m3*sin(t2) - b*dt2 - L1*L2*dt1^2*mB*sin(t2) - (L1*L3*dt1^2*m3*sin(t2 + t3))/2 + L2*L3*dt1*dt3*m3*sin(t3) + L2*L3*dt2*dt3*m3*sin(t3);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    - b*dt3 - (L3*g*m3*cos(t1 + t2 + t3))/2 - (L2*L3*dt1^2*m3*sin(t3))/2 - (L2*L3*dt2^2*m3*sin(t3))/2 - (L1*L3*dt1^2*m3*sin(t2 + t3))/2 - L2*L3*dt1*dt2*m3*sin(t3)];

    massMatrix = [(L1^2*m1)/4 + L1^2*m2 + L1^2*m3 + (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L1^2*mA + L1^2*mB + L2^2*mB + L1*L2*m2*cos(t2) + 2*L1*L2*m3*cos(t2) + L2*L3*m3*cos(t3) + 2*L1*L2*mB*cos(t2) + L1*L3*m3*cos(t2 + t3), (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L2^2*mB + (L1*L2*m2*cos(t2))/2 + L1*L2*m3*cos(t2) + L2*L3*m3*cos(t3) + L1*L2*mB*cos(t2) + (L1*L3*m3*cos(t2 + t3))/2, (L3^2*m3)/4 + (L2*L3*m3*cos(t3))/2 + (L1*L3*m3*cos(t2 + t3))/2;
                                                                   (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L2^2*mB + (L1*L2*m2*cos(t2))/2 + L1*L2*m3*cos(t2) + L2*L3*m3*cos(t3) + L1*L2*mB*cos(t2) + (L1*L3*m3*cos(t2 + t3))/2,                                                                                                         (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L2^2*mB + L2*L3*m3*cos(t3),                                   (m3*L3^2)/4 + (L2*m3*cos(t3)*L3)/2;
                                                                                                                                                                                 (m3*(L3^2/2 + L2*L3*cos(t3) + L1*L3*cos(t2 + t3)))/2,                                                                                                                                          (m3*(L3^2/2 + L2*cos(t3)*L3))/2,                                                             (L3^2*m3)/4];
    dXdt_bot3 = inv(massMatrix)*(E + u);

    dXdt = [X(4);X(5);X(6);dXdt_bot3];
end