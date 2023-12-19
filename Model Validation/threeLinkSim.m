clear all; clc;
param.L1 = 1;             % [m]
param.L2 = 1;             % [m]
param.L3 = 1;             % [m]
param.m1 = 1;             % [kg]
param.mA = 1;             % [kg]
param.m2 = 1;             % [kg]
param.mB = 1;             % [kg]
param.m3 = 1;             % [kg]
param.g = 9.81;           % [m/s^2]

%DAMPED SIMULATION SETTINGS
param.b = 3;              % [N*m/(rad/s)]
name = 'Damped';

%UNDAMPED SIMULATION SETTINGS
% param.b = 0;              % [N*m/(rad/s)]
% name = 'Undamped';

t1_0 = 0;           % [rad]
t2_0 = pi/10;           % [rad]
t3_0 = pi/12;           % [rad]
dt1_0 = 0;          % [rad/sec]
dt2_0 = 0;          % [rad/sec]
dt3_0 = 0;          % [rad/sec]
theta_inital = [t1_0;t2_0;t3_0;dt1_0;dt2_0;dt3_0];

timeList = linspace(0,10,1000);

[t,output] = ode45(@(t,output) threeArmSim(t,output,param),timeList,theta_inital);

theta1 = output(:,1);
theta2 = output(:,2);
theta3 = output(:,3);
dtheta1 = output(:,4);
dtheta2 = output(:,5);
dtheta3 = output(:,6);

U = ones(length(theta1),1)*1000;
T = ones(length(theta1),1)*1000;
for i = 1:length(U)
    [Ui,Ti] = totalEnergy(theta1(i),theta2(i),theta3(i),dtheta1(i),dtheta2(i),dtheta3(i),param);
    U(i) = Ui;
    T(i) = Ti;
end
    
%PLOT ENERGY OVER TIME
figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;
plot(t,U,'LineWidth',2)
plot(t,T,'LineWidth',2)
plot(t,U+T,'LineWidth',2)
ylabel('Energy (J)')
xlabel('Time (s)')
% title([name,' Model'])
legend('Potential','Kinetic','Total','Location','best')

fig = figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;h.DataAspectRatio=[1 1 1];
xlabel('X')
ylabel('Y')
title([name,' Model'])
xlim([-3.5 3.5])
ylim([-3.5,3.5])

%SET INITAL POSITIONS
origin = [0,0];
link1end = [param.L1*cos(theta1(1)),param.L1*sin(theta1(1))];
link2end = link1end + [param.L2*cos(theta1(1) + theta2(1)),param.L2*sin(theta1(1) + theta2(1))];
link3end = link2end + [param.L3*cos(theta1(1) + theta2(1) + theta3(1)),param.L3*sin(theta1(1) + theta2(1) + theta3(1))];

%GENERATE INITAL DRAWINGS
link1 = plot([origin(1),link1end(1)],[origin(2),link1end(2)],'r','LineWidth',5);
link2 = plot([link1end(1),link2end(1)],[link1end(2),link2end(2)],'r','LineWidth',5);
link3 = plot([link2end(1),link3end(1)],[link2end(2),link3end(2)],'r','LineWidth',5);
motor1 = plot(origin(1),origin(2),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);
motor2 = plot(link1end(1),link1end(2),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);
motor3 = plot(link2end(1),link2end(2),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);

frameFreq = 2;
for i = 1:frameFreq:length(t)
    %UPDATE TIME
    subtitle(h,['t = ',num2str(t(i))],"FontSize",11)
    %UPDATE POSITIONS
    %SET INITAL POSITIONS
    origin = [0,0];
    link1end = [param.L1*cos(theta1(i)),param.L1*sin(theta1(i))];
    link2end = link1end + [param.L2*cos(theta1(i) + theta2(i)),param.L2*sin(theta1(i) + theta2(i))];
    link3end = link2end + [param.L3*cos(theta1(i) + theta2(i) + theta3(i)),param.L3*sin(theta1(i) + theta2(i) + theta3(i))];
    
    %UPDATE DRAWINGS
    set(link1,'XData',[origin(1),link1end(1)],'YData',[origin(2),link1end(2)])
    set(motor2,'XData',link1end(1),'YData',link1end(2))
    set(link2,'XData',[link1end(1),link2end(1)],'YData',[link1end(2),link2end(2)])
    set(motor3,'XData',link2end(1),'YData',link2end(2))
    set(link3,'XData',[link2end(1),link3end(1)],'YData',[link2end(2),link3end(2)])
    pause(.001)
    drawnow;
    frame = getframe(fig);
    im{i} = frame2im(frame);
end

nImages = length(t);
tDelay = mean(diff(t)) / 3;
filename = [name,'ModelAnimation.gif']; % Specify the output file name
for idx = 1:frameFreq:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",tDelay);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",tDelay);
    end
end

function outputVect = threeArmSim(~,inputVect,param)
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

    t1 = inputVect(1);
    t2 = inputVect(2);
    t3 = inputVect(3);
    dt1 = inputVect(4);
    dt2 = inputVect(5);
    dt3 = inputVect(6);

    E = [(L1*L2*dt2^2*m2*sin(t2))/2 - (L2*g*m2*cos(t1 + t2))/2 - L2*g*m3*cos(t1 + t2) - L2*g*mB*cos(t1 + t2) - (L1*g*m1*cos(t1))/2 - L1*g*m2*cos(t1) - L1*g*m3*cos(t1) - L1*g*mA*cos(t1) - L1*g*mB*cos(t1) - (L3*g*m3*cos(t1 + t2 + t3))/2 - b*dt1 + L1*L2*dt2^2*m3*sin(t2) + (L2*L3*dt3^2*m3*sin(t3))/2 + L1*L2*dt2^2*mB*sin(t2) + (L1*L3*dt2^2*m3*sin(t2 + t3))/2 + (L1*L3*dt3^2*m3*sin(t2 + t3))/2 + L1*L2*dt1*dt2*m2*sin(t2) + 2*L1*L2*dt1*dt2*m3*sin(t2) + L2*L3*dt1*dt3*m3*sin(t3) + L2*L3*dt2*dt3*m3*sin(t3) + 2*L1*L2*dt1*dt2*mB*sin(t2) + L1*L3*dt1*dt2*m3*sin(t2 + t3) + L1*L3*dt1*dt3*m3*sin(t2 + t3) + L1*L3*dt2*dt3*m3*sin(t2 + t3);
                                                                                                                                                                                                                                                                                                                                                                     (L2*L3*dt3^2*m3*sin(t3))/2 - (L2*g*m2*cos(t1 + t2))/2 - L2*g*m3*cos(t1 + t2) - L2*g*mB*cos(t1 + t2) - (L3*g*m3*cos(t1 + t2 + t3))/2 - (L1*L2*dt1^2*m2*sin(t2))/2 - L1*L2*dt1^2*m3*sin(t2) - b*dt2 - L1*L2*dt1^2*mB*sin(t2) - (L1*L3*dt1^2*m3*sin(t2 + t3))/2 + L2*L3*dt1*dt3*m3*sin(t3) + L2*L3*dt2*dt3*m3*sin(t3);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    - b*dt3 - (L3*g*m3*cos(t1 + t2 + t3))/2 - (L2*L3*dt1^2*m3*sin(t3))/2 - (L2*L3*dt2^2*m3*sin(t3))/2 - (L1*L3*dt1^2*m3*sin(t2 + t3))/2 - L2*L3*dt1*dt2*m3*sin(t3)];

    massMatrix = [(L1^2*m1)/4 + L1^2*m2 + L1^2*m3 + (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L1^2*mA + L1^2*mB + L2^2*mB + L1*L2*m2*cos(t2) + 2*L1*L2*m3*cos(t2) + L2*L3*m3*cos(t3) + 2*L1*L2*mB*cos(t2) + L1*L3*m3*cos(t2 + t3), (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L2^2*mB + (L1*L2*m2*cos(t2))/2 + L1*L2*m3*cos(t2) + L2*L3*m3*cos(t3) + L1*L2*mB*cos(t2) + (L1*L3*m3*cos(t2 + t3))/2, (L3^2*m3)/4 + (L2*L3*m3*cos(t3))/2 + (L1*L3*m3*cos(t2 + t3))/2;
                                                                   (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L2^2*mB + (L1*L2*m2*cos(t2))/2 + L1*L2*m3*cos(t2) + L2*L3*m3*cos(t3) + L1*L2*mB*cos(t2) + (L1*L3*m3*cos(t2 + t3))/2,                                                                                                         (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L2^2*mB + L2*L3*m3*cos(t3),                                   (m3*L3^2)/4 + (L2*m3*cos(t3)*L3)/2;
                                                                                                                                                                                 (m3*(L3^2/2 + L2*L3*cos(t3) + L1*L3*cos(t2 + t3)))/2,                                                                                                                                          (m3*(L3^2/2 + L2*cos(t3)*L3))/2,                                                             (L3^2*m3)/4];
    dx1 = inputVect(4:6);
    dx2 = inv(massMatrix)*E;
    outputVect = [dx1;dx2];
end