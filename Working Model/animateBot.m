function animateBot(theta1,theta2,theta3,theta1d,theta2d,theta3d,param,t,name,gifName,timeList)

%ANIMATE ROBOT ARM
fig = figure;
hold on;
h=gca;h.LineWidth=2;h.FontSize=18;h.DataAspectRatio=[1 1 1];
xlabel('X')
ylabel('Y')
title(name)
xlim([-3.5 3.5])
ylim([-3.5,3.5])

%SET INITAL POSITIONS
origin = [0,0];
link1end = [param.L1*cos(theta1(1)),param.L1*sin(theta1(1))];
link2end = link1end + [param.L2*cos(theta1(1) + theta2(1)),param.L2*sin(theta1(1) + theta2(1))];
link3end = link2end + [param.L3*cos(theta1(1) + theta2(1) + theta3(1)),param.L3*sin(theta1(1) + theta2(1) + theta3(1))];

link1end_des = [param.L1*cos(theta1d(1)),param.L1*sin(theta1d(1))];
link2end_des = link1end_des + [param.L2*cos(theta1d(1) + theta2d(1)),param.L2*sin(theta1d(1) + theta2d(1))];
link3end_des = link2end_des + [param.L3*cos(theta1d(1) + theta2d(1) + theta3d(1)),param.L3*sin(theta1d(1) + theta2d(1) + theta3d(1))];

%GENERATE INITIAL DESIRED DRAWINGS
link1_des = plot([origin(1),link1end_des(1)],[origin(2),link1end_des(2)],'LineWidth',5,'Color',[1,0,0,.25]);
link2_des = plot([link1end_des(1),link2end_des(1)],[link1end_des(2),link2end_des(2)],'LineWidth',5,'Color',[1,0,0,.25]);
link3_des = plot([link2end_des(1),link3end_des(1)],[link2end_des(2),link3end_des(2)],'LineWidth',5,'Color',[1,0,0,.25]);
motor1_des = plot(origin(1),origin(2),'o','MarkerFaceColor','k','MarkerSize',5);
motor2_des = plot(link1end_des(1),link1end_des(2),'o','MarkerFaceColor','k','MarkerSize',5);
motor3_des = plot(link2end_des(1),link2end_des(2),'o','MarkerFaceColor','k','MarkerSize',5);

%GENERATE INITAL DRAWINGS
link1 = plot([origin(1),link1end(1)],[origin(2),link1end(2)],'r','LineWidth',5);
link2 = plot([link1end(1),link2end(1)],[link1end(2),link2end(2)],'r','LineWidth',5);
link3 = plot([link2end(1),link3end(1)],[link2end(2),link3end(2)],'r','LineWidth',5);
motor1 = plot(origin(1),origin(2),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);
motor2 = plot(link1end(1),link1end(2),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);
motor3 = plot(link2end(1),link2end(2),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);

frameFreq = 2;
for i = 1:frameFreq:length(theta1)
    %UPDATE TIME
    subtitle(h,['t = ',num2str(t(i))],"FontSize",11)

    %UPDATE POSITIONS
    origin = [0,0];
    link1end = [param.L1*cos(theta1(i)),param.L1*sin(theta1(i))];
    link2end = link1end + [param.L2*cos(theta1(i) + theta2(i)),param.L2*sin(theta1(i) + theta2(i))];
    link3end = link2end + [param.L3*cos(theta1(i) + theta2(i) + theta3(i)),param.L3*sin(theta1(i) + theta2(i) + theta3(i))];
    
    link1end_des = [param.L1*cos(theta1d(i)),param.L1*sin(theta1d(i))];
    link2end_des = link1end_des + [param.L2*cos(theta1d(i) + theta2d(i)),param.L2*sin(theta1d(i) + theta2d(i))];
    link3end_des = link2end_des + [param.L3*cos(theta1d(i) + theta2d(i) + theta3d(i)),param.L3*sin(theta1d(i) + theta2d(i) + theta3d(i))];

    %UPDATE DRAWINGS
    set(link1,'XData',[origin(1),link1end(1)],'YData',[origin(2),link1end(2)])
    set(motor2,'XData',link1end(1),'YData',link1end(2))
    set(link2,'XData',[link1end(1),link2end(1)],'YData',[link1end(2),link2end(2)])
    set(motor3,'XData',link2end(1),'YData',link2end(2))
    set(link3,'XData',[link2end(1),link3end(1)],'YData',[link2end(2),link3end(2)])

    %UPDATE DESIRED DRAWINGS
    set(link1_des,'XData',[origin(1),link1end_des(1)],'YData',[origin(2),link1end_des(2)])
    set(motor2_des,'XData',link1end_des(1),'YData',link1end_des(2))
    set(link2_des,'XData',[link1end_des(1),link2end_des(1)],'YData',[link1end_des(2),link2end_des(2)])
    set(motor3_des,'XData',link2end_des(1),'YData',link2end_des(2))
    set(link3_des,'XData',[link2end_des(1),link3end_des(1)],'YData',[link2end_des(2),link3end_des(2)])
    pause(.01)
    drawnow;
    frame = getframe(fig);
    im{i} = frame2im(frame);
end

nImages = length(theta1);
tDelay = mean(diff(timeList));
filename = gifName; % Specify the output file name
for idx = 1:frameFreq:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",tDelay);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",tDelay);
    end
end