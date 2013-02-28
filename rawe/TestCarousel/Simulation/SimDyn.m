clear all
close all
clc

Ncvp = 20;
P.param = 1;
P.w = 0;
P.Tf = 25;

Type = 'A'
switch Type
    case 'A'
        load imuDyn.dat
        load LEDposDyn.dat
        imu = imuDyn;
        LEDpos = LEDposDyn;
        RPM = 64;
        P.tu = [imu(:,1) zeros(length(imu),2) 4e-2*imu(:,8) zeros(length(imu),1)];
        P.Tf = 50;
    case 'E1'
        load imu_elevator.dat
        load LEDpos_elevator.dat
        imu = imu_elevator;
        LEDpos = LEDpos_elevator;
        RPM = 62;
        P.tu = [imu(:,1) zeros(length(imu),3) imu(:,10)];
    case 'E2'
        load imu_elev_3e3_5e3.dat
        load LEDpos_elev_3e3_5e3.dat
        imu = imu_elev_3e3_5e3;
        LEDpos = LEDpos_elev_3e3_5e3;
        RPM = 62;
        P.tu = [imu(:,1) zeros(length(imu),3) imu(:,10)];
        P.Tf = 100;
    case 'N'
        load imu.dat
        load LEDpos.dat
        RPM = 65;
        P.tu = [imu(:,1) zeros(length(imu),4)];
        P.Tf = 25;
end


time = linspace(0,P.Tf,Ncvp)';
% X = [sqrt(L^2-R^2)*ones(1,Ncvp);
%      R*cos(2*pi*time);
%      R*sin(2*pi*time)];
% tx = [time X']

Parameters
delta = pi/2;ddelta = 2*pi*RPM/60;dddelta = 0;
ArmKin
xA = XA(1);yA = XA(2);zA = XA(3);
dxA = dXA(1);dyA = dXA(2);dzA = dXA(3);

TetherAngle = 0;
r = 1.2;dr = 0;ZT = 0;
x = XA(1);
y = XA(2) + cos(TetherAngle)*(r-ZT);
z = XA(3) + sin(TetherAngle)*(r-ZT);
dx = -y*ddelta;dy = 0;dz = 0;

e1 = [dx;dy;dz];e1 = e1/norm(e1);
e3 = [x;y;z] - XA;e3 = e3/norm(e3);
e2 = cross(e3,e1);
R = [e1 e2 e3];
e11 = R(1,1);e12 = R(2,1);e13 = R(3,1);
e21 = R(1,2);e22 = R(2,2);e23 = R(3,2);
e31 = R(1,3);e32 = R(2,3);e33 = R(3,3);

w = R'*[0;0;ddelta];
w1 = w(1);w2 = w(2);w3 = w(3);

% p = [x;y;z];dp = [dx;dy;dz];
% 
% XTE = [0;0;ZT] %Tether attachement point in E (Aircraft reference frame)
% 
% %TETHER CONSTRAINT
% Xc = p + R*XTE;
% dXc = dp + R*cross(w,XTE);
% 
% cvec = Xc-XA;
% dcvec = dXc - dXA;
% 
% c = 0.5*(cvec'*cvec - r^2)
% 
% return

ConstFile
c
dc
X0 = [x y z reshape(R,1,9) dx dy dz w1 w2 w3 r dr delta ddelta]';

%

ax = 2;
%DisplayStatic(X0,tu(1,:),ax,X0')

%SIMULATE


intoptions.AbsTol = 1e-4;
intoptions.RelTol = 1e-4;
intoptions.MaxStep = 1/10;
intoptions.Events = 'off';

tic
%tu = [time zeros(Ncvp,3)];
[time,state]=ode45('Model_integ',[0 P.Tf],X0,intoptions,P);

[time,state]=ode45('Model_integ',[0 P.Tf],state(end,:),intoptions,P);%state(end,:)


ET = toc;
RelTime = num2str(P.Tf/ET);
display(['Simulation time ratio: ',RelTime]);

for k = 1:length(time)
    Out(:,k) = Model_integ(time(k),state(k,:),'output',P);
end

figure(33)
plot(time,state(:,3))
ylabel('z')
drawnow

timeDS = linspace(0,P.Tf,P.Tf*60)';
stateDS = [];
for k = 1:length(timeDS)
    stateDS(k,:) = linint([time state],timeDS(k))';
end

figure(34);clf
timeIMU = linspace(0,length(imu)/800,length(imu));
for k = 1:6
    subplot(2,3,k)
    hold on
    plot(timeIMU,imu(:,k+1),'','linewidth',2)
end

figure(34);
Labels = {'a_x','a_y','a_z','w_1','w_2','w_3'};
for k = 1:3
    subplot(2,3,k+3);hold on
    plot(time,Out(k+1,:),'r','linewidth',2)
    ylabel(Labels{k})
    xlim([0 time(end)])
    grid
end
for k = 4:6
    subplot(2,3,k-3);hold on
    plot(time,Out(k+1,:),'r','linewidth',2)
    ylabel(Labels{k})
     xlim([0 time(end)])
     grid
end

figure(36);clf
for k = 1:length(timeIMU)
    Nw(k) = norm(imu(k,2:4));
    Na(k) = norm(imu(k,5:7));
end
for k = 1:length(time)
    Nwmodel(k) = norm(Out(5:7,k));
    Namodel(k) = norm(Out(2:4,k));
end
subplot(2,1,1)
plot(timeIMU,Nw,'b');hold on
plot(time,Nwmodel,'r')
subplot(2,1,2)
plot(timeIMU,Na,'b');hold on
plot(time,Namodel,'r')

figure(38)
plot(time,Out(1,:))

display(['Averaged measure w :',num2str(60*mean(Nw)/2/pi),' RPM'])
display(['Averaged model w :',num2str(60*mean(Nwmodel)/2/pi),' RPM'])

display(['Averaged measure a :',num2str(mean(Na))])
display(['Averaged model a :',num2str(mean(Namodel))])
% 
figure(37);clf;clc
subplot(2,1,1)
scatter(LEDpos(:,4),-LEDpos(:,3),'or');hold on
scatter(LEDpos(:,6),-LEDpos(:,5),'og');hold on
scatter(LEDpos(:,8),-LEDpos(:,7),'ob');hold on

scatter(Out(9,:),-Out(8,:),'+g');hold on
scatter(Out(11,:),-Out(10,:),'+b');hold on
scatter(Out(13,:),-Out(12,:),'+r');hold on
axis equal

subplot(2,1,2)
scatter(LEDpos(:,10),-LEDpos(:,9),'or');hold on
scatter(LEDpos(:,12),-LEDpos(:,11),'og');hold on
scatter(LEDpos(:,14),-LEDpos(:,13),'ob');hold on

scatter(Out(15,:),-Out(14,:),'+g');hold on
scatter(Out(17,:),-Out(16,:),'+b');hold on
scatter(Out(19,:),-Out(18,:),'+r');hold on
axis equal


figure(39);clf
LED = [5 6 7 8 3 4 11 12 13 14 9 10];
for k = 1:6
    subplot(3,2,k)
    plot(LEDpos(:,1),LEDpos(:,LED(k)),'b','linewidth',2);hold on
    plot(time,Out(7+k,:),'r','linewidth',2)
    ylim([-1600 1600])
    grid
end

figure(40);clf
LED = [5 6 7 8 3 4 11 12 13 14 9 10];
for k = 7:12
    subplot(3,2,k-6)
    plot(LEDpos(:,1),LEDpos(:,LED(k)),'b','linewidth',2);hold on
    plot(time,Out(7+k,:),'r','linewidth',2)
    ylim([-1600 1600])
    grid
end

figure(41);clf
plot(time,Out(end,:),'b','linewidth',2)
