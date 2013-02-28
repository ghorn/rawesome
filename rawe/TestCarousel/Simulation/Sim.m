clear all
close all
clc

P.Tf = 10;
Ncvp = 20;
P.param = 1;
P.w = 0;
RPM = 64;

time = linspace(0,P.Tf,Ncvp)';
tu = [time zeros(Ncvp,3)];
% X = [sqrt(L^2-R^2)*ones(1,Ncvp);
%      R*cos(2*pi*time);
%      R*sin(2*pi*time)];
% tx = [time X']

Parameters
delta = 0;ddelta = 2*pi*RPM/60;dddelta = 0;
ArmKin
xA = XA(1);yA = XA(2);zA = XA(3);
dxA = dXA(1);dyA = dXA(2);dzA = dXA(3);

TetherAngle = 0;
r = 1.2;dr = 0;ZT = 0;
x = r;
y = 0;
z = 0;
dx = 0;dy = 0;dz = 0;

e1 = [0;1;0];e2 = [0;0;1];e3 = [1;0;0];
R = [e1 e2 e3];
e11 = R(1,1);e12 = R(2,1);e13 = R(3,1);
e21 = R(1,2);e22 = R(2,2);e23 = R(3,2);
e31 = R(1,3);e32 = R(2,3);e33 = R(3,3);

w = R'*[0;0;ddelta];
w1 = w(1);w2 = w(2);w3 = w(3);


KinFile
c
dc
X0 = [x y z reshape(R,1,9) dx dy dz w1 w2 w3 delta ddelta]'
return
% eq0 = X0(1:end-2);
% X1 = X0;
% %
% eq = fsolve(@equilibrium,eq0)
% 
% X0 = [eq; delta; ddelta];
% X0(3) = 0;
% X0(6) = 0;

% load Xf
% X0 = Xf;

ax = 2;
%DisplayStatic(X0,tu(1,:),ax,X0')

%SIMULATE


intoptions.AbsTol = 1e-4;
intoptions.RelTol = 1e-4;
intoptions.MaxStep = 1/10;
intoptions.Events = 'off';

P.tu =  [0 0 0 0 0;
         P.Tf 0 0 0 0];
tic
%for k = 1:Ncvp
P.r = r;
[time,state]=ode45('Model_integ',[0 P.Tf],X0,intoptions,P);

%[time,state]=ode45('Model_integ',[0 P.Tf],state(end,:),intoptions,P);

%end
ET = toc;
RelTime = num2str(P.Tf/ET);
display(['Simulation time ratio: ',RelTime]);

% for k = 1:length(time)
%     Out(:,k) = Model_integ(time(k),state(k,:),'output',P);
% 
% end

figure(1)
Lab = {'x','y','z'};
for k = 1:3
    subplot(2,2,k)
    plot(time,state(:,k))
    ylabel(Lab{k})
end
drawnow


figure(2)
for k = 4:12
    subplot(3,3,k-3)
    plot(time,state(:,k))
end

figure(3)
for k = 16:18
    subplot(3,1,k-15)
    plot(time,state(:,k))
end
drawnow