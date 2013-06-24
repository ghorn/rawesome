clear all
close all
clc

adat = aerodata();
%a=aerodata();
n = length(adat.Alpha);

% names = [%{'CL_tot','CL'};
%     	 %{'CD_tot','CD'};
% %      	 {'Cl_tot','Cl'};
% %     	 {'Cm_tot','Cm'};
%      	 {'Cn_tot','Cn'};
% %          {'CX_tot','Cx'};
% %          {'CY_tot','Cy'};
% %          {'CZ_tot','Cz'};
%         ];
% 
% [rows,~] = size(names);
% 
% for k=1:rows
% %    subplot(3,2,k);
%     figure()
%     plot3(adat.Alpha,adat.Beta,adat.(char(names(k,1))),'.')
%     xlabel('alpha') adat.Alpha(k) >= -4.5 && adat.Alpha(k) <= 8.5
%     ylabel('beta')
%     zlabel(char(names(k,2)))
% end
% return


% fit linear CL alpha
M = [adat.Alpha*pi/180, ones(n,1)];
res = M\adat.CL_tot;
fprintf('CL = %f*alpha + %f\n', res)

figure()
plot3(adat.Alpha, adat.Beta, adat.CL_tot,'b.')
hold on
plot3(adat.Alpha, adat.Beta, M*res,'r.')
xlabel('alpha')
ylabel('beta')
zlabel('CL')


% fit quadratic CD alpha/beta
alpha = zeros(0,0);
beta  = zeros(0,0);
CD    = zeros(0,0);
for k=1:n
    if adat.Alpha(k) >= -4.5 && adat.Alpha(k) <= 9.5
        alpha = [alpha; adat.Alpha(k)];
        beta = [beta; adat.Beta(k)];
        CD = [CD; adat.CD_tot(k)];
    end
end
M = [(alpha*pi/180).^2, (alpha*pi/180), (beta*pi/180).^2, ones(length(alpha),1)];
res = M\CD;
fprintf('CD = %f*alpha^2 + %f*alpha + %f*beta^2 + %f\n', res)

figure()
plot3(adat.Alpha, adat.Beta, adat.CD_tot,'b.')
hold on
M = [(adat.Alpha*pi/180).^2, (adat.Alpha*pi/180), (adat.Beta*pi/180).^2, ones(n,1)];
plot3(adat.Alpha, adat.Beta, M*res,'r.')
xlabel('alpha')
ylabel('beta')
zlabel('CD')


% fit quadratic CY alpha/beta
alpha = zeros(0,0);
beta  = zeros(0,0);
CY    = zeros(0,0);
for k=1:n
    if adat.Alpha(k) >= -4.5 && adat.Alpha(k) <= 9.5
        alpha = [alpha; adat.Alpha(k)];
        beta = [beta; adat.Beta(k)];
        CY = [CY; adat.CY_tot(k)];
    end
end
M = [(alpha*pi/180).^2, (alpha*pi/180), (beta*pi/180), ones(length(alpha),1)];
res = M\CY;
fprintf('CY = %f*alpha^2 + %f*alpha + %f*beta + %f\n', res)

figure()
plot3(adat.Alpha, adat.Beta, adat.CY_tot,'b.')
hold on
M = [(adat.Alpha*pi/180).^2, (adat.Alpha*pi/180), (adat.Beta*pi/180), ones(n,1)];
plot3(adat.Alpha, adat.Beta, M*res,'r.')
xlabel('alpha')
ylabel('beta')
zlabel('CY')


% fit linear Cm alpha
M = [adat.Alpha*pi/180, ones(n,1)];
res = M\adat.Cm_tot;
fprintf('Cm = %f*alpha + %f\n', res)

figure()
plot3(adat.Alpha, adat.Beta, adat.Cm_tot,'b.')
hold on
plot3(adat.Alpha, adat.Beta, M*res,'r.')
xlabel('alpha')
ylabel('beta')
zlabel('Cm')


% fit Cl alpha beta
M = [adat.Beta*pi/180, adat.Alpha*pi/180.*adat.Beta*pi/180, ones(n,1)];
res = M\adat.Cl_tot;
fprintf('Cl = %f*beta + %f*alpha*beta + %f\n', res)

figure()
plot3(adat.Alpha, adat.Beta, adat.Cl_tot,'b.')
hold on
plot3(adat.Alpha, adat.Beta, M*res,'r.')
xlabel('alpha')
ylabel('beta')
zlabel('Cl')


% fit Cn alpha
alpha = zeros(0,0);
beta  = zeros(0,0);
Cn    = zeros(0,0);
for k=1:n
    if adat.Alpha(k) >= -4.5 && adat.Alpha(k) <= 8.5
        alpha = [alpha; adat.Alpha(k)];
        beta = [beta; adat.Beta(k)];
        Cn = [Cn; adat.Cn_tot(k)];
    end
end
M = [beta*pi/180, alpha*pi/180.*beta*pi/180, ones(length(alpha),1)];
res = M\Cn;
fprintf('Cn = %f*beta + %f*alpha*beta + %f\n', res)

figure()
plot3(adat.Alpha, adat.Beta, adat.Cn_tot,'b.')
hold on
M = [adat.Beta*pi/180, adat.Alpha*pi/180.*adat.Beta*pi/180, ones(n,1)];
plot3(adat.Alpha, adat.Beta, M*res,'r.')
xlabel('alpha')
ylabel('beta')
zlabel('Cn')



% fit linear Cd_dflaps
M = [adat.Flap, adat.Alpha, ones(n,1)];
res = M\adat.CD_ff_d1;
fprintf('CD_dflaps = %g*flaps + %g*alpha + %g\n', res)
fprintf('CD_flaps = %g*flaps^2 + %g*alpha*delta + %g*delta\n',res(1)*0.5, res(2), res(3))

figure()
plot3(adat.Alpha, adat.Flap, adat.CD_ff_d1,'b.')
hold on
plot3(adat.Alpha, adat.Flap, M*res,'r.')
xlabel('alpha')
ylabel('flaps')
zlabel('CD_dflaps')

% integrate previous fit analytically, plot result
alphas = [];
deltas = [];
drags = [];
for alpha=linspace(-4,8,20)
    for delta=linspace(-6,6,20);
        alphas = [alphas, alpha];
        deltas = [deltas, delta];
        drags = [drags,res(1)*0.5*delta.^2 + res(2)*alpha*delta + res(3)*delta];
    end
end
figure()
plot3(alphas,deltas,drags,'b.')
xlabel('alpha');ylabel('delta');zlabel('CD_flaps');
grid on


% fit linear Cd_delev
adat.Elevator = -adat.Elevator;
adat.CD_ff_d3 = -adat.CD_ff_d3;
M = [adat.Elevator, adat.Alpha, ones(n,1)];
res = M\adat.CD_ff_d3;
fprintf('CD_dElevator = %g*Elevator + %g*alpha + %g\n', res)
fprintf('CD_Elevator = %g*Elevator^2 + %g*alpha*delta + %g*delta\n',res(1)*0.5, res(2), res(3))

figure()
plot3(adat.Alpha, adat.Elevator, adat.CD_ff_d3,'b.')
hold on
plot3(adat.Alpha, adat.Elevator, M*res,'r.')
xlabel('alpha')
ylabel('elevator')
zlabel('CD_dElevator')

% integrate previous fit analytically, plot result
alphas = [];
deltas = [];
drags = [];
for alpha=linspace(-4,8,20)
    for delta=linspace(-6,6,20);
        alphas = [alphas, alpha];
        deltas = [deltas, delta];
        drags = [drags,res(1)*0.5*delta.^2 + res(2)*alpha*delta + res(3)*delta];
    end
end
figure()
plot3(alphas,deltas,drags,'b.')
xlabel('alpha');ylabel('delta');zlabel('CD_Elevator');
grid on

% fit linear Cd_dail
adat.Aileron=-adat.Aileron;
adat.CD_ff_d2=-adat.CD_ff_d2;
M = [adat.Aileron, adat.Beta, ones(n,1)];
res = M\adat.CD_ff_d2;
fprintf('CD_dAileron = %g*Aileron + %g*beta + %g\n', res)
fprintf('CD_Aileron = %g*Aileron^2 + %g*beta*delta + %g*delta\n',res(1)*0.5, res(2), res(3))

figure()
plot3(adat.Beta, adat.Aileron, adat.CD_ff_d2,'b.')
hold on
plot3(adat.Beta, adat.Aileron, M*res,'r.')
xlabel('Beta')
ylabel('Aileron')
zlabel('CD_dAileron')

% integrate previous fit analytically, plot result
alphas = [];
deltas = [];
drags = [];
for alpha=linspace(-8,8,20)
    for delta=linspace(-6,6,20);
        alphas = [alphas, alpha];
        deltas = [deltas, delta];
        drags = [drags,res(1)*0.5*delta.^2 + res(2)*alpha*delta + res(3)*delta];
    end
end
figure()
plot3(alphas,deltas,drags,'b.')
xlabel('beta');ylabel('delta');zlabel('CD_Aileron');
grid on

% fit linear Cd_drud
M = [adat.Rudder, adat.Beta, ones(n,1)];
res = M\adat.CD_ff_d4;
fprintf('CD_dRudder = %g*Rudder + %g*beta + %g\n', res)
fprintf('CD_Rudder = %g*Rudder^2 + %g*beta*delta + %g*delta\n',res(1)*0.5, res(2), res(3))

figure()
plot3(adat.Beta, adat.Rudder, adat.CD_ff_d4,'b.')
hold on
plot3(adat.Beta, adat.Rudder, M*res,'r.')
xlabel('Beta')
ylabel('Rudder')
zlabel('CD_dRudder')

% integrate previous fit analytically, plot result
alphas = [];
deltas = [];
drags = [];
for alpha=linspace(-8,8,20)
    for delta=linspace(-6,6,20);
        alphas = [alphas, alpha];
        deltas = [deltas, delta];
        drags = [drags,res(1)*0.5*delta.^2 + res(2)*alpha*delta + res(3)*delta];
    end
end
figure()
plot3(alphas,deltas,drags,'b.')
xlabel('beta');ylabel('delta');zlabel('CD_Rudder');
grid on