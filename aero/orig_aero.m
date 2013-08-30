clear all
close all
clc

adat = aerodata

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


% fit quadratic CD alpha
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
M = [(alpha*pi/180).^2, (alpha*pi/180), ones(length(alpha),1)];
res = M\CD;
fprintf('CD = %f*alpha^2 + %f*alpha + %f\n', res)

figure()
plot3(adat.Alpha, adat.Beta, adat.CD_tot,'b.')
hold on
M = [(adat.Alpha*pi/180).^2, (adat.Alpha*pi/180), ones(n,1)];
plot3(adat.Alpha, adat.Beta, M*res,'r.')
xlabel('alpha')
ylabel('beta')
zlabel('CD')


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