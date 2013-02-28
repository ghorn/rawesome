function MHE = MHEmeasurements(MHE,Sim)

%%% Generate measurements
MHE.t_index = 1;


MHE.state_index = [16:18];
MHE.state_index2 = [19 21:22];


P.W0 = 0;           % Measurements always built with no wind
P.r = Sim.r;

% keyboard


imu = load('Measurements/IMUData.dat').';
LEDpos = load('Measurements/absData.dat').';
controls = load('Measurements/controlData.dat').';

LEDpos(65,2:end) = LEDpos(64,2:end);
LEDpos(75,2:end) = LEDpos(74,2:end);
LEDpos(87,2:end) = LEDpos(86,2:end);
LEDpos(121,2:end) = LEDpos(120,2:end);

% tIMU = imu(:,1)-imu(1,1);
w1 = imu(:,1);
w2 = imu(:,2);
w3 = imu(:,3);
ax = imu(:,4);
ay = imu(:,5);
az = imu(:,6);

% tCAM = LEDpos(:,1)-LEDpos(1,1);
delta = LEDpos(:,1);
C1L1 = LEDpos(:,2:3);
C1L2 = LEDpos(:,4:5);
C1L3 = LEDpos(:,6:7);
C2L1 = LEDpos(:,8:9);
C2L2 = LEDpos(:,10:11);
C2L3 = LEDpos(:,12:13);
Markers = LEDpos(:,2:13);

% figure
% plot(LEDpos(:,2),LEDpos(:,3),'r.')
% figure
% plot(LEDpos(:,4),LEDpos(:,5),'g.')
% figure
% plot(LEDpos(:,6),LEDpos(:,7),'b.')

% Markers = [];
% for k = MHE.t_index:MHE.Nref-1+MHE.t_index
%     
%     P.tu = [0 MHE.Uref(k,:)];
%     P.delta = MHE.Xref(k,21);
%     P.ddelta = MHE.Xref(k,22);
%     P.dddelta = MHE.Uref(k,1);
%     P.r = Sim.r;
%     
%     out = Model_integ_ACADO( 0, MHE.Xref(k,:), 'markers', P );
%     Markers = [Markers;out(1:12).';];
% end

% figure()
% clf
% hold on
% plot3(Markers(:,1),Markers(:,2),Markers(:,3),'.k')
% plot3(Markers(:,4),Markers(:,5),Markers(:,6),'.b')
% plot3(Markers(:,7),Markers(:,8),Markers(:,9),'.r')
% axis equal
% grid on
% keyboard

timeIMU = (0:size(w1,1)-1)*0.002;

index = [0:size(Markers,1)-1]*50+1;

time = (0:size(Markers,1)-1)*0.1;
% timeIMU(index) - time

IMU0 = [w1 w2 w3 ax ay az];

% keyboard
IMU = IMU0(1,:);
for k = 2:length(index)-1
    IMU1 = sum(IMU0(index(k)-10:index(k)+10,:))/21;
    IMU = [IMU; IMU1];
end
IMU = [IMU; IMU0(index(end),:)];
% IMU = IMU0(index,:);

delta = delta - delta(1);

MHE.Measurements = [time.' Markers IMU delta zeros(size(delta,1),2) zeros(size(delta,1),3)];
% measurements = [Markers MHE.Xref(MHE.t_index:MHE.Nref-1+MHE.t_index,MHE.state_index) MHE.ddXref(MHE.t_index:MHE.Nref-1+MHE.t_index,:) MHE.Xref(MHE.t_index:MHE.Nref-1+MHE.t_index,MHE.state_index2) MHE.Uref(MHE.t_index:MHE.Nref-1+MHE.t_index,1:3)];
% measurements = [MHE.ddXref(MHE.t_index:MHE.Ncvp+MHE.t_index,:) MHE.Xref(MHE.t_index:MHE.Ncvp+MHE.t_index,MHE.state_index) MHE.Uref(MHE.t_index:MHE.Ncvp+MHE.t_index,1:4)];
% measurements = [MHE.ddXref(MHE.t_index:MHE.Ncvp+MHE.t_index,:) Markers MHE.Xref(MHE.t_index:MHE.Ncvp+MHE.t_index,MHE.state_index) MHE.Uref(MHE.t_index:MHE.Ncvp+MHE.t_index,1:4)];

% isNoise = Sim.Noise.is;
% factorNoise = Sim.Noise.factor;
% 
% MHE.Scaling = [ones(1,12)*1e3  ones(1,3) ones(1,3)*10 0 0 0 0e-4 0e-4 0e-4];
% % MHE.Scaling = [ ones(1,3)*10 ones(1,12)*1e3  ones(1,3) 0 0 0 0 0 0 0e-4 0e-4 0e-4 0e-4];
% 
% Noise = randn(size(measurements))*factorNoise*diag(MHE.Scaling);
% measurements = measurements + Noise*isNoise;
% 
% MHE.Measurements = [MHE.Tref(1:MHE.Nref,:) measurements];
% % MHE.Measurements = [MHE.Tref(1:MHE.Nref,:) measurements zeros(MHE.Nref,6)]; % add the ficticious forces and torques
% % MHE.Measurements = [MHE.Tref(1:MHE.Ncvp+1,:) measurements zeros(MHE.Ncvp+1,6)]; % add the ficticious forces and torques

formMeas = '';
for k = 1:size(MHE.Measurements,2)
    formMeas = [formMeas,' %6.16e'];
end
MHE.formMeas = [formMeas,' \n'];

!rm Measurements.txt
fid = fopen('Measurements.txt', 'w');
fprintf(fid,MHE.formMeas,MHE.Measurements(1:MHE.Ncvp+1,:)');
fclose(fid);

end

