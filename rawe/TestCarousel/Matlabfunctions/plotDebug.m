function [] = plotDebug(MPC)


File = [MPC.FileState,'.txt'];
fid = fopen(File);
Out = textscan(fid,MPC.StateForm);
fclose(fid);

for k = 1:size(Out,2)
    eval([MPC.StateLabels{k},' = Out{k};']);
end
% timeProcess = t;
% stateProcess = [x  y  z  dx dy dz e11 e12 e13 e21 e22 e23 e31 e32 e33 w1 w2 w3 r  dr delta ddelta ur up]; 



Fontsize = 18;

figure(101);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
    subplot(2,1,1)
        plot(t,x,'linewidth',2,'color','k');hold on
        plot(t,y,'linewidth',2,'color','b');hold on
        plot(t,z,'linewidth',2,'color','r');hold on
        plot(MPC.Reference(:,1),MPC.Reference(:,2),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,3),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,4),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        
        
        xlabel('t [s]','FontSize',Fontsize);ylabel('position [m]','FontSize',Fontsize)
        legend('x','y','z','x_r','y_r','z_r','location','NorthEastOutside')
          
    subplot(2,1,2)
        plot(t,dx,'linewidth',2,'color','k');hold on
        plot(t,dy,'linewidth',2,'color','b');hold on
        plot(t,dz,'linewidth',2,'color','r');hold on
        plot(MPC.Reference(:,1),MPC.Reference(:,5),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,6),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,7),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        
        xlabel('t [s]','FontSize',Fontsize);ylabel('position [m]','FontSize',Fontsize)
        legend('dx','dy','dz','dx_r','dy_r','dz_r','location','NorthEastOutside')

Fontsize = 18;
figure(102);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
    subplot(3,1,1)
        plot(t,e11,'linewidth',2,'color','k');hold on
        plot(t,e12,'linewidth',2,'color','b');hold on
        plot(t,e13,'linewidth',2,'color','r');hold on
        plot(MPC.Reference(:,1),MPC.Reference(:,8),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,9),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,10),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        
        xlabel('Time (s)','FontSize',Fontsize);ylabel('e1 [-]','FontSize',Fontsize)
        legend('e11','e12','e13','e11_r','e12_r','e13_r','location','NorthEastOutside')

    subplot(3,1,2)
        plot(t,e21,'linewidth',2,'color','k');hold on
        plot(t,e22,'linewidth',2,'color','b');hold on
        plot(t,e23,'linewidth',2,'color','r');hold on
        plot(MPC.Reference(:,1),MPC.Reference(:,11),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,12),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,13),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        
        xlabel('Time (s)','FontSize',Fontsize);ylabel('e2 [-]','FontSize',Fontsize)
        legend('e21','e22','e23','e21_r','e22_r','e23_r','location','NorthEastOutside')

    subplot(3,1,3)
        plot(t,e31,'linewidth',2,'color','k');hold on
        plot(t,e32,'linewidth',2,'color','b');hold on
        plot(t,e33,'linewidth',2,'color','r');hold on
        plot(MPC.Reference(:,1),MPC.Reference(:,14),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,15),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,16),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        
        xlabel('Time (s)','FontSize',Fontsize);ylabel('e3 [-]','FontSize',Fontsize)
        legend('e31','e32','e33','e31_r','e32_r','e33_r','location','NorthEastOutside')


Fontsize = 18;
figure(103);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
        plot(t,w1,'linewidth',2,'color','k');hold on
        plot(t,w2,'linewidth',2,'color','b');hold on
        plot(t,w3,'linewidth',2,'color','r');hold on
        plot(MPC.Reference(:,1),MPC.Reference(:,17),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,18),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.Reference(:,1),MPC.Reference(:,19),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        
        xlabel('Time (s)','FontSize',Fontsize);ylabel('w [-]','FontSize',Fontsize)
        legend('w1','w2','w3','w1_r','w2_r','w3_r','location','NorthEastOutside')

Fontsize = 18;
figure(104);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
% subplot(3,1,1)
%         plot(t,r,'linewidth',2,'color','k');hold on
%         plot(t,dr,'linewidth',2,'color','b');hold on
%         plot(MPC.Reference(:,1),MPC.Reference(:,20),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
%         plot(MPC.Reference(:,1),MPC.Reference(:,21),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
%         
%         xlabel('Time (s)','FontSize',Fontsize);
%         legend('r','dr','r_r','dr_r','location','NorthEastOutside')
        
subplot(2,1,1)
        plot(t,delta,'linewidth',2,'color','k');hold on
        plot(MPC.Reference(:,1),MPC.Reference(:,20),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        
        xlabel('Time (s)','FontSize',Fontsize);ylabel('delta [rad]','FontSize',Fontsize)
        legend('\delta','\delta_r','location','NorthEastOutside')
        
subplot(2,1,2)
        plot(t,ddelta,'linewidth',2,'color','k');hold on
        plot(MPC.Reference(:,1),MPC.Reference(:,21),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        
        xlabel('Time (s)','FontSize',Fontsize);ylabel('ddelta [rad/s]','FontSize',Fontsize)
        legend('d\delta','d\delta_r','location','NorthEastOutside')
