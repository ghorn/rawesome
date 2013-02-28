function [] = plotFigures(MHE,MPC,Sim)

%%%%%%%%%%%% PAPER PLOTS %%%%%%%%%%%%%
    
        
Fontsize = 18;
figure(1);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
    subplot(2,1,1)
        plot(Sim.time,Sim.state(:,1),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.state(:,2),'linewidth',2,'color','b');hold on
        plot(Sim.time,Sim.state(:,3),'linewidth',2,'color','r');hold on
        plot(MPC.refT,MPC.refX(:,1),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,2),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,3),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        plot(MHE.time,MHE.state(:,1),'linestyle','none','marker','o','color','k','markersize',8);hold on
        plot(MHE.time,MHE.state(:,2),'linestyle','none','marker','o','color','b','markersize',8);hold on
        plot(MHE.time,MHE.state(:,3),'linestyle','none','marker','o','color','r','markersize',8);hold on
%         plot(MHE.measT,MHE.measM(:,4),'linewidth',2,'color','k','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,5),'linewidth',2,'color','b','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,6),'linewidth',2,'color','r','linestyle','none','marker','*','markersize',5);
        
        grid on;
        xlabel('t [s]','FontSize',Fontsize);ylabel('position [m]','FontSize',Fontsize)
%         legend('x','y','z','x_r','y_r','z_r','x_{MHE}','y_{MHE}','z_{MHE}','x_m','y_m','z_m','location','NorthEastOutside')
        legend('x','y','z','x_r','y_r','z_r','x_{MHE}','y_{MHE}','z_{MHE}','location','NorthEastOutside')
        
    subplot(2,1,2)
        plot(Sim.time,Sim.state(:,4),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.state(:,5),'linewidth',2,'color','b');hold on
        plot(Sim.time,Sim.state(:,6),'linewidth',2,'color','r');hold on
        plot(MPC.refT,MPC.refX(:,4),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,5),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,6),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        plot(MHE.time,MHE.state(:,4),'linestyle','none','marker','o','color','k','markersize',8);hold on
        plot(MHE.time,MHE.state(:,5),'linestyle','none','marker','o','color','b','markersize',8);hold on
        plot(MHE.time,MHE.state(:,6),'linestyle','none','marker','o','color','r','markersize',8);hold on
%         plot(MHE.measT,MHE.measM(:,7),'linewidth',2,'color','k','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,8),'linewidth',2,'color','b','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,9),'linewidth',2,'color','r','linestyle','none','marker','*','markersize',5);
        
        grid on;
        xlabel('t [s]','FontSize',Fontsize);ylabel('speed [m]','FontSize',Fontsize)
%         legend('dx','dy','dz','dx_r','dy_r','dz_r','dx_{MHE}','dy_{MHE}','dz_{MHE}','dx_m','dy_m','dz_m','location','NorthEastOutside')
        legend('dx','dy','dz','dx_r','dy_r','dz_r','dx_{MHE}','dy_{MHE}','dz_{MHE}','location','NorthEastOutside')

        
Fontsize = 18;
figure(2);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
    subplot(3,1,1)
        plot(Sim.time,Sim.state(:,7),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.state(:,8),'linewidth',2,'color','b');hold on
        plot(Sim.time,Sim.state(:,9),'linewidth',2,'color','r');hold on
        plot(MPC.refT,MPC.refX(:,7),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,8),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,9),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        plot(MHE.time,MHE.state(:,7),'linestyle','none','marker','o','color','k','markersize',8);hold on
        plot(MHE.time,MHE.state(:,8),'linestyle','none','marker','o','color','b','markersize',8);hold on
        plot(MHE.time,MHE.state(:,9),'linestyle','none','marker','o','color','r','markersize',8);hold on
%         plot(MHE.measT,MHE.measM(:,10-3),'linewidth',2,'color','k','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,11-3),'linewidth',2,'color','b','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,12-3),'linewidth',2,'color','r','linestyle','none','marker','*','markersize',5);
        
        grid on;
        xlabel('Time (s)','FontSize',Fontsize);ylabel('e1 [-]','FontSize',Fontsize)
%         legend('e11','e12','e13','e11_r','e12_r','e13_r','e11_{MHE}','e12_{MHE}','e13_{MHE}','e11_m','e12_m','e13_m','location','NorthEastOutside')
        legend('e11','e12','e13','e11_r','e12_r','e13_r','e11_{MHE}','e12_{MHE}','e13_{MHE}','location','NorthEastOutside')

    subplot(3,1,2)
        plot(Sim.time,Sim.state(:,10),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.state(:,11),'linewidth',2,'color','b');hold on
        plot(Sim.time,Sim.state(:,12),'linewidth',2,'color','r');hold on
        plot(MPC.refT,MPC.refX(:,10),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,11),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,12),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        plot(MHE.time,MHE.state(:,10),'linestyle','none','marker','o','color','k','markersize',8);hold on
        plot(MHE.time,MHE.state(:,11),'linestyle','none','marker','o','color','b','markersize',8);hold on
        plot(MHE.time,MHE.state(:,12),'linestyle','none','marker','o','color','r','markersize',8);hold on
%         plot(MHE.measT,MHE.measM(:,13-3),'linewidth',2,'color','k','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,14-3),'linewidth',2,'color','b','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,15-3),'linewidth',2,'color','r','linestyle','none','marker','*','markersize',5);
        
        grid on;
        xlabel('Time (s)','FontSize',Fontsize);ylabel('e2 [-]','FontSize',Fontsize)
%         legend('e21','e22','e23','e21_r','e22_r','e23_r','e21_{MHE}','e22_{MHE}','e23_{MHE}','e21_m','e22_m','e23_m','location','NorthEastOutside')
        legend('e21','e22','e23','e21_r','e22_r','e23_r','e21_{MHE}','e22_{MHE}','e23_{MHE}','location','NorthEastOutside')

    subplot(3,1,3)
        plot(Sim.time,Sim.state(:,13),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.state(:,14),'linewidth',2,'color','b');hold on
        plot(Sim.time,Sim.state(:,15),'linewidth',2,'color','r');hold on
        plot(MPC.refT,MPC.refX(:,13),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,14),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,15),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
        plot(MHE.time,MHE.state(:,13),'linestyle','none','marker','o','color','k','markersize',8);hold on
        plot(MHE.time,MHE.state(:,14),'linestyle','none','marker','o','color','b','markersize',8);hold on
        plot(MHE.time,MHE.state(:,15),'linestyle','none','marker','o','color','r','markersize',8);hold on
%         plot(MHE.measT,MHE.measM(:,16-3),'linewidth',2,'color','k','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,17-3),'linewidth',2,'color','b','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,18-3),'linewidth',2,'color','r','linestyle','none','marker','*','markersize',5);
        
        grid on;
        xlabel('Time (s)','FontSize',Fontsize);ylabel('e3 [-]','FontSize',Fontsize)
%         legend('e31','e32','e33','e31_r','e32_r','e33_r','e31_{MHE}','e32_{MHE}','e33_{MHE}','e31_m','e32_m','e33_m','location','NorthEastOutside')
        legend('e31','e32','e33','e31_r','e32_r','e33_r','e31_{MHE}','e32_{MHE}','e33_{MHE}','location','NorthEastOutside')


Fontsize = 18;
figure(3);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
%         plot(Sim.time,Sim.state(:,16),'linewidth',2,'color','k');hold on
%         plot(Sim.time,Sim.state(:,17),'linewidth',2,'color','b');hold on
%         plot(Sim.time,Sim.state(:,18),'linewidth',2,'color','r');hold on
        plot(Sim.time,Sim.ddX(:,1),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.ddX(:,2),'linewidth',2,'color','b');hold on
        plot(Sim.time,Sim.ddX(:,3),'linewidth',2,'color','r');hold on
        plot(MPC.refT,MPC.refX(:,16),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,17),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,18),'linewidth',2,'color','r','linestyle','none','marker','.','markersize',15);
%         plot(MHE.time,MHE.state(:,16),'linestyle','none','marker','o','color','k','markersize',8);hold on
%         plot(MHE.time,MHE.state(:,17),'linestyle','none','marker','o','color','b','markersize',8);hold on
%         plot(MHE.time,MHE.state(:,18),'linestyle','none','marker','o','color','r','markersize',8);hold on
        plot(MHE.time,MHE.ddXIMU(:,1),'linestyle','none','marker','o','color','k','markersize',8);hold on
        plot(MHE.time,MHE.ddXIMU(:,2),'linestyle','none','marker','o','color','b','markersize',8);hold on
        plot(MHE.time,MHE.ddXIMU(:,3),'linestyle','none','marker','o','color','r','markersize',8);hold on
        plot(MHE.measT,MHE.measM(:,13),'linewidth',2,'color','k','linestyle','none','marker','*','markersize',5);
        plot(MHE.measT,MHE.measM(:,14),'linewidth',2,'color','b','linestyle','none','marker','*','markersize',5);
        plot(MHE.measT,MHE.measM(:,15),'linewidth',2,'color','r','linestyle','none','marker','*','markersize',5);
        
        grid on;
        xlabel('Time (s)','FontSize',Fontsize);ylabel('w [-]','FontSize',Fontsize)
        legend('w1','w2','w3','w1_r','w2_r','w3_r','w1_{MHE}','w2_{MHE}','w3_{MHE}','w1_m','w2_m','w3_m','location','NorthEastOutside')

Fontsize = 18;
figure(4);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
subplot(3,1,1)
        plot(Sim.time,Sim.state(:,19),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.state(:,20),'linewidth',2,'color','b');hold on
        plot(MPC.refT,MPC.refX(:,19),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MPC.refT,MPC.refX(:,20),'linewidth',2,'color','b','linestyle','none','marker','.','markersize',15);
        plot(MHE.time,MHE.state(:,19),'linestyle','none','marker','o','color','k','markersize',8);hold on
        plot(MHE.time,MHE.state(:,20),'linestyle','none','marker','o','color','b','markersize',8);hold on
        plot(MHE.measT,MHE.measM(:,19),'linewidth',2,'color','k','linestyle','none','marker','*','markersize',5);
%         plot(MHE.measT,MHE.measM(:,23-3),'linewidth',2,'color','b','linestyle','none','marker','*','markersize',5);
        
        grid on;
        xlabel('Time (s)','FontSize',Fontsize);
        legend('\delta','d\delta','\delta_r','d\delta_r','\delta_{MHE}','d\delta_{MHE}','\delta_m','location','NorthEastOutside')
        
subplot(3,1,2)
        plot(Sim.time,Sim.state(:,21),'linewidth',2,'color','k');hold on
        plot(MPC.refT,MPC.refX(:,21),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MHE.time,MHE.state(:,21),'linestyle','none','marker','o','color','k','markersize',8);hold on
        
        grid on;
        xlabel('Time (s)','FontSize',Fontsize);ylabel('delta [rad]','FontSize',Fontsize)
        legend('ur','ur_r','ur_{MHE}','location','NorthEastOutside')
        
subplot(3,1,3)
        plot(Sim.time,Sim.state(:,22),'linewidth',2,'color','k');hold on
        plot(MPC.refT,MPC.refX(:,22),'linewidth',2,'color','k','linestyle','none','marker','.','markersize',15);
        plot(MHE.time,MHE.state(:,22),'linestyle','none','marker','o','color','k','markersize',8);hold on
        
        grid on;
        xlabel('Time (s)','FontSize',Fontsize);ylabel('ddelta [rad/s]','FontSize',Fontsize)
        legend('up','up_r','up_{MHE}','location','NorthEastOutside')

        
Fontsize = 18;
figure(5);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
    subplot(6,1,1)
        plot(MHE.time,MHE.constraints(:,1),'linewidth',2,'color','k','linestyle','none','marker','o');hold on
        plot(MHE.time,MHE.constraints(:,2),'linewidth',2,'color','b','linestyle','none','marker','o');hold on
        
        grid on;
%         xlabel('Time (s)','FontSize',Fontsize);
        legend('Const','dConst','location','NorthEastOutside')
        title('Constraint violation (MHE)')
        
    subplot(6,1,2)
        plot(MHE.time,MHE.constraints(:,3),'linewidth',2,'color','k','linestyle','none','marker','o');hold on
        plot(MHE.time,MHE.constraints(:,4),'linewidth',2,'color','b','linestyle','none','marker','o');hold on
        plot(MHE.time,MHE.constraints(:,5),'linewidth',2,'color','r','linestyle','none','marker','o');hold on
        plot(MHE.time,MHE.constraints(:,6),'linewidth',2,'color','y','linestyle','none','marker','o');hold on
        plot(MHE.time,MHE.constraints(:,7),'linewidth',2,'color','g','linestyle','none','marker','o');hold on
        plot(MHE.time,MHE.constraints(:,8),'linewidth',2,'color','m','linestyle','none','marker','o');hold on
        
        grid on;
%         xlabel('Time (s)','FontSize',Fontsize);
        legend('e_1e_1 - 1','e_1e_2','e_1e_3','e_2e_2 - 1','e_2e_3','e_3e_3 - 1','location','NorthEastOutside')
        
    subplot(6,1,3)
        semilogy(MPC.time,abs(MPC.constraints(:,1)),'linewidth',2,'color','k','linestyle','none','marker','*');hold on
        semilogy(MPC.time,abs(MPC.constraints(:,2)),'linewidth',2,'color','b','linestyle','none','marker','*');hold on
        
        grid on;
        ylim([1e-10 1e-0])
%         xlabel('Time (s)','FontSize',Fontsize);
        legend('Const','dConst','location','NorthEastOutside')
        title('Constraint violation (MPC)')
        
    subplot(6,1,4)
        semilogy(MPC.time,abs(MPC.constraints(:,3)),'linewidth',2,'color','k','linestyle','none','marker','*');hold on
        semilogy(MPC.time,abs(MPC.constraints(:,4)),'linewidth',2,'color','b','linestyle','none','marker','*');hold on
        semilogy(MPC.time,abs(MPC.constraints(:,5)),'linewidth',2,'color','r','linestyle','none','marker','*');hold on
        semilogy(MPC.time,abs(MPC.constraints(:,6)),'linewidth',2,'color','y','linestyle','none','marker','*');hold on
        semilogy(MPC.time,abs(MPC.constraints(:,7)),'linewidth',2,'color','g','linestyle','none','marker','*');hold on
        semilogy(MPC.time,abs(MPC.constraints(:,8)),'linewidth',2,'color','m','linestyle','none','marker','*');hold on
        
        grid on;
%         xlabel('Time (s)','FontSize',Fontsize);
        legend('e_1e_1 - 1','e_1e_2','e_1e_3','e_2e_2 - 1','e_2e_3','e_3e_3 - 1','location','NorthEastOutside')

    subplot(6,1,5)
        plot(Sim.time,Sim.constraints(:,1),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.constraints(:,2),'linewidth',2,'color','b');hold on
        
        grid on;
%         xlabel('Time (s)','FontSize',Fontsize);
        legend('Const','dConst','location','NorthEastOutside')
        title('Constraint violation (Simulation)')
        
    subplot(6,1,6)
        plot(Sim.time,Sim.constraints(:,3),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.constraints(:,4),'linewidth',2,'color','b');hold on
        plot(Sim.time,Sim.constraints(:,5),'linewidth',2,'color','r');hold on
        plot(Sim.time,Sim.constraints(:,6),'linewidth',2,'color','y');hold on
        plot(Sim.time,Sim.constraints(:,7),'linewidth',2,'color','g');hold on
        plot(Sim.time,Sim.constraints(:,8),'linewidth',2,'color','m');hold on
        
        grid on;
%         xlabel('Time (s)','FontSize',Fontsize);
        legend('e_1e_1 - 1','e_1e_2','e_1e_3','e_2e_2 - 1','e_2e_3','e_3e_3 - 1','location','NorthEastOutside')

Fontsize = 18;
figure(6);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
        plot(Sim.time,Sim.ddX(:,4),'linewidth',2,'color','k');hold on
        plot(Sim.time,Sim.ddX(:,5),'linewidth',2,'color','b');hold on
        plot(Sim.time,Sim.ddX(:,6),'linewidth',2,'color','r');hold on
        plot(MHE.time,MHE.ddXIMU(:,4),'linestyle','none','marker','o','color','k','markersize',8);hold on
        plot(MHE.time,MHE.ddXIMU(:,5),'linestyle','none','marker','o','color','b','markersize',8);hold on
        plot(MHE.time,MHE.ddXIMU(:,6),'linestyle','none','marker','o','color','r','markersize',8);hold on
        plot(MHE.measT,MHE.measM(:,16),'linewidth',2,'color','k','linestyle','none','marker','*','markersize',5);
        plot(MHE.measT,MHE.measM(:,17),'linewidth',2,'color','b','linestyle','none','marker','*','markersize',5);
        plot(MHE.measT,MHE.measM(:,18),'linewidth',2,'color','r','linestyle','none','marker','*','markersize',5);

        grid on;
        xlabel('Time (s)','FontSize',Fontsize);
        legend('ddx','ddy','ddz','ddxIMU_{MHE}','ddyIMU_{MHE}','ddzIMU_{MHE}','ddxIMU_m','ddyIMU_m','ddzIMU_m','location','NorthEastOutside')

rA = 1.085;
Fontsize = 18;
% 3D trajectory plot
figure(7);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
        plot3(Sim.state(:,1),Sim.state(:,2),Sim.state(:,3),'linewidth',2,'color','k');hold on
        plot3(MHE.state(:,1),MHE.state(:,2),MHE.state(:,3),'linestyle','none','marker','o','color','k','markersize',10);hold on
        plot3(0,0,0,'color','k','linewidth',1,'marker','.','color','k','markersize',20);hold on
%         plot3(-rA*sin(Sim.state(:,19)),rA*cos(Sim.state(:,19)),0*Sim.state(:,3),'color','k','linewidth',1);hold on

%         plot3(MPC.Xguess(:,1),MPC.Xguess(:,2),MPC.Xguess(:,3),'linewidth',2,'color','b');
%         plot3(MHE.Xref0(:,1),MHE.Xref0(:,2),MHE.Xref0(:,3),'color','k','linestyle','none','marker','.','markersize',15);
grid on;
axis equal;
view(-114,10)
xlabel('x','FontSize',Fontsize);
ylabel('y','FontSize',Fontsize);
zlabel('z','FontSize',Fontsize);
% legend('xyz','xyz_{MHE}','xyz_{arm tip}','location','NorthEastOutside')
        
export_fig ./trajectory.pdf -pdf -transparent -painters

    %%%%%%%%%%%% END PAPER PLOTS %%%%%%%%%%%%%
            
    drawnow
    
    
    
    if 0
            
Fontsize = 18;
figure(5);clf
whitebg([1.0 1.0 1.0])
set(gcf,'Color',[1 1 1])
    subplot(2,1,1)
        semilogy(MHE.time,abs(MHE.constraints(:,1)),'linewidth',2,'color','k','linestyle','none','marker','o','markersize',10);hold on
        semilogy(MHE.time,abs(MHE.constraints(:,2)),'linewidth',2,'color','r','linestyle','none','marker','o','markersize',10);hold on
        semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints(2:end,1)),'linewidth',2,'color','k','linestyle','none','marker','*');hold on
        semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints(2:end,2)),'linewidth',2,'color','b','linestyle','none','marker','*');hold on
        
        
%         semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints0(2:end,1)),'linewidth',2,'color','k','linestyle','none','marker','*');hold on
%         semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints0(2:end,2)),'linewidth',2,'color','b','linestyle','none','marker','*');hold on
        
        grid on;
%         xlabel('Time (s)','FontSize',Fontsize);
%         legend('Const_{MHE}','dConst_{MHE}','Const_{MPC}','dConst_{MPC}','location','NorthEastOutside')
%         title('Constraint violation (MHE)')
        
    subplot(2,1,2)
        semilogy(MHE.time,abs(MHE.constraints(:,3)),'linewidth',2,'color','k','linestyle','none','marker','o','markersize',10);hold on
        semilogy(MHE.time,abs(MHE.constraints(:,4)),'linewidth',2,'color','k','linestyle','none','marker','o','markersize',10);hold on
        semilogy(MHE.time,abs(MHE.constraints(:,5)),'linewidth',2,'color','k','linestyle','none','marker','o','markersize',10);hold on
        semilogy(MHE.time,abs(MHE.constraints(:,6)),'linewidth',2,'color','k','linestyle','none','marker','o','markersize',10);hold on
        semilogy(MHE.time,abs(MHE.constraints(:,7)),'linewidth',2,'color','k','linestyle','none','marker','o','markersize',10);hold on
        semilogy(MHE.time,abs(MHE.constraints(:,8)),'linewidth',2,'color','k','linestyle','none','marker','o','markersize',10);hold on
        semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints(2:end,3)),'linewidth',2,'color','k','linestyle','none','marker','*');hold on
        semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints(2:end,4)),'linewidth',2,'color','b','linestyle','none','marker','*');hold on
        semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints(2:end,5)),'linewidth',2,'color','r','linestyle','none','marker','*');hold on
        semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints(2:end,6)),'linewidth',2,'color','y','linestyle','none','marker','*');hold on
        semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints(2:end,7)),'linewidth',2,'color','g','linestyle','none','marker','*');hold on
        semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints(2:end,8)),'linewidth',2,'color','m','linestyle','none','marker','*');hold on
        
        
%         semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints0(2:end,3)),'linewidth',2,'color','k','linestyle','none','marker','*');hold on
%         semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints0(2:end,4)),'linewidth',2,'color','b','linestyle','none','marker','*');hold on
%         semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints0(2:end,5)),'linewidth',2,'color','r','linestyle','none','marker','*');hold on
%         semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints0(2:end,6)),'linewidth',2,'color','y','linestyle','none','marker','*');hold on
%         semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints0(2:end,7)),'linewidth',2,'color','g','linestyle','none','marker','*');hold on
%         semilogy(MPC.time(2:end)-MPC.time(2),abs(MPC.constraints0(2:end,8)),'linewidth',2,'color','m','linestyle','none','marker','*');hold on
        
        grid on;
%         xlabel('Time (s)','FontSize',Fontsize);
%         legend('e_1e_1 - 1','e_1e_2','e_1e_3','e_2e_2 - 1','e_2e_3','e_3e_3 - 1','location','NorthEastOutside')
        
    end