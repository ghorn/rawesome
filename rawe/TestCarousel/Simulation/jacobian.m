function [A,B,C,D] = jacobian(x,P)

dx=1e-6;
A=[];
for i=1:length(x)
    dxv=zeros(size(x));
    dxv(i)=dx;
    A(:,i) = (  Model_integ_ACADO(0,x+dxv/2,'LQR',P) - Model_integ_ACADO(0,x-dxv/2,'LQR',P)  )/dx;
%     Jx(:,i)=(kite0ode(x+dxv/2,u,t)-kite0ode(x-dxv/2,u,t))/dx;
end

u = P.tu(1,2:end);
% keyboard

du=1e-6;
B=[];
for i=1:length(u)
    duv=zeros(size(u));
    duv(i)=du;
    P.tu(1,2:end) = u+duv/2;
    f1 = Model_integ_ACADO(0,x,'LQR',P);
    P.tu(1,2:end) = u-duv/2;
    f2 = Model_integ_ACADO(0,x,'LQR',P);
    B(:,i) = (  f1 - f2  )/du;
end

dx=1e-6;
C=[];
for i=1:length(x)
    dxv=zeros(size(x));
    dxv(i)=dx;
    C(:,i) = (  Model_integ_ACADO(0,x+dxv/2,'kalman',P) - Model_integ_ACADO(0,x-dxv/2,'kalman',P)  )/dx;
%     Jx(:,i)=(kite0ode(x+dxv/2,u,t)-kite0ode(x-dxv/2,u,t))/dx;
end

% u = P.tu(1,2:end);
% keyboard

du=1e-6;
D=[];
for i=1:length(u)
    duv=zeros(size(u));
    duv(i)=du;
    P.tu(1,2:end) = u+duv/2;
    f1 = Model_integ_ACADO(0,x,'kalman',P);
    P.tu(1,2:end) = u-duv/2;
    f2 = Model_integ_ACADO(0,x,'kalman',P);
    D(:,i) = (  f1 - f2  )/du;
end


end