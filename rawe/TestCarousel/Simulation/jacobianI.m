function [A,B,C,D] = jacobianI(x,P)

dx=1e-100;
A=[];
for k=1:length(x)
    dxv=zeros(size(x));
    dxv(k)=1i*dx;
    A(:,k) = imag((  Model_integ_ACADO(0,x+dxv,'',P)  )/(dx));
%     Jx(:,i)=(kite0ode(x+dxv/2,u,t)-kite0ode(x-dxv/2,u,t))/dx;
end

u = P.tu(1,2:end);
% keyboard

du=1e-100;
B=[];
for k=1:length(u)
    duv=zeros(size(u));
    duv(k)=1i*du;
    P.tu(1,2:end) = u+duv;
    f1 = Model_integ_ACADO(0,x,'',P);
    B(:,k) = imag((  f1  )/(du));
end


dx=1e-100;
C=[];
% for k=1:length(x)
%     dxv=zeros(size(x));
%     dxv(k)=1i*dx;
%     C(:,k) = imag((  Model_integ_ACADO(0,x+dxv,'kalman',P)  )/(dx));
% %     Jx(:,i)=(kite0ode(x+dxv/2,u,t)-kite0ode(x-dxv/2,u,t))/dx;
% end

% u = P.tu(1,2:end);
% keyboard

du=1e-100;
D=[];
% for k=1:length(u)
%     duv=zeros(size(u));
%     duv(k)=1i*du;
%     P.tu(1,2:end) = u+duv;
%     f1 = Model_integ_ACADO(0,x,'kalman',P);
%     D(:,k) = imag((  f1  )/(du));
% end


end
