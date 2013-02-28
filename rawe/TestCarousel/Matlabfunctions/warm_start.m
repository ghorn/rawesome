function MHE = warm_start(MHE)

% Take the equilibrium point and rotate it back so as to build measurements
% to warm start the MHE

Nref = MHE.Nref;
Ts = MHE.Ts;

x = MHE.X0(1);
y = MHE.X0(2);
z = MHE.X0(3);

dx = MHE.X0(4);
dy = MHE.X0(5);
dz = MHE.X0(6);

e11 = MHE.X0(7);
e12 = MHE.X0(8);
e13 = MHE.X0(9);
e21 = MHE.X0(10);
e22 = MHE.X0(11);
e23 = MHE.X0(12);
e31 = MHE.X0(13);
e32 = MHE.X0(14);
e33 = MHE.X0(15);

w1 = MHE.X0(16);
w2 = MHE.X0(17);
w3 = MHE.X0(18);

r = MHE.X0(19);
dr = MHE.X0(20);

delta = MHE.X0(21);
ddelta = MHE.X0(22);

ur = MHE.X0(23);
up = MHE.X0(24);


R = [e11 e21 e31;
     e12 e22 e32;
     e13 e23 e33];

delta = delta - (Nref-1)*Ts*ddelta;

Xref = zeros(Nref,24);
Uref = zeros(Nref,4);
Tref = zeros(Nref,1);
for k = 1:Nref
    
    Xk = [x;y;z];
    dXk = [dx;dy;dz];
    X0k = [Xk;dXk;reshape(R,9,1);w1;w2;w3;r;dr;delta;ddelta;ur;up].';
    
    
    Xref(k,:) = X0k;
    Uref(k,:) = MHE.U0;
    Tref(k) = Ts*(k-1);
    
    delta = delta + ddelta*Ts;
end

MHE.Tref = Tref;
MHE.Xref = Xref;
MHE.Uref = Uref;

