function MPC = generate_reference(MPC)

%% Read the values and build up an initial guess
z = MPC.Ref.z;
r = MPC.Ref.r;
delta = MPC.Ref.delta;
delta0 = MPC.Ref.delta0;
ddelta = MPC.Ref.ddelta;
Ts = MPC.Ts;
Nref = MPC.Nref;

rA = 1.085; % Arm length
MPC.Params.rA = rA; % Store the value in the structure

r1 = sqrt(r^2 - z^2);
x = r1;
y = 0;

dx = 0; 
dy = 0; 
dz = 0;

e1 = [0;1;0];e2 = [0;0;1];e3 = [1;0;0]; 
R = [e1 e2 e3].';                       

% Rotational speed in aricraft's referential
wE = R.'*[0;0;ddelta];
w1 = wE(1);
w2 = wE(2);
w3 = wE(3);

dr = 0;

ur = 0;
up = 0;

X0 = [x;y;z;dx;dy;dz;reshape(R,1,9).';w1;w2;w3;delta;ddelta;ur;up].';


U0 = zeros(1,4);
X0U0 = [0 X0 U0];

%% Save the file, compile and run the code
!rm eq/EQ_init.txt
save('eq/EQ_init.txt','X0U0','-ascii')

delta = 0;
cd eq
% !make equilibrium
eval(['!/usr/bin/env -u LD_LIBRARY_PATH ./equilibrium ',num2str(z),'     ',num2str(r),'     ',num2str(delta),'     ',num2str(ddelta)]);
cd ..

%% Extract the values and create the reference

formEQ = '[';
for k = 1:size(X0U0,2)
    formEQ = [formEQ,' %f'];
end
formEQ = [formEQ,' ]\n'];
File = 'eq/EQ_params.txt';
fid = fopen(File);
Out = textscan(fid,formEQ);
fclose(fid);


EQLabels = {'t', 'x', 'y', 'z', 'dx', 'dy', 'dz', 'e11', 'e12', 'e13', 'e21', 'e22', 'e23', 'e31', 'e32', 'e33', 'w1', 'w2', 'w3', 'delta', 'ddelta','ur','up', 'dddelta',' ddr', 'dur', 'dup'};

for k = 1:size(Out,2)
    eval([EQLabels{k},' = Out{k};']);
end

P.tu =  [0 0 0 0 0;
         1 0 0 0 0];
P.r = r;

% Stack values in matrices that are gonna be saved in files for ACADO
delta = delta0;
Xref = zeros(Nref,24);
Uref = zeros(Nref,4);
Tref = zeros(Nref,1);
for k = 1:Nref
    
    X0k = [x y z  dx dy dz e11 e12 e13 e21 e22 e23 e31 e32 e33 w1 w2 w3 r dr delta ddelta ur up].';

    Xref(k,:) = X0k;
    Uref(k,:) = zeros(1,4);
    Tref(k) = Ts*(k-1);
    
    delta = delta + ddelta*Ts;
end

MPC.Tref = Tref;
MPC.Xref = Xref;
MPC.Uref = Uref;

