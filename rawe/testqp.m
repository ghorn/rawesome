clc
clear all
s = load('hessL.mat');

eig(s.hessL);
H = s.hessL;
f = s.gradF;
uba = s.uba;
lba = s.lba;
lb = s.lbx;
ub = s.ubx;
x0 = s.x0;

beq = [];
Aeq = [];
for k=1:length(x0)
    if lb(k)==ub(k)
        beq = [beq;lb(k)];
        blah = zeros(1,length(x0));
        blah(k) = 1;
        Aeq = [Aeq;blah];
        lb(k) = -1e6;
        ub(k) = 1e6;
    end
end

Aineq = [];
bineq = [];
neqs = 0;
nineqs = 0;
for k=1:size(s.jacobG,1)
    if lba(k)==uba(k)
        Aeq = [Aeq; s.jacobG(k,:)];
        beq = [beq; uba(k)];
        neqs = neqs + 1;
    else
        Aineq = [Aineq;s.jacobG(k,:);-s.jacobG(k,:)];
        bineq = [bineq; uba(k); -lba(k)];
        nineqs = nineqs + 1;
    end
end

options = optimset;
% options.MaxIter = 10000;
options.Algorithm = 'active-set';
dx = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
%dx = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,dx,options);
% (dx - s.xopt)'
dx' - s.xopt'

disp 'length(x0)'
length(x0)

disp 'size(Aineq)'
size(Aineq)

disp 'size(Aeq)'
size(Aeq)

fprintf('neqs: %d, nineqs: %d\n',neqs,nineqs)
