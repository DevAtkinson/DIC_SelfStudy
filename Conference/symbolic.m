syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 z
P = sym('P', [1 8])
% P=[P1 P2 P3 P4 P5 P6]
% W=[x0+dx.*(1+P(2))+P(3).*dy+P(1);
% y0+dy.*(1+P(6))+P(5).*dx+P(4)]

% dP(:,1)=diff(W,P1)
% dP(:,2)=diff(W,P2)
% dP(:,3)=diff(W,P3)
% dP(:,4)=diff(W,P4)
% dP(:,5)=diff(W,P5)
% dP(:,6)=diff(W,P6)

% WW=[(1+P(2)), P(3), P(1);
% P(5), (P(6)+1), P(4);
% 0 0 1]
% X=[dx; dy; 1]
% I=WW*X

% diff(I,P1)

% dI(:,1)=diff(I,P1)
% dI(:,2)=diff(I,P2)
% dI(:,3)=diff(I,P3)
% dI(:,4)=diff(I,P4)
% dI(:,5)=diff(I,P5)
% dI(:,6)=diff(I,P6)

% A=[x0+dx.*(1)+P(3)*sin(dx);
% y0+dy]

% diff(A,P3)

B=[(1+P(2)), P(3),P7,0, P(1);
P(5), (P(6)+1),0,P8, P(4);
0 0 1 0 0;
0 0 0 1 0;
0 0 0 0 1]
eqn=B(1)==z
zz=solve(eqn,P2)

X=[dx;dy;sin(dx);sin(dy);1];

symbolic_warp(B,X)
% aa=matlabFunction(a)
% bb=matlabFunction(b)
% bbb=b(0.5,0.25)