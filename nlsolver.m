%Demo non-linear Poisson solver
% we are solving:
% d^2phi/(dh^2) = -rho/eps0*(ni-n0*exp((phi-phi0)/kTe))
% with phi=-5, 5 on left and right boundary
% See http://www.particleincell.com/2012/nonlinear-poisson-solver/ 
% for additional info

%clear screen
clc

%constants
EPS0 = 8.85418782e-12;	
Q = 1.60217646e-19;

%setup coefficients
den0=1e16;
phi0=0;
kTe=5;

%precomputed values
lambda_D = sqrt(EPS0*kTe/(den0*Q));	%for kTe  in eV
dh = lambda_D;						%cell spacing
C = -Q/EPS0*dh*dh;

%setup matrixes
nn=15;			%number of nodes
A=zeros(nn,nn);
fixed_node = zeros(nn,1);
b0=zeros(nn,1);
	
%left boundary
A(1,1)=1;
b0(1)=-5;
fixed_node(1)=1;

%right boundary
A(nn,nn)=1;
b0(nn)=5;
fixed_node(nn)=1;

%internal nodes
for n=2:nn-1
	%FD stencil
	A(n,n-1)=1;
	A(n,n)=-2;
	A(n,n+1)=1;
	
	fixed_node(n)=false;
	b0(n)=C*den0;
end

%initial values
bx = zeros(nn,1);
P = zeros(nn,1);
J = zeros(nn,nn);
x = zeros(nn,1);
y = zeros(nn,1);

%--- Newton Solver ----
for it=1:20

	%1) compute bx
	for n=1:nn
		if (fixed_node(n))
			bx(n)=0;
		else
			bx(n) = -C*den0*exp((x(n)-phi0)/kTe);
		end
	end
	
	%2) b=b0+bx
	b = b0 + bx;
	
	%3) F=Ax-b
	F = A*x-b;
	
	%4) compute P=dbx(i)/dx(j), used in computing J
	for n=1:nn
		if (fixed_node(n))
			P(n)=0;
		else
			P(n) = -C*den0*exp((x(n)-phi0)/kTe)/kTe;
		end
	end
	
	%5) J=df(i)/dx(j), J=A+diag(dbx(i)/dx(j))
	J = A - diag(P);
	
	%6) Solve Jy=F
	y = J\F;
	
	%7) update x
	x = x - y;
	
	%8) compute norm;
	l2 = norm(y);
	if (l2<1e-6)
		disp(sprintf("Converged in %d iterations with norm %g\n",it,l2));	
		break;
	end
end

disp(x');
plot(x);
	
