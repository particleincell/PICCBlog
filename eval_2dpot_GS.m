%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potential solver for a particle-in-cell example program
% Based on the Gauss-Seidel method
%
% For more, visit http://www.particleincell.com/2010/es-pic-method/
% and http://www.particleincell.com/2011/particle-in-cell-example/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x]=eval_2dpot_GS(phi)
global A den n0 phi0 Te phi_p box EPS0 QE;

tol = 0.1;      %solver tolerance

%get nx from size of density
nx = size(den,1);
ny = size(den,2);
nn = numel(den);

%convert density and potential into column vectors
b0 = reshape(den, numel(den),1);
x = reshape(phi, numel(phi),1);

%solve
for it=1:2000
    
    %recalculate rhs
    b = b0 - n0*exp((x-phi0)/Te); %add boltzmann term for electrons
    b = -b*QE/EPS0;
 
    %set boundaries
    b(1:nx) = 0;                 %zero electric field on y=0;
    b(nn-nx+1:nn) = 0;           %zero electric field on y=L;
    b(nx:nx:nn) = 0;             %zero electric field on x=L;
    b(1:nx:nn) = phi0;           %fixed potential on x=0;

    %set potential on fixed nodes
    for j=box(2,1):box(2,2)
        b([box(1,1):box(1,2)]+(j-1)*nx)=ones(box(1,2)-box(1,1)+1,1)*phi_p;      %wall potential
    end

    %update nodes
	for i=1:nn
		x(i)=(b(i) - A(i,1:i-1)*x(1:i-1) - A(i,i+1:nn)*x(i+1:nn))/A(i,i);
    end
    
    %compute residue to check for convergence, do only every 10 iterations
    if mod(it,10)==0
        R=norm(b-A*x);      %residue
        if (R<=tol)         %converged
            %disp(sprintf('  GS converged in %i iterations with norm %g',it,R));
            break;
        end
    end
end

%check if the solver converged to the specified tolerance
if (R>tol)
    disp('  GS failed to converge!!');
end

%return solution as a nx*ny array
x=reshape(x,nx,ny);
