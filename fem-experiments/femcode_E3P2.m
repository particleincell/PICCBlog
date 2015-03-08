% femcode_P2.m
 
% [p,t,b] from distmesh tool
% make sure your matlab path includes the directory where distmesh is installed.
 
fd=@(p) sqrt(sum(p.^2,2))-1;
[p,t]=distmesh2d(fd,@huniform,0.20,[-1,-1;1,1],[]);
be=boundedges(p,t);
b=unique(be);
Nnorm = size(p,1);
 
f=vectorize(inline('0','x','y'));
k=vectorize(inline('1','x','y'));
 
u=vectorize(inline('3*y*x^2-y^3','x','y'));
ux=vectorize(inline('6*y*x','x','y'));
uy=vectorize(inline('3*x^2-3*y^2','x','y'));
 
 
% for quadratic element we need 6 nodes per triangle, so add 1 node to each triangle edge midpoint.
N=size(p,1);T=size(t,1);
M=sparse(N,N);
cnt = N+1;
for e=1:T
  nodes=t(e,:); % row of t = node numbers of the 3 corners of triangle e
  if (M(nodes(1,1),nodes(1,2)) == 0)
     M(nodes(1,1),nodes(1,2)) = cnt;
     M(nodes(1,2),nodes(1,1)) = cnt;
     p(cnt,:) = mean(p([nodes(1,1) nodes(1,2)],:));
     cnt += 1;
  end
  if (M(nodes(1,2),nodes(1,3)) == 0)
     M(nodes(1,2),nodes(1,3)) = cnt;
     M(nodes(1,3),nodes(1,2)) = cnt;
     p(cnt,:) = mean(p([nodes(1,2) nodes(1,3)],:));
     cnt += 1;
  end
  if (M(nodes(1,1),nodes(1,3)) == 0)
     M(nodes(1,1),nodes(1,3)) = cnt;
     M(nodes(1,3),nodes(1,1)) = cnt;
     p(cnt,:) = mean(p([nodes(1,1) nodes(1,3)],:));
     cnt += 1;
  end
  t(e,4) = M(nodes(1,1),nodes(1,2));
  t(e,5) = M(nodes(1,2),nodes(1,3));
  t(e,6) = M(nodes(1,1),nodes(1,3));
 
end
 
% add boundary edge midpoint node to list of boundary points
for i=1:size(be)
  nodes = be(i,:);
  if (M(nodes(1,1),nodes(1,2)) != 0)
    b = [b; M(nodes(1,1),nodes(1,2))];
  end
end
 
 
% [K,F] = assemble(p,t) % K and F for any mesh of triangles: quadratic phi's
N=size(p,1);T=size(t,1); % number of nodes, number of triangles
% p lists x,y coordinates of N nodes, t lists triangles by 3 node numbers
K=sparse(N,N); % zero matrix in sparse format: zeros(N) would be "dense"
F=zeros(N,1); % load vector F to hold integrals of phi's times load f(x,y)
 
for e=1:T  % integration over one triangular element at a time
  nodes=t(e,:); % row of t = node numbers of the 6 mesh nodes of triangle e
  Pe=[ones(6,1),p(nodes,:),p(nodes,1).^2,p(nodes,1).*p(nodes,2),p(nodes,2).^2]; % 6 by 6 matrix with rows=[1 x y x^2 xy y^2]
  Area=abs(det(Pe(1:3,1:3)))/2; % area of triangle e = half of parallelogram area
  % now compute 6 by 6 Ke and 6 by 1 Fe for element e
  % perform 3 point quadrature over triangle
  lambda = zeros(3,3);
  lambda(1,1) = 2/3;
  lambda(1,2) = 1/6;
  lambda(1,3) = 1/6;
  lambda(2,1) = 1/6;
  lambda(2,2) = 2/3;
  lambda(2,3) = 1/6;
  lambda(3,1) = 1/6;
  lambda(3,2) = 1/6;
  lambda(3,3) = 2/3;
  qp = lambda*p(nodes(1,[1:3]),:);      % three gaussian quadrature points inside triangle element
 
  kvals = feval(k,qp(:,1),qp(:,2));
 
  C=[ones(3,1),qp(:,:),qp(:,1).^2,qp(:,1).*qp(:,2),qp(:,2).^2]; % 3 by 6 matrix with rows=[1 x y x^2 xy y^2] of 3 quadrature nodes on triangle
  V=(Pe'\C')';
  Cx = [zeros(3,1),ones(3,1),zeros(3,1),qp(:,1)*2,qp(:,2),zeros(3,1)];
  Cy = [zeros(3,1),zeros(3,1),ones(3,1),zeros(3,1),qp(:,1),qp(:,2)*2];
  gradV = (Pe'\[Cx' Cy'])';
  w = ones(3,1)/3;
  wk = w.*kvals;
  Ke = zeros(size(gradV,2));
  qpsize=size(qp,1);
  for i=1:qpsize
    Ke = Ke + wk(i)*(gradV(i,:)'*gradV(i,:));                   % Vx
    Ke = Ke + wk(i)*(gradV(i+qpsize,:)'*gradV(i+qpsize,:));     % Vy
  end
  fval = feval(f,qp(:,1),qp(:,2));
  Fe=V'*(w.*fval); % integral of phi over triangle is volume of pyramid: f(x,y)=4
  % multiply Fe by f at quadrature nodes for load f(x,y): three-point quadrature!
  K(nodes,nodes)=K(nodes,nodes)+Ke*Area; % add Ke to 36 entries of global K
  F(nodes)=F(nodes)+Fe*Area; % add Fe to 3 components of load vector F
end   % all T element matrices and vectors now assembled into K and F
 
 
% [Kb,Fb] = dirichlet(K,F,b) % assembled K was singular! K*ones(N,1)=0
% Implement Dirichlet boundary conditions U(b)=0 at nodes in list b
K(b,:)=0; [theta,rho] = cart2pol(p(b,:)); F(b)=sin(theta.*3); % set rows in K and F for boundary nodes.
K(b,b)=speye(length(b),length(b)); % put I into boundary submatrix of K
Kb=K; Fb=F; % Stiffness matrix Kb (sparse format) and load vector Fb
 
% Solving for the vector U will produce U(b)=0 at boundary nodes
U=Kb\Fb;  % The FEM approximation is U_1 phi_1 + ... + U_N phi_N
 
% Plot the FEM approximation U(x,y) with values U_1 to U_N at the nodes
trisurf(t(:,1:3),p(:,1),p(:,2),0*p(:,1),U,'edgecolor','k','facecolor','interp');
view(2),axis([-1 1 -1 1]),axis equal,colorbar
 
% Calculate Energy Norms
enorm = 0;
for e=1:T  % integration over one triangular element at a time
  nodes=t(e,:); % row of t = node numbers of the 6 mesh nodes of triangle e
  Pe=[ones(6,1),p(nodes,:),p(nodes,1).^2,p(nodes,1).*p(nodes,2),p(nodes,2).^2]; % 6 by 6 matrix with rows=[1 x y x^2 xy y^2]
  Area=abs(det(Pe(1:3,1:3)))/2; % area of triangle e = half of parallelogram area
  % perform 3 point quadrature over triangle
  lambda = zeros(3,3);
  lambda(1,1) = 2/3;
  lambda(1,2) = 1/6;
  lambda(1,3) = 1/6;
  lambda(2,1) = 1/6;
  lambda(2,2) = 2/3;
  lambda(2,3) = 1/6;
  lambda(3,1) = 1/6;
  lambda(3,2) = 1/6;
  lambda(3,3) = 2/3;
  qpts = lambda*p(nodes(1,[1:3]),:);      % three gaussian quadrature points inside triangle element
  w = ones(3,1)/3;                        % three quadrature weights
 
  kvals = feval(k,qpts(:,1),qpts(:,2));
 
  uxvals = feval(ux,qpts(:,1),qpts(:,2));
  uyvals = feval(uy,qpts(:,1),qpts(:,2));
 
  enorm = enorm + Area*((w.*kvals)'*(uxvals.^2 + uyvals.^2));
 
end   % all T element matrices and vectors now assembled into K and F
 
nrm = sqrt(enorm);
 
 
enorm = 0;
for e=1:T  % integrate one triangle at a time
  nodes=t(e,:); % row of t = node numbers of the 6 mesh nodes of triangle e
  Pe=[ones(6,1),p(nodes,:),p(nodes,1).^2,p(nodes,1).*p(nodes,2),p(nodes,2).^2]; % 6 by 6 matrix with rows=[1 x y x^2 xy y^2]
  Area=abs(det(Pe(1:3,1:3)))/2; % area of triangle e = half of parallelogram area
 
  % perform 3 point quadrature over triangle
  lambda = zeros(3,3);
  lambda(1,1) = 2/3;
  lambda(1,2) = 1/6;
  lambda(1,3) = 1/6;
  lambda(2,1) = 1/6;
  lambda(2,2) = 2/3;
  lambda(2,3) = 1/6;
  lambda(3,1) = 1/6;
  lambda(3,2) = 1/6;
  lambda(3,3) = 2/3;
  qpts = lambda*p(nodes(1,[1:3]),:);      % three gaussian quadrature points inside triangle element
  w = ones(3,1)/3;
 
  kvals = feval(k,qpts(:,1),qpts(:,2));
 
%  C=[ones(3,1),qpts(:,:),qpts(:,1).^2,qpts(:,1).*qpts(:,2),qpts(:,2).^2]; % 3 by 6 matrix with rows=[1 x y x^2 xy y^2] of 3 quadrature nodes in triangle
  Cx = [zeros(3,1),ones(3,1),zeros(3,1),qpts(:,1)*2,qpts(:,2),zeros(3,1)];
  Cy = [zeros(3,1),zeros(3,1),ones(3,1),zeros(3,1),qpts(:,1),qpts(:,2)*2];
 
  Vx = (Pe'\Cx')'*U(nodes);
  Vy = (Pe'\Cy')'*U(nodes);
 
  uxvals = feval(ux,qpts(:,1),qpts(:,2));
  uyvals = feval(uy,qpts(:,1),qpts(:,2));
 
  enorm = enorm + Area*((w.*kvals)'*((uxvals - Vx).^2 + (uyvals - Vy).^2));
 
end   % all T element matrices and vectors now assembled into K and F
 
err = sqrt(enorm);
 
disp(['Energy Norm: ',num2str(err)]);
disp(['Energy Norm of exact solution: ',num2str(nrm)]);
disp(['Relative energy norm error in computed solution: ',num2str(err/nrm)])
 
% Calculate L2 Norm, which is the displacement error
u2 = feval(u,p(1:Nnorm,1),p(1:Nnorm,2));
u3=U(1:Nnorm,:) - u2(1:Nnorm,:);
L2exact = norm(u2,2);
L2 = norm(u3,2);
 
disp(['L2 Norm: ',num2str(L2)]);
disp(['L2 Norm of exact solution: ',num2str(L2exact)]);
disp(['Relative L2 norm error in computed solution: ',num2str(L2/L2exact)])