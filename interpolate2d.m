% interpolate2d.m
%
% illustrates scatter / gather operation using a polygonal element
% See http://www.particleincell.com/2012/quad-interpolation for additional information
% (c) 2012 Lubos Brieda, lubos.brieda@particleincell.com
%
% -- tested with Octave 3.2.4 on Windows 

%main function
function []=interpolate2d()

	%create our polygon
	px = [-1, 8, 13, -4];
	py = [-1, 3, 11, 8];

	%compute coefficients
	A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0];
	AI = inv(A);
	a = AI*px';
	b = AI*py';

	%plot random internal points
	plot_points_in_poly(px,py,a,b);

	%classify points as internal or external
	plot_internal_and_external_points(px,py,a,b);

	%scatter data 
	scatter_data(px,py,a,b);
	
	%gather_data
	gather_data(px,py,a,b);
	
end


%plots random points inside the polygon
function []=plot_points_in_poly(px,py,a,b)

	%pick random logical coordinate and convert to physical
	for i=1:100
		l=rand();
		m=rand();
		x(i) = a(1) + a(2)*l + a(3)*m + a(4)*l*m;
		y(i) = b(1) + b(2)*l + b(3)*m + b(4)*l*m;	
	end

	%plot
	figure(1);
	hold off;
	plot([px px(1)],[py py(1)]);
	hold on;
	plot (x,y,'x');
end

%picks random points and colors them based on internal or external
function []=plot_internal_and_external_points(px,py,a,b)
	
	%bounding box
	x0 = min(px);
	lx = max(px)-min(px);
	y0 = min(py);
	ly = max(py)-min(py);
	
	%convert random physical coordinates to logical
	for i=1:200
	
		%pick a random point in the bounding box
		x(i) = x0+rand()*lx;
		y(i) = y0+rand()*ly;
		
		%calculate logical coordinates
		[l,m] = XtoL(x(i),y(i),a,b);		
	
		%is the point outside the quad?
		if (m<0 || m>1 || l<0 || l>1)
			type(i)=1;
		else
			type(i)=0;
		end
	end

	%plot
	figure(2);
	hold off;
	plot([px px(1)],[py py(1)]);
	hold on;
	scatter(x,y,4,type,'x');
end

%scatters counts to the four grid nodes
function []=scatter_data(px,py,a,b)
	
	%init
	c = [0 0 0 0];
	
	%pick random logical coordinate and convert to physical
	for i=1:100
		
		%random point inside the quad
		l=rand();
		m=rand();
		
		%deposit counts, note here l,m=[0:1]. In a real 2d mesh
		%the index we would first need to subtract the integral
		%cell index to get the cell-local distance, dl = l - (int)l;
		dl = l;
		dm = m;
		
		c(1) += (1-dl)*(1-dm);
		c(2) += dl *(1-dm);
		c(3) += dl * dm;
		c(4) += (1-dl)*dm;
	end

	%scale data, uniform distribution = 10
	c = 40*(c/i);
	
	%output
	disp(c);
	
	%plot
	figure(3);
	hold off;
	plot([px px(1)],[py py(1)]);
	hold on;
	scatter(px,py,c,[1 0 0], 'o');
end

%interpolates data from the polygon vertices onto a grid
function []=gather_data(px,py,a,b)
	
	%bounding box
	x0 = min(px);
	lx = max(px)-min(px);
	y0 = min(py);
	ly = max(py)-min(py);
	
	%set spacing such that we have 30x30 grid (i.e. 29 cells)
	ni = 30;
	nj = 15;
	
	dx = lx/(ni-1);
	dy = ly/(nj-1);
	
	%node positions (for plotting)
	for i=1:ni
		x(i) = x0 + (i-1)*dx;
	end
	for j=1:nj
		y(j) = y0 + (j-1)*dy;
	end
	
	%node values
	c = [1 2 3 4];
	
	val=zeros(nj,ni);
	
	for i=1:ni
		for j=1:nj
		
			%convert to logical coordinates
			[l m] = XtoL(x(i),y(j),a,b);	
			
			%evaluate only if we are inside
			if (l>0 && l<=1 && m>=0 && m<=1)
				
				%again, if we had multiple cells, we would set dl = l - (int)l
				dl = l;
				dm = m;
				
				val(j,i) = (1-dl)*(1-dm)*c(1) + ...
						   dl*(1-dm)*c(2) + ...
						   dl*dm*c(3) +...
						   (1-dl)*dm*c(4);		
			end
		end
	end
	
	%plot
	figure(4);
	hold off;
	contourf(x,y,val);
	hold on;
	plot([px px(1)],[py py(1)],'w','LineWidth',4);
end

% converts physical (x,y) to logical (l,m)
function [l, m] = XtoL(x,y,a,b)
	%quadratic equation coeffs, aa*mm^2+bb*m+cc=0
	aa = a(4)*b(3) - a(3)*b(4);
	bb = a(4)*b(1) -a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + x*b(4) - y*a(4);
	cc = a(2)*b(1) -a(1)*b(2) + x*b(2) - y*a(2);
	
	%compute m = (-b+sqrt(b^2-4ac))/(2a)
	det = sqrt(bb*bb - 4*aa*cc);
	m = (-bb+det)/(2*aa);
	
	%compute l
	l = (x-a(1)-a(3)*m)/(a(2)+a(4)*m);
end
