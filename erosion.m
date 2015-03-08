function [max_erosion_mag max_erosion_deg points] = erosion(time, numpoints, plotting)

%EROSION.m 
%Author: Alex Barrie

%This code performs an erosion calculation on a circle. There is no element splittting
%or merging, so don't erode for too long! A simple cosine model is used for
%angular dependence. 

% INPUTS:
% time:           time to erode in seconds
% numpoints:      number of points to make the circleh
% plotting:       1 or 0 to turn the plotting on or off
% 
% OUTPUTS:
% max_erosion_mag             magnitude of maximum erosion
% max_erosion_deg             location of maximum erosion
% points                      array of points making up eroded circle      
% 

%material params
clc;
normal_yield=1;                 %atoms sputtered / incident ion @ normal inc
mol_wt = 26.98;                 %Aluminum
density = 2700;                 %Aluminum kg/m^3

%sim setup
r=0.005;
maxflux = 1.0e16;               %num particles per sq meter per sec
fluxdir=[0 -1];
dt=time/300;
av_num=6.0221415e23;

close all;

%make surface

rad=ones(numpoints+1,1)'*r;
theta=linspace(pi,-pi,numpoints+1);

%we use cartesian coords since it is easier to think about this way once we
%start moving surfaces around 
[points(:,1) points(:,2)]=pol2cart(theta, rad);
points=points(1:end-1,:);

points_orig=points;

function [point]=P(i)
    
%This function acts like a circular reference - letting us do +1 or -1

    while i > length(points)
        i=i-length(points);
    end
    while i < 1
        i=i+length(points);
    end
    
    point = i;
    
end

function [angle]=fix(a)

%get angle into right range

    while (a < 0)
       a = a + pi/2;
    end
    while (a > pi/2)
       a = a - pi/2;
    end

    angle=a;
end

for i=1:dt:time
    
    if plotting == true
        plot([points(:,1); points(1, 1)], [points(:,2); points(1, 2)],'-');
        axis equal;
        hold on;
    end
    
    for j=1:numpoints
        
        
        %get surface vecs for left and right panels
        vl=points(P(j-1),:)-points(j,:);
        vr=points(P(j),:)-points(P(j+1),:);
        
        %area is just a fraction of unit area visible (ie dot)
        areal=abs(dot(fluxdir,vl));
        arear=abs(dot(fluxdir,vr));
        
        vl=vl/norm(vl);
        vr=vr/norm(vr);
        
        %surface normals of neigboring panels
        nl=cross([0 0 1], [vl 0]);
        nr=cross([0 0 1], [vr 0]);
        
        nl=-nl(1:2);
        nr=-nr(1:2);
        
        %can we go along flux line?
        dl=0;
        dr=0;
        
        if dot(fluxdir, nl) <= 0
           dl=1;
        end
        if dot(fluxdir, nr) <= 0
           dr=1;
        end
        
        %boolean math to choose appropriate erosion direction
        erodedir=(dl*dr*fluxdir) + (1-dl)*vl - (1-dr)*vr;
        erodedir=erodedir/norm(erodedir);
        
        %grab half area of each panel
        anglel=acos(dot(fluxdir,vl/norm(vl)));
        angler=acos(dot(fluxdir,vr/norm(vr)));
        
        anglel=fix(anglel);
        angler=fix(angler);
        
        fluxl=maxflux*dt*areal;
        fluxr=maxflux*dt*arear;

        %this checks if the flux is hitting the underside of the panel
        %instead of the top, ie the bottom half of the circle
        fluxl=(fluxl>0) * fluxl;
        fluxr=(fluxr>0) * fluxr;

        flux=((dl+dr)>0) * (fluxl+fluxr);

        %get sputter yield f(f(0), theta)
        yieldl = normal_yield*cos(anglel-pi/4);
        yieldr = normal_yield*cos(angler-pi/4);   

        yield=yieldl+yieldr;

        %get the volume of material to remove
        erode_vol=(yield*flux*mol_wt)/(av_num*density);

        %now get dx (divide by dy)
        dx=2*erode_vol/(areal+arear);

        hold on;
        axis equal;

        %below are some more plots you may want to turn on (vectors and
        %such)
       
%         plot(points(P(j-1),1), points(P(j-1),2),'ok');
%         plot(points(j,1), points(j,2),'xk');
%         plot(points(P(j+1),1), points(P(j+1),2),'+k');
%         plot([points(j,1), points(j,1) + 3e-3*nl(1)], [points(j,2), points(j,2) + 3e-3*nl(2)], '-r');
%         plot([points(j,1), points(j,1) + 3e-3*nr(1)], [points(j,2), points(j,2) + 3e-3*nr(2)], '-r');
%         plot([points(j,1), points(j,1) + 5e-3*erodedir(1)], [points(j,2), points(j,2) + 5e-3*erodedir(2)], '-k');
%         plot([points(P(j+1),1), points(P(j+1),1) + 5e-3*erode2dir(1)], [points(P(j+1),2), points(P(j+1),2) + 5e-3*erode2dir(2)], '-r');
% 


        points(j,:)=points(j,:) + dx*erodedir;
    end
 
end

%final plots (thick line)
plot([points(:,1); points(1, 1)], [points(:,2); points(1, 2)],'-xg', 'LineWidth', 2);        
plot([points_orig(:,1); points_orig(1, 1)], [points_orig(:,2); points_orig(1, 2)],'-or', 'LineWidth', 2);
        
%now convert back to polar
[theta_eroded, rad_eroded] = cart2pol(points(:,1), points(:,2));

[val, ind]=max(r-rad_eroded);

max_erosion_mag = val;
max_erosion_deg = theta_eroded(ind)*180/pi;
points = [theta_eroded, rad_eroded];

end