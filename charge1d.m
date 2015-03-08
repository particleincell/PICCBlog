% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charge1d.m - 1d charging calculation of static dissipative or dielectric materials
% 
% AUTHOR: Alexander C Barrie
%
% LICENSING: This file is intended for educational use. It is licensed
% under GPL (GNU Public License)
%
% This code calculates the charging of a semiconductive strip of material
% measuring 10cm by 1cm with a 3mm thickness. There is an input flux at one
% end and a ground at the other end. This case simulates carbon infused materials
% common to space applications where a dielectric material would typically be used,
% but some level of charge dissipation is required - ESD materials, teflon tubing, etc.
% 
% USE: [dq pot]=charge1d(perm, cond, V0, flux, time, ploton)
% The code is intended to be used in two ways: Setting an initial potential
% with little or no flux will show you how long it takes for a surface
% charge to dissipate to the ground. Having a high flux with little or no
% initial voltage will show you how high that flux will charge the surface
% to. Of course, you can have a high flux and high initial charge to make
% cool movies as well. 
% 
% INPUTS:
% perm - permittivity of material
% cond - conductivity of material
% V0 - initial voltage to put on end of strip
% flux - input flux (#/sq m) on end of strip
% tol - tolerance of qin to qout, which should be equal at steady state
% ploton - 1=turn plotting on, 0=turn plotting off (plotting as in make movie)
% 
% OUTPUTS:
% dq - charge flow between elements. At steady state this will be the same at the 
%      two ends of the strip and ~ 0 in between
% pot - potential on each element of the strip
% mov - an array of frames that you can play with the movie command. This
%       will just be null if ploton is turned off. 
%
% NOTES:
% - The sim is much faster with plotting turned off
% - Modelling conductive materials is impractical due to small timestamp
%   so keep conductivity below ~ 1e-7
% - The number of cells can be changed near the top of the code if need be
%
% EXAMPLE: 
% [dq pot frames]=charge1d(1e-9, 1e-9, 0, 1e14, 1e-3, 0);
% Converged to 0.001000 tolerance with a leading edge potential of 1067.314345V after 3284.356691s and after 18814 iterations 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dq pot frames]=charge1d(perm, cond, V0, flux, tol, ploton)

cells=50;

pot=zeros(1,cells+1);                   %set initial potentials to 0

dx=0.1 / cells;                         %square dimensions of each cell (m)
dy=0.01;
dz=0.003;                               %cell thickness (m)

C=perm*dx*dy/dz;                         %parallel plate based capacitance of each cell

pot(1:round(cells/10))=V0;

eq=1.602e-19;                           %electron charge

%zero out some initial stuff
t=0;
i=zeros(1,length(pot));
dq=zeros(1,length(pot));
dt=zeros(1,length(pot));
dt(end)=1e5;

%set up figure
if ploton==1
    frm=0;
    hFig = figure();
    set(hFig, 'Position', [0 500 1000 100]);

    axis tight
    set(hFig,'nextplot','replacechildren');  
end

x=0:dx:cells*dx;
x=x*100;                                %cm

plotcounter=0;

it=0;

converged=0;

% main loop, finish when specified time is reached
while converged==0
    it=it+1;
    i(:)=0;
    q=C*pot;                            %get initial charges
    % loop through each element
    for cnt=1:length(pot)-1             
        E=-(pot(cnt)-pot(cnt+1))/dx;    %E field
        j=E*cond;                       %current density
        
        % now get the current - negative applied to this cell, positive
        % applied to next cell
        i(cnt)=i(cnt)+j*dy*dz;          
        i(cnt+1)=i(cnt+1)-j*dy*dz;      
        
        % calculate the timestep based on the amount of charge we want to
        % move
        if i(cnt) < 0
            % we only want to move a small amount of charge
            % so we check for 1% movement. If you make this
            % too big, you can overshoot.
            dt(cnt)=abs(0.01*q(cnt)/i(cnt));
        else
            %moving charge TO an element, we don't care so much (set it
            %high)
            dt(cnt)=1;
        end
        
    end
    
    % for the overall timestep, use the smallest value
    dtmin=min(dt);
    
    % now calculate the new charge and potential values
    dq=i*dtmin;
    dqin=flux*dx*dy*dtmin*eq;           %q in due to flux
    dq(1)=dq(1)+dqin;
    dpot=dq/C;
    pot=pot+dpot;
    pot(end)=0;
   
    t=t+dtmin;
    
    %calculate convergence based on current in vs out
    conv=abs(dqin/dq(end) -1 );
    
    % make the new movie frame every 100 timesteps
    if mod(plotcounter,25)==0
        if ploton==1
            frm=frm+1;
            contourf(x,[0 1],[pot;pot], 10);
            set(gca, 'ytick',[]);
            set(gca, 'clim', [0 1000]);
            lbl=sprintf('X (cm), t=%1.6fs',t);
            xlabel(lbl);
            h=colorbar;
            colormap(jet(11));
            set(get(h,'ylabel'),'String', 'Potential (V)');
            frames(frm)=getframe(hFig);
            
        end
        fprintf('time: %fs, dt: %f, conv: %f\n', t, dtmin, conv );
        
    end
    plotcounter=plotcounter+1;

    % check for steady state based on the input and output current
    if conv < tol
        converged=1;
    end
    
    if ploton==0
        frames=0;
    end
    
end

fprintf('Converged to %f tolerance with a leading edge potential of %fV after %fs and after %i iterations\n', conv, pot(1), t, it);


%now get analytic answer to compare
R=0.1/(cond*dz*dy);
I=flux*dx*dy*eq;
V=I*R;

fprintf('Analytical solution yields potential of %fV\n', V);


