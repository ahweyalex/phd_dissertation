% Alexander Moreno
% 07-19-2018
% ------------------------------------------------------------------------
% Description:
% A)
% 1) Generates points in cartesian plane (x,y,& z) for a 3D linear helical
% spiral or rectangular loop. The 3D linear helical spiral has a variable 
% width (i.e. x and y), height (z), number of turns (N), and number of turns 
% wrapped (Nw).
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Input Parameters:
%   N
%   Description: Number of turns
%   Units: <scalar>
% ------------------------------------------------------------------------
%   Nw
%   Description: Number of wraps (in x & y direction)
%   Units: <scalar>
% ------------------------------------------------------------------------
%   wT
%   Description: wire thickness
%   Units: [m] <scalar>
% ------------------------------------------------------------------------
%   rx
%   Description: distance in x-direction 
%   Units: [m] <scalar>
% ------------------------------------------------------------------------
%   ry
%   Description: distance in y-direction 
%   Units: [m] <scalar>
% ------------------------------------------------------------------------
%   hz
%   Description: starting and ending points (height) 
%   Units: [m] <array>
%   Ex: hz = [start, end];
% ------------------------------------------------------------------------
%   type
%   Description: determines type of wire antenna. 
%   Values: 'c' for circular/ellipsoid 
%           'r' for rectangular 
%   Units: <char>
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%   Ouput Parameters:
%
%   xS
%   Description: B-Fields in x-direction
%   Units:
% ------------------------------------------------------------------------
%   yS
%   Description: B-Fields in y-direction
%   Units:
% ------------------------------------------------------------------------
%   zS
%   Description: B-Fields in z-direction
%   Units:
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [xS,yS,zS] = constrWireAnt(N,Nw,wT,rx,ry,hz,type)
    if(type=='c' || type=='C') % circular/ellipsoid structure
%% Construct Circular Wire Antenna Structure
        % Starting and ending points for theta
        % Determines the number of turns
        % Begins at 0[rad] ends at 2*N [rad]
        deltaS = 200;
        start=-pi; fin=(2*N*pi)-pi;
        theta=linspace(start,fin,deltaS);            % theta angles
        xS0 = rx*sin(theta'); yS0 = -ry*cos(theta'); % set ups loop's x and y points
        zS0 = (linspace(hz(1),hz(end),deltaS))';     % set ups loop's z points
        % iterates through the turns in the width (x&y plane)
        if (Nw>=0)
            for wn=1:Nw
                xSwn = (wn*wT + rx)*sin(theta'); ySwn = -(wn*wT + ry)*cos(theta');
                xS0 = [xS0;xSwn]; yS0 = [yS0;ySwn]; zS0 = [zS0;(linspace(hz(1),hz(end),deltaS))'];
            end
        end
        % Adds "thickness" to the wire
        S0 = [ xS0,      yS0,      zS0+wT/2;
               xS0,      yS0,      zS0-wT/2;
               xS0+wT/2, yS0+wT/2, zS0;
               xS0-wT/2, yS0-wT/2, zS0;];
        xS = S0(:,1); yS = S0(:,2); zS = S0(:,3);
        
    elseif(type=='r' || type=='R') % rectangular structure
%% Construct Rect Wire Antenna Structure   
    deltaS = 200;
    xSwn=[]; ySwn=[]; zSwn=[];
    xS0 = (linspace(0,rx/2,ceil(deltaS/8)))';  % first iteration
    xSN = (linspace(-rx/2,0,ceil(deltaS/8)))'; % last iteration
    
   
    xSA = (linspace(-rx/2,rx/2,ceil(deltaS/4)))';      % in-between
    xSB = ((rx/2)*ones(1,ceil(deltaS/4)))';
    xSC = flipud(xSA);
    xSD = -xSB;
    
    ySA = ((ry/2)*ones(1,ceil(deltaS/4)))';
    ySB = (linspace(ry/2,-ry/2,ceil(deltaS/4)))';
    ySC =  -ySA;
    ySD =  flipud(ySB);
    
    z0 = (ones(1,ceil(deltaS/4)))';
    % iterate number of turns
    for nn=1:N
        if(nn==1) % first iteration
            co   = (1.5*wT*nn);
            zn   = co*z0;
            zz   = (co*linspace(nn,(nn+1),ceil(deltaS/4)))'; 
            xSwn = [xS0; xSB; xSC; xSD; xSA];
            ySwn = [ySA(1:numel(xS0)); ySB; ySC; ySD; ySA];
            zSwn = [z0(1:numel(xS0));   z0;  z0; z0;   zz];
        elseif(nn==N) % last iteration
            co   = (1.5*wT*nn);
            zn   = co*z0;
            xSwn = [xSwn; xSN; xSB; xSC; xSD];
            ySwn = [ySwn; ySA(1:numel(xS0)); ySB; ySC; ySD];
            zSwn = [zSwn; zn(1:numel(xS0));  zn;  zn;  zn];
        else % in-betwen iteration(s) 
            co   = (1.5*wT*nn);
            zn   = co*z0;
            zz   = (co*linspace(nn,(nn+1),ceil(deltaS/4)))'; 
            xSwn = [xSwn; xSA; xSB; xSC; xSD];
            ySwn = [ySwn; ySA; ySB; ySC; ySD];
            zSwn = [zSwn;  zz; zn;  zn;   zn];
        end
    end
    %xS=xSwn; yS=ySwn; zS=zSwn;
 
    xS0=xSwn; yS0=ySwn; zS0=zSwn;
    S0 = [xS0,     yS0,      zS0+wT/2;
         xS0,      yS0,      zS0-wT/2;
         xS0+wT/2, yS0+wT/2, zS0;
         xS0-wT/2, yS0-wT/2, zS0;];
    xS = S0(:,1); yS = S0(:,2); zS = S0(:,3);
    else
%% Construct Arbitrary Wire Antenna Structure
    end

end % end 