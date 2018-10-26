% Alexander Moreno
% 07-17-2018
% ------------------------------------------------------------------------
% Description:
% A)
% 1) Computes B-Fields for a given current path. 
% 2) Calucates B-Fields from antenna structure
% 3) Saves the Bx,By,Bz corresponding to each cartesian points (x,y,z)
% ------------------------------------------------------------------------
% [[[Input Parameters]]]
%   I
%   Description: Current through the wire 
%   Units: [A] <scalar>
% ------------------------------------------------------------------------
%   xS
%   Description: Points in the zx-plane where the structure
%   Units: [m] <1D array>
% ------------------------------------------------------------------------
%   yS
%   Description: Points in the y-plane where the structure
%   Units: [m] <1D array>
% ------------------------------------------------------------------------
%   zS
%   Description: Points in the z-plane where the structure
%   Units: [m] <1D array>
% ------------------------------------------------------------------------
%   bBox
%   Description: Boundary box, the space in which the B-Fields will be
%   computed. The boundary box's values are the distance away from the 
%   furthest point of the arbitrary antenna wire structure.
%   Units: [m] <1x3 vector>
%   Example: 
%   bBox = [xdist, ydist, zdist]
% ------------------------------------------------------------------------
%   Ns
%   Description: Number of segments in space in which the B-Fields will be
%   computed
%   Units: [none] <scalar>
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%   [[[Ouput Parameters]]]
%
%   Bx
%   Description: B-Fields in x-direction
%   Units:
% ------------------------------------------------------------------------

function [X,Y,Z,BX,BY,BZ,normB] = CalcBFields_Wire_Antenna(I,xS,yS,zS,bBox,Ns)
    %% Initialize variables
    % Ns = 20;         % Total segments in space <scalar> 
    mu0 = 4*pi*1e-7; % free space permeability <scalar> [H/m]
    dS = 1e9;        % filament max discretization step [m]
    xdis = bBox(1);
    ydis = bBox(2);
    zdis = bBox(3);
    % Calc Boundaries for cartesian points of interest
    xi=ceil(min(xS)); xf=ceil(max(xS)); Nx=Ns; % xdis=abs(xi-xf);
    yi=ceil(min(yS)); yf=ceil(max(yS)); Ny=Ns; % ydis=abs(yi-yf);
    zi=ceil(min(zS)); zf=ceil(max(zS)); Nz=Ns; % zdis=abs(zi-zf);
    % discrete points in space 
    x_M = linspace( xi-xdis/4, xf+xdis/4, Nx);
    y_M = linspace( yi-ydis/4, yf+ydis/4, Ny);
    z_M = linspace( zi-zdis/4, zf+zdis/4, Nz);
    [X,Y,Z]=meshgrid(x_M,y_M,z_M);
    % Initialize B-Fields matrices
    BX = zeros(Ny,Nx,Nz); 
    BY = zeros(Ny,Nx,Nz); 
    BZ = zeros(Ny,Nx,Nz);
    
    %% Discretization of source
    xP=[]; yP=[]; zP=[];
    Ns = numel(xS)-1;
    % Source of current flow along the wire antenna
    S = [xS;yS;zS]; 
    for sn=1:Ns
        L_Si = norm(S(sn,:)-S(sn+1,:));
        NP  = ceil(L_Si/dS);
        
        xP = [xP, linspace(xS(sn),xS(sn+1), NP)];
        yP = [yP, linspace(yS(sn),yS(sn+1), NP)];
        zP = [zP, linspace(zS(sn),zS(sn+1), NP)];
    end
    %% Compute B-Fields
    for yn=1:size(X,1)          % iterate y-points (points of interest)
        for xn=1:size(X,2)      % iterate x-points (points of interest)
            for zn=1:size(X,3)  % iterate z-points (points of interest)
                % M is point of interest for Bfield
                xM = X(yn,xn,zn);
                yM = Y(yn,xn,zn);
                zM = Z(yn,xn,zn);
                for n=1:length(xP)-1    % iterate through Source (points)
                    R = (sqrt((xM-xP(n))^2 + (yM-yP(n))^2 + (zM-zP(n))^2))^3; % source to point of interest
                    % cross product(s)
                    % reference: Biot-Savart Law
                    % (dl x R)/ R^3
                    dBx(n) = ((yP(n+1)-yP(n))*(zM-zP(n)) - (zP(n+1)-zP(n))*(yM-yP(n)))/R; % cross product x
                    dBy(n) = ((zP(n+1)-zP(n))*(xM-xP(n)) - (xP(n+1)-xP(n))*(zM-zP(n)))/R; % cross product y
                    dBz(n) = ((xP(n+1)-xP(n))*(yM-yP(n)) - (yP(n+1)-yP(n))*(xM-xP(n)))/R; % cross product z
                end
                % B-Fields at that point in space, summed all source point calcs 
                % reference: Biot-Savart Law
                % (I*mu0/4pi)*((dl x R)/ R^3) [T]
                BX(yn,xn,zn) = BX(yn,xn,zn) + mu0*I/4/pi*sum(dBx);
                BY(yn,xn,zn) = BY(yn,xn,zn) + mu0*I/4/pi*sum(dBy);
                BZ(yn,xn,zn) = BZ(yn,xn,zn) + mu0*I/4/pi*sum(dBz);
            end
        end
    end
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    nBX = BX./normB;
    nBY = BY./normB;
    nBZ = BZ./normB;
end