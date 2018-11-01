clear all; close all; clc;

% Number of grids in Z and Y-axis
Nz=51; Ny=51;
% Number of grids in the coil (X-Y plane)
N=25;
% radius of coil in X-Y plane
Ra=6;
% Current in the coil
I=3;
% permitivity
u0=1;
phi = -(1/2)*pi:(2/(N-1))*pi:(3/2)*pi;
% coordinates of coil
Xc = Ra*cos(phi);
Yc = Ra*sin(phi);
% coordinates of the plane
yp = -25:1:25;
zp = 0:1:50;
% Array for 1D to 2D conversion of coordinates
Y(1:Ny,1:Nz)=0; 
Z(1:Ny,1:Nz)=0;
for i=1:Ny
    Y(i,:)=yp(i); % all y-coordinates value in 2-d form
end
for i=1:Nz
    Z(:,i)=zp(i);% all z-coordinates value in 2-d form
end

%% Calc stuff
for yn=1:Ny
    for zn=1:Nz
       %------------------------------------------------------------------
       % Calc R-vector from the coil (XY plane) to point of interest (YZ Plane)
       % and dl-vector along the coil where current is flowing
       % NOTE:
       % R is source (XY plane) to point of interest (YZ plane)
       % dl is current element vector which makes up the coil
       %------------------------------------------------------------------
       for i=1:N-1
           % [Alex] not understanding why Rx,Ry,Rz are being computed in
           % such a way
           Rx(i)  = -0.5*(Xc(i)+Xc(i+1));
           Ry(i)  = (Y(yn,zn)) - (0.5*(Yc(i) + Yc(i+1)));
           Rz(i)  = Z(yn,zn);
           dlx(i) = Xc(i+1) - Xc(i);
           dly(i) = Yc(i+1) - Yc(i);
       end
       Rx(N)  = -0.5*(Xc(N)+Xc(1));
       Ry(N)  = Y(yn,zn) - (0.5*(Yc(N) + Yc(1)));
       Rz(N)  = Z(yn,zn);
       dlx(N) = -Xc(N) + Xc(1);
       dly(N) = -Yc(N) + Yc(1);
       %%
       %------------------------------------------------------------------
       % for all elements along coil, calc dl cross R
       % dl cross R is the curl of vector dl and R
       % XCross is X-component of the curl of dl and R, so on for Y and Z
       %------------------------------------------------------------------
       for i=1:N
            Xcross(i) =  dly(i).*Rz(i);
            Ycross(i) = -dlx(i).*Rz(i);
            Zcross(i) = (dlx(i).*Ry(i)) - (dly(i).*Rx(i));
            R(i) = sqrt(Rx(i).^2 + Ry(i).^2 + Rz(i).^2);
       end
       %%
       %------------------------------------------------------------------
       % Biot Savart Law Equ
       %------------------------------------------------------------------
       Bx1 = (I*u0)./(4*pi*(R.^3)).*Xcross;
       By1 = (I*u0)./(4*pi*(R.^3)).*Ycross;
       Bz1 = (I*u0)./(4*pi*(R.^3)).*Zcross;
       %%
       %------------------------------------------------------------------
       % sum magnetic field to get total magnetic field
       %------------------------------------------------------------------       
       BX(yn,zn)=0; BY(yn,zn)=0; BZ(yn,zn)=0; % initialize sum mag field
       % loop over all currents elements along coil
       for i=1:N
        BX(yn,zn) = BX(yn,zn)+Bx1(i);
        BY(yn,zn) = BY(yn,zn)+By1(i);
        BZ(yn,zn) = BZ(yn,zn)+Bz1(i);
       end
    end
end

%% Figures