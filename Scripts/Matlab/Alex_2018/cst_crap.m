% Alexander Moreno
% 
clear all; close all; clc;

% fileName = 'bfield_3.txt';
n=500;
nx=n; ny=1; nz=n;
fileName = 'bfp.txt';
CSTB = importdata(fileName);
xm = reshape(CSTB.data(:,1),[nx,ny,nz]);
ym = reshape(CSTB.data(:,2),[nx,ny,nz]);
zm = reshape(CSTB.data(:,3),[nx,ny,nz]);

xR = reshape(CSTB.data(:,4),[nx,ny,nz]);
yR = reshape(CSTB.data(:,5),[nx,ny,nz]);
zR = reshape(CSTB.data(:,6),[nx,ny,nz]);
xI = reshape(CSTB.data(:,7),[nx,ny,nz]);
yI = reshape(CSTB.data(:,8),[nx,ny,nz]);
zI = reshape(CSTB.data(:,9),[nx,ny,nz]);
%%
figure(1)
%normB=sqrt(xR.^2+yR.^2+zR.^2);
%xIn=xR./normB; zIn=zR./normB; yIn=yR./normB;

normB=sqrt(xI.^2+yI.^2+zI.^2);
xIn=xI./normB; zIn=zI./normB; yIn=yI./normB;

%quiver3(xm,ym,zm,xIn,yIn,zIn,'b');
quiver(ym,zm,yIn,zIn,'b');
view(90,0);
%%
figure(2)
xIn=abs(1j.*xI+xR); zIa=abs(1j.*zI+zR); yIa=abs(1j.*yI+yR);
quiver(ym,zm,yIa,zIa,'b');