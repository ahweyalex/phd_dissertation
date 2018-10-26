clear all; close all; clc;

h=10; ra=10; ri=10; phi=90; N=4; O=1; type='c'; wT=1;%8.2515e-3;
[xS,yS,zS] = constrWireAnt_10_25_2018(h,ra,ri,phi,N,O,type,wT/2);

figure(1)
h=plot3(xS,yS,zS);
xlabel('x'); ylabel('y'); zlabel('z');
grid on; axis equal; %axis tight;
%view(0,90)