clear all; close all; clc;
N=1; Nw=0; wT=1; rx=5;ry=5;hz=[0,0.5]; type='c';

[xS,yS,zS] = constrWireAnt(N,Nw,wT,rx,ry,hz,type);
%%
figure(1)
h=plot3(xS,yS,zS,'o');
xlabel('x'); ylabel('y'); zlabel('z');
grid on; axis equal; %axis tight;
%%
view(90,0)
%%
%{
figure(2)
h=plot3(xS(1),yS(1),zS(1),'o',xS(201),yS(201),zS(201),'o',xS(401),yS(401),zS(401),'o',xS(601),yS(601),zS(601),'o',xS(1),5,0,'d');
xlabel('x'); ylabel('y'); zlabel('z');
zlim([-0.6 0.6]); ylim([4.4 5.6])
grid on; %axis equal; %axis tight;
title('YZ Cut-Plane: Single Segment');
get(h,'MarkerSize',3);
view(90,0)
%}