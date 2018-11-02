% Alexander Moreno
% 11-01-2018
%
% Description:
% Testing the functions 
% 1) CalcBFields_Wire_Antenna 
% 2) constrWireAnt
%
clear all; close all; clc;

%% testing constrWireAnt_10_25_2018
% h=10; 
ra=1; ri=1; phi=10; N=4; O=1; wT=0.25; h=(1.1)*(2*wT*N);
[xS,yS,zS] = constrWireAnt(h,ra,ri,phi,N,O,wT);
%%
%
figure(1)
h=plot3(xS,yS,zS,'-');
xlabel('x'); ylabel('y'); zlabel('z');
grid on; axis equal; %axis tight;
%view(0,90)
%}

%% testing CalcBFields_Wire_Antenna
I = 1; Ns = 50;
xminb=h; yminb=h; zminb=h;
xmaxb=h; ymaxb=h; zmaxb=h;
bBox = [xminb,yminb,zminb; xmaxb,ymaxb,zmaxb];
[X,Y,Z,BX,BY,BZ,normB] = CalcBFields_Wire_Antenna(I,xS,yS,zS,bBox,Ns);

%% Plot B-Fields
plotBFields(X,Y,Z,BX,BY,BZ)