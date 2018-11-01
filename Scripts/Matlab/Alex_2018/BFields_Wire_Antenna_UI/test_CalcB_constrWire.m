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
ra=10; ri=10; phi=90; N=4; O=1; type='c'; wT=1; h=(1.1)*(2*wT*N);
[xS,yS,zS] = constrWireAnt_10_25_2018(h,ra,ri,phi,N,O,type,wT/2);
figure(1)
h=plot3(xS,yS,zS);
xlabel('x'); ylabel('y'); zlabel('z');
grid on; axis equal; %axis tight;
%view(0,90)
%% testing CalcBFields_Wire_Antenna

% [X,Y,Z,BX,BY,BZ,normB] = CalcBFields_Wire_Antenna(I,xS,yS,zS,bBox,Ns)