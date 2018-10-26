clear all; close all; clc;
N  = 200;
Nf = floor(N/4);
W = 4;
L = 1;
xn1=-1; xp1=1;
yn1=-1; yp1=1;
z0 = zeros(Nf,1);
A2Bx = W.*xn1.*ones(Nf,1);          A2By = L.*linspace(yp1,yn1,Nf)';
B2Cx = W.*linspace(xn1,xp1,Nf)';    B2Cy = L.*yn1*ones(Nf,1);
C2Dx = W.*xp1*ones(Nf,1);           C2Dy = L.*linspace(yn1,yp1,Nf)';
D2Ax = W.*linspace(xp1,xn1,Nf)';    D2Ay = L.*yp1.*ones(Nf,1);
Source = [ A2Bx, A2By, z0;
           B2Cx, B2Cy, z0;
           C2Dx, C2Dy, z0;
           D2Ax, D2Ay, z0;];
%A11 = ones(floor(N/4),1);
%A1n = linspace(1,-1,floor(N/4))';
%Ann = -1*ones(floor(N/4),1);
%An1 = linspace(-1,1,floor(N/4))';
%z0 = zeros(floor(N/4),1);
%cI=2;
%Source = [A11,An1,z0;
%          A1n,A11,z0;
%          Ann,A1n,z0;
%          An1,Ann,z0];
      
%% Plot Antenna (Structure)
figure(1)
    H=plot3(Source(:,1),Source(:,2),Source(:,3),'.-r');
    xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
    title('Antenna Structure');
    axis tight; grid on; box on; axis equal;