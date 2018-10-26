clear all; close all; clc;

antCase=2; coil=1;
N  = 1;    % number of turns in z-direction
Nw = 2;   % number of turns (width) in xy-direction
wT = 0.00162814; %[m] is 14 AWG (0.06410in)
I=1;

dSource = 1e9; % filament max discretization step [m]
deltaS=200;
%%
if(antCase==1) % circular loop
    radiusX=1; radiusY=1; heightZ=1;
    % calc number of turns aka N
    if(coil==1) 
        % coiled
        N=4; %number of turns, needs to be calc in way
        z0=ones(deltaS,1);
        start=-(2*N)*pi; fin=(2*N)*pi; theta=linspace(start,fin,deltaS);
        Source = [radiusX*cos(theta'),radiusY*sin(theta'),heightZ*(theta/10)']; % x,y,z [m,m,m]
    else 
        % no coil
        z0=zeros(deltaS,1);
        start=-(1/2)*pi; fin=(3/2)*pi; theta=linspace(start,fin,deltaS);
        Source = [radiusX*cos(theta'),radiusY*sin(theta'),(heightZ*z0)]; % x,y,z [m,m,m]
    end
elseif(antCase==2) % rectangle loop
    W=4; L=1; Nf=floor(deltaS/4);
    if(coil==1) % coiled
        
        N=1; %number of turns,  needs to be calc in way
        xn1=-1/2; xp1=1/2;
        yn1=-1/2; yp1=1/2;
        z0 = zeros(Nf,1);
        A2Bx = W.*xn1.*ones(Nf,1);          A2By = L.*linspace(yp1,yn1,Nf)';
        B2Cx = W.*linspace(xn1,xp1,Nf)';    B2Cy = L.*yn1*ones(Nf,1);
        C2Dx = W.*xp1*ones(Nf,1);           C2Dy = L.*linspace(yn1,yp1,Nf)';
        D2Ax = W.*linspace(xp1,xn1,Nf)';    D2Ay = L.*yp1.*ones(Nf,1);
        Source = [];
        climb = linspace(0,1,Nf)';
        
        N=5; heightZ=1;
       for nn=0:N
           if nn==0
               Source = [ Source;
                          D2Ax(end:-1:floor(numel(D2Ax)/2)), D2Ay(end:-1:floor(numel(D2Ax)/2)), z0(end:-1:floor(numel(D2Ax)/2))+nn*heightZ;
                          A2Bx, A2By, z0+nn*heightZ;
                          B2Cx, B2Cy, z0+nn*heightZ;
                          C2Dx, C2Dy, z0+nn*heightZ;
                          D2Ax, D2Ay, z0+heightZ*(climb+nn);];
           elseif(nn==N)
               Source = [ Source;
                          D2Ax(1:floor(numel(D2Ax)/2)), D2Ay(1:floor(numel(D2Ax)/2)), z0(1:floor(numel(D2Ax)/2))+(N-1)*heightZ];
           elseif(nn==N-1)
               Source = [ Source;
                          A2Bx, A2By, z0+nn*heightZ;
                          B2Cx, B2Cy, z0+nn*heightZ;
                          C2Dx, C2Dy, z0+nn*heightZ;];
                          %D2Ax, D2Ay, z0+heightZ*(climb+nn);];
           else
               Source = [ Source;
                          A2Bx, A2By, z0+nn*heightZ;
                          B2Cx, B2Cy, z0+nn*heightZ;
                          C2Dx, C2Dy, z0+nn*heightZ;
                          D2Ax, D2Ay, z0+heightZ*(climb+nn);];
           end
       end
        
    else 
        % no coil
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
    end
elseif(antCase==3)  
end
%}

%%
%theta = linspace(-(1/2)*pi,(3/2)*pi,200);
%z0 = zeros(numel(theta),1);
%Source = [cos(theta'),sin(theta'),z0]; % x,y,z [m,m,m]

%% Plot Antenna (Structure)
x = Source(:,1); y = Source(:,2); z = Source(:,3);
figure(1)
%H=plot3(Source(:,1),Source(:,2),Source(:,3),'.-r');
H=plot3(x,y,z);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
title('Antenna Structure');
    grid on; box on; axis tight; axis equal;
%% Point of interets
xi=-1.5; xf=1.5; Nx=21;
yi=-2;   yf=2;   Ny=22;
zi=-2;   zf=2;   Nz=23;

x_M = linspace(xi,xf,Nx);    % x [m]
y_M = linspace(yi,yf,Ny);    % y [m]
z_M = linspace(zi,zf,Nz);    % z [m]
[X,Y,Z]=meshgrid(x_M,y_M,z_M);

%% Calc B-Fields
mu0 = 4*pi*1e-7;
BX = zeros(Ny,Nx,Nz); 
BY = zeros(Ny,Nx,Nz); 
BZ = zeros(Ny,Nx,Nz);
%% Discretization of source
xP=[]; yP=[]; zP=[];
% Total source points
Ns = size(Source,1)-1;
for sn=1:Ns
    L_Si = norm(Source(sn,:)-Source(sn+1,:));
    NP  = ceil(L_Si/dSource);
    xP = [xP, linspace(Source(sn,1),Source(sn+1,1),NP)];
    yP = [yP, linspace(Source(sn,2),Source(sn+1,2),NP)];
    zP = [zP, linspace(Source(sn,3),Source(sn+1,3),NP)];
end
%%

for yn=1:size(X,1)%Nx             % iterate x-points 
    for xn=1:size(X,2)%Ny         % iterate y-points 
        for zn=1:size(X,3)%Nz     % iterate z-points     
            % M is point of interest for Bfield
            xM = X(yn,xn,zn);
            yM = Y(yn,xn,zn);
            zM = Z(yn,xn,zn);
            for n=1:length(xP)-1    % iterate through Source
                % Source to Point of Interest
                %R = ((sqrt(xM-xP(nn)))^2 + (sqrt(yM-yP(nn)))^2 + (sqrt(zM-zP(nn)))^2 )^3;
                %dBx(nn) = ((yP(nn+1)-yP(nn))*(zM-zP(nn)) - (zP(nn+1)-zP(nn))*(yM-yP(nn)))/R;
                %dBy(nn) = ((zP(nn+1)-zP(nn))*(xM-xP(nn)) - (xP(nn+1)-xP(nn))*(zM-zP(nn)))/R;
                %dBz(nn) = ((xP(nn+1)-xP(nn))*(yM-yP(nn)) - (yP(nn+1)-yP(nn))*(xM-xP(nn)))/R;
                PkM3 = (sqrt((xM-xP(n))^2 + (yM-yP(n))^2 + (zM-zP(n))^2))^3; % source to point of interest
                dBx(n) = ((yP(n+1)-yP(n))*(zM-zP(n)) - (zP(n+1)-zP(n))*(yM-yP(n)))/PkM3; % cross product x
                dBy(n) = ((zP(n+1)-zP(n))*(xM-xP(n)) - (xP(n+1)-xP(n))*(zM-zP(n)))/PkM3; % cross product y
                dBz(n) = ((xP(n+1)-xP(n))*(yM-yP(n)) - (yP(n+1)-yP(n))*(xM-xP(n)))/PkM3; % cross product z
            end
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
%%
figure(44)
 normB=sqrt(BX.^2+BY.^2+BZ.^2);
quiver(Y,Z,BY./normB,BZ./normB,'b')
%% Plot B/|B|
%figure(2)
%quiver3(X,Y,Z,nBX,nBY,nBZ,'b');