%---------------------------------------------------
%  NAME:      example_3D_solenoid_filament.m
%  WHAT:      Calculation of the magnetic field of a solenoid
%             on a volume (+ 3D plot).
%  REQUIRED:  BSmag Toolbox 20150407
%  AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%  COPYRIGHT: 2015, Loic Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause).
%----------------------------------------------------

%% Initialize
clear all; close all; clc;
BSmag = BSmag_init(); % Initialize BSmag analysis

%% Source points (where there is a current source)
%theta = linspace(-2*2*pi,2*2*pi,2*100);
%theta = linspace(-(1/2)*pi,(3/2)*pi,200); 
%z0 = zeros(200,1);
%Gamma = [cos(theta'),sin(theta'),z0]; % x,y,z [m,m,m]
%{
Gamma   = [ 1,-1,0; % x,y,z [m,m,m]
            1,1,0;
           -1,1,0;
           -1,-1,0;
            1,-1,0];
%}
N=200;
A11 = ones(floor(N/4),1);
A1n = linspace(1,-1,floor(N/4))';
Ann = -1*ones(floor(N/4),1);
An1 = linspace(-1,1,floor(N/4))';
z0 = zeros(floor(N/4),1);
Gamma = [A11,An1,z0;
          A1n,A11,z0;
          Ann,A1n,z0;
          An1,Ann,z0];
I = 1; % filament current [A]
dGamma = 1e9; % filament max discretization step [m]
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);
%% Field points (where we want to calculate the field)
x_M = linspace(-1.5,1.5,21); % x [m]
y_M = linspace(-2,2,22); % y [m]
z_M = linspace(-2,2,23); % z [m]
[X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points volume

%% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);
%%
normB=sqrt(BX.^2+BY.^2+BZ.^2);
BXn=BX./normB; BYn=BY./normB; BZn=BZ./normB;
%% Plot B/|B|
figure(1)
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    quiver3(X,Y,Z,BXn,BYn,BZn,'b');
%axis tight
%%
figure(9)
quiver(Y,X,BYn,BXn);
%% Plot Bz on the volume
figure(2), hold on, box on, grid on
	plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
	slice(X,Y,Z,BZ,[0],[],[-1,0,1]), colorbar % plot Bz
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
view(3), axis equal, axis tight
caxis([-0.5,0.5]*1e-5)

%% Plot some flux tubes
figure(3), hold on, box on, grid on
	plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
	[X0,Y0,Z0] = ndgrid(-1.5:0.5:1.5,-1.5:0.5:1.5,-2); % define tubes starting point        
	htubes = streamtube(stream3(X,Y,Z,BX,BY,BZ,X0,Y0,Z0), [0.2 10]);
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Some flux tubes')
view(3), axis equal, axis tight
set(htubes,'EdgeColor','none','FaceColor','c') % change tube color
camlight left % change tube light