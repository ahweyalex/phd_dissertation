%% Initialize
%clear all; close all; clc;
BSmag = BSmag_init(); % Initialize BSmag analysis

%% Source points (where there is a current source)
%theta = linspace(-2*2*pi,2*2*pi,200);
theta = linspace(-(1/2)*pi,(3/2)*pi,200);
z0 = zeros(numel(theta),1);
Gamma = [cos(theta'),sin(theta'),z0]; % x,y,z [m,m,m]
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
%function [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X,Y,Z)
%[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);
X=X_M; Y=Y_M; Z=Z_M;

mu0 = 4*pi*1e-7; % vacuum permeability [N/A^2]

%B = zeros( Ny, Nx, Nz)
BX = zeros(size(X,1),size(X,2),size(X,3));
BY = zeros(size(X,1),size(X,2),size(X,3));
BZ = zeros(size(X,1),size(X,2),size(X,3));

Gamma = BSmag.filament.Gamma;
dGamma = BSmag.filament.dGamma;
I = BSmag.filament.I;
%%
% Discretization of Gamma
x_P = []; y_P = []; z_P = [];
N = size(Gamma,1)-1; % Number of points defining Gamma
for i = 1:N % Loop on the segments defining gamma
    L_Gamma_i = norm(Gamma(i,:)-Gamma(i+1,:));
    NP = ceil(L_Gamma_i/dGamma); % Number of points required to have a discretization step smaller than dGamma
    x_P = [x_P,linspace(Gamma(i,1), Gamma(i+1,1), NP)]; % discretization of Gamma for x component
    y_P = [y_P,linspace(Gamma(i,2), Gamma(i+1,2), NP)]; % discretization of Gamma for y component
    z_P = [z_P,linspace(Gamma(i,3), Gamma(i+1,3), NP)]; % discretization of Gamma for z component
end

%%
% Add contribution of each source point P on each field point M (where we want to calculate the field)
for m = 1:size(X,1)           % y
    for n = 1:size(X,2)       % x
        for p = 1:size(X,3)    % z

        % M is the field point
        x_M = X(m,n,p);
        y_M = Y(m,n,p);
        z_M = Z(m,n,p);

        % Loop on each discretized segment of Gamma PkPk+1
        % ALEX: cross product divided by mag R?
        % ALEX: x_P,y_P,z_P correspond to Gamma (source) pts
        %       x_M, y_M, z_M correspond to point of interest (B location)
        %length_XP= length(x_P)-1
        for k = 1:length(x_P)-1
            PkM3 = (sqrt((x_M-x_P(k))^2 + (y_M-y_P(k))^2 + (z_M-z_P(k))^2))^3; % source to point of interest
            DBx(k) = ((y_P(k+1)-y_P(k))*(z_M-z_P(k)) - (z_P(k+1)-z_P(k))*(y_M-y_P(k)))/PkM3; % cross product x
            DBy(k) = ((z_P(k+1)-z_P(k))*(x_M-x_P(k)) - (x_P(k+1)-x_P(k))*(z_M-z_P(k)))/PkM3; % cross product y
            DBz(k) = ((x_P(k+1)-x_P(k))*(y_M-y_P(k)) - (y_P(k+1)-y_P(k))*(x_M-x_P(k)))/PkM3; % cross product z
        end
        % Sum
        BX(m,n,p) = BX(m,n,p) + mu0*I/4/pi*sum(DBx);
        BY(m,n,p) = BY(m,n,p) + mu0*I/4/pi*sum(DBy);
        BZ(m,n,p) = BZ(m,n,p) + mu0*I/4/pi*sum(DBz);

        end
    end
end
%end
  
%% Plot B/|B|
figure(1)
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')
%axis tight
%%
figure(44)
quiver(Y,Z,BY./normB,BZ./normB,'b')
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