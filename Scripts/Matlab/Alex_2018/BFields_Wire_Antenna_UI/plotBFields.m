
function plotBFields(X,Y,Z,BX,BY,BZ)

%% Plot B/|B|
figure(2)
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    BXn=BX./normB; BYn=BY./normB; BZn=BZ./normB;
    quiver3(X,Y,Z,BXn,BYn,BZn,'b');
    %axis tight
    
%%
    %{
figure(3)
    quiver(Y,X,BYn,BXn);

%% Plot Bz on the volume
figure(4), hold on, box on, grid on
	plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
	slice(X,Y,Z,BZ,[0],[],[-1,0,1]), colorbar % plot Bz
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
    view(3), axis equal, axis tight
    caxis([-0.5,0.5]*1e-5)
    
%% Plot some flux tubes
figure(5), hold on, box on, grid on
	plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
	[X0,Y0,Z0] = ndgrid(-1.5:0.5:1.5,-1.5:0.5:1.5,-2); % define tubes starting point        
	htubes = streamtube(stream3(X,Y,Z,BX,BY,BZ,X0,Y0,Z0), [0.2 10]);
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Some flux tubes')
    view(3), axis equal, axis tight
    set(htubes,'EdgeColor','none','FaceColor','c') % change tube color
    camlight left % change tube light
%}
end