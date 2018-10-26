clear all; close all; clc;
u0 = 4*pi*10^(-7);
I=1;
b=0.5e-3;
N=60;
R  = 10e-3;%linspace(-100e-3,100e-3,N);
th = linspace(0,pi,N);

%{
%B(:,:) = [((u0.*I.*b.^2)./(4.*R(:).^3)).*2.*cos(th(:)), ((u0.*I.*b.^2)./(4.*R(:).^3)).*sin(th(:))];

for tn=1:N
    for rn=1:numel(R)
        Br(rn,tn) = ((u0*I*b^2)/(4*(R(rn))^3)) *2*cos(th(tn));
        Bt(rn,tn) = ((u0*I*b^2)/(4*(R(rn))^3)) *sin(th(tn));
    end
end
[x,y] = meshgrid(R.*cos(th),R.*sin(th));
%quiver(x,y,B(1,:),B(:,1));
%}

[Br,Bt] = calcB(I,b,R,N);
[x,y] = pol2cart(th,R);
%x=x+b;
%y=y+b;
%q=quiver(x/1e-3,y/1e-3,Br,Bt);

h=quiver(x,y,Br,Bt);

%q=quiver(x,y,Br,Bt);

%[Br2,Bt2] = calcB(-I,b,R,N);

%{
for rn=1:N
    [x0(:,rn),y0(:,rn)] = pol2cart(th,R(rn));
    [Bx(:,rn),By(:,rn)] = pol2cart(Bt(:,rn),Br(:,rn));
    %[Btx(:,rn),Bty(:,rn)] = pol2cart(Bt,Br(rn));
    
end
%[x1,y1] = pol2cart(th,R(1));

%
%[x0 y0] = meshgrid(x,y);
%[Br0,Bt0] = meshgrid(Br,Bt);
%figure(1)
%q1=quiver(x0,y0,Bx,By);
%figure(2)
%q1=quiver(x0,y0,Br0,Bt0);

figure(3)
p=polar(th,Br)
figure(4)
p1=polar(th,Bt);
%}