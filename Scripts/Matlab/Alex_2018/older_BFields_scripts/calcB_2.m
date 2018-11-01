function [Br,Bz,r,z] = calcB_2(I,b,N)
    % u0 is permeability constant/magnetic constant of free space
    % units: H
    u0 = 4*pi*10^(-7);
    
    z = linspace(-10e-3,10e-3,N);
    r = linspace(-10e-3,10e-3,N);
    
    for zn=1:N
        for rn=1:N
            mag(zn,rn) = (u0*I)/(2*(z(zn)^2 + (r(rn)-b)^2 )^(3/2)); 
        end
    end
  
    for zn=1:N
        for rn=1:N
            % mag(zn,rn) = (u0*I) / 2*((z(zn))^2 + (r(rn)-b)^2)^(3/2);
            % mag(zn,rn) = (u0*I)/(2*(z(zn)^2 + (r(rn)-b)^2 )^(3/2)); 
            Br(zn,rn) = mag(zn,rn)*(-z(zn)*b);
            Bz(zn,rn) = mag(zn,rn)*(b*(r(rn))-b);
        end
    end
    
    figure(2)
    q=quiver(r,z,Br,Bz,8);%This quiver plots both matrices as vectors
    xlabel('r-axis')
    ylabel('z-axis')
    %axis([-0.01 0.01 -0.01 0.01])
end