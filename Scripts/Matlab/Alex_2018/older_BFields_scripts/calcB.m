% Alexander Moreno
% 3-31-2018
%
% Calc B fields of loop antenna 
%
% Input Parameters:
% --
% I
% description: current through the loop
% units: A
% scalar
%--
% b
% description: radius of loop 
% untits: m
%
%--
% R
% description: origin to the point of interest
% units: m
% 1d vector

function [Br,Bt] = calcB(I,b,R,N)
    % N is the number of segments (aka resolution)
    % N = 360;
    % th is theta  
    % units: rads
    th = linspace(0,pi,N);
    % u0 is permeability constant/magnetic constant of free space
    % units: H
    u0 = 4*pi*10^(-7);
    
    for tn=1:numel(th)
        for rn=1:numel(R)
            Br(rn,tn) = ((u0*I*b^2)/(4*(R(rn))^3)) *2*cos(th(tn));
            Bt(rn,tn) = ((u0*I*b^2)/(4*(R(rn))^3)) *sin(th(tn));
        end
    end
    
end