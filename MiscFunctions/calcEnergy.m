function [] = calcEnergy(E,x,y,nrOfFrames)
%Determines the evolution of the magnitude of the electric field in a given
%point i.f.o. time
energy=0;

for i=2:nrOfFrames
   energy=[energy sqrt(E(x,y,i)^2)];
end

%Calculate evolution of max E-field in a fixed direction (x or y)
E=abs(E);
s=size(E);

%Excuses voor de for loops Mathias maar zo is het veel beter te volgen
for j=x:s(1)
    maxesX(j-x+1)=max(E(j,y,:));%These values are in abs  
end
%Works for the moment only for rectangular windows
for i=y:s(1)
    maxesY(i-y+1)=max(E(x,i,:));
end
% numel(1:s(1)-x)
% numel(maxesX)
% 
% maxesX=permute(maxesX,[3,2,1]);
% maxesY=permute(maxesY,[3,2,1]);

figure
plot(energy)
%Coordinates in title moeten nog aangepast worden naar de genormaileerde
%waarden.
title(['Norm of the E-field in 1 fixed point with coordinates: x ='...
                        num2str(x) ' and y = ' num2str(y)])
figure
yyaxis left
plot(1:s(1)-x+1,maxesX)
yyaxis right
plot(1:s(1)-y+1,maxesY)
title(['Evolution of the norm of E-field starting in (x,y)=(' num2str(x)...
                                   ',' num2str(y) ') to right and bottom'])
legend('X-direction (?)','Y-direction (?)')
end

