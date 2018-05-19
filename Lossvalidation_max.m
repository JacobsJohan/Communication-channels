clc
clear
close all
addpath(genpath(pwd))

% NOTE:
% For comments and info on the general working of the code please check the
% file 'main_commented.m'. The structure is a bit different but the principles
% are the same. Not all general comments are present in this file.

%% Initialising the fields for simulation 
conf.fmax           = 900e6;
conf.x_length       = 1;
conf.y_length       = 1;
conf.nrOfFrames     = 100;
conf.Resolution_X   = 300;
conf.Resolution_Y   = 300;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

%% Initialising the different sources 
Source = struct;

f = 900e6;
sX=conf.x_length/2;
sY=conf.y_length/2;
% Source = addSource( Source,conf,sX,sY,f,sin(2*pi*f*T) );
Source = addSource( Source,conf,sX,sY,f,3*exp(-((T-30*conf.deltat)./10./conf.deltat).^2) );


%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];

%% Prep video
% v = VideoWriter('Output','MPEG-4');
% v.Quality = 100; 
% open(v);
% figure()
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])

EpsRel = field(1).EpsRel();
MuRel = field(1).MuRel();

[M,N,T] = size(EpsRel);
[X2,Y2]         = meshgrid(      linspace(0,conf.x_length,M),...
                                    linspace(0,conf.y_length,N)...
                                    );


[Xq2,Yq2]       = meshgrid(     linspace(0,conf.x_length,conf.Resolution_X),...
                       linspace(0,conf.y_length,conf.Resolution_Y)); 

EPSrelalpha     = (interp2(X2,Y2,EpsRel,Xq2,Yq2) -1) /   (max(EpsRel(:))-1)/2.5;
MUrelalpha      = (interp2(X2,Y2,MuRel,Xq2,Yq2)  -1) /   (max(MuRel(:))-1)/2.5;


%% Simulate field

deltat  = conf.deltat;
delta   = conf.delta;
c       = 3e8;              %Speed of light
mu0     = pi*4e-7;          %Permeabillity of free space 
eps0    = 1/mu0/c/c;        %Permitivity of free space 
Z0      = sqrt(mu0/eps0);   %Free space impedance 
Sc      = c*deltat/delta;  	%Courant number

MUXrel  = field(1).MuRel(1:2:end,2:2:end);
MUYrel  = field(1).MuRel(2:2:end,1:2:end);
SIG     = field(1).Sig(1:2:end,1:2:end);
SIGm    = field(1).SigM(1:2:end,1:2:end);
EPSrel  = field(1).EpsRel(1:2:end,1:2:end);

CHXH = (1-SIGm(:,end-1).*deltat/2./MUXrel)./(1+SIGm(:,end-1).*deltat/2./MUXrel);
CHXE = (Sc/Z0)./MUXrel * 1./(1+SIGm(:,end-1).*deltat/2./MUXrel);
CHYH = (1-SIGm(end-1,:).*deltat/2./MUYrel)./(1+SIGm(end-1,:).*deltat/2./MUYrel);
CHYE = (Sc/Z0)./MUYrel * 1./(1+SIGm(end-1,:).*deltat/2./MUYrel);
CEZE = (1-deltat./2*SIG./EPSrel) ./ (1+deltat/2*SIG./EPSrel);
CEZH = 1./(1+deltat/2*SIG./EPSrel) * Z0*Sc./EPSrel;

results=struct;
results.Ez=field(1).Ez;
results.Hx=field(1).Hx;
results.Hy=field(1).Hy;

prev=struct;
prev.Ez=field(1).Ez;
prev.Hx=field(1).Hx;
prev.Hy=field(1).Hy;

E_temp=[]; %Save E-field values
Hx_temp=[]; %Save Hx-field values
Hy_temp=[]; %Save Hy-field values
sz=size(prev.Ez);


%% CheckMax
checking=1;%Boolean 
checkMax = 1;
sourceMax = 1;

prevval = -50;
prevsrc = -50;

for i=1:conf.nrOfFrames-1
%%Calculate new fields
    disp([num2str(i),' / ',num2str(conf.nrOfFrames)])
    results.Hx =    CHXH.*prev.Hx -...
                    CHXE.*(prev.Ez(:,(1:end-1)+1) - prev.Ez(:,1:end-1)  );
    results.Hy =    CHYH.*prev.Hy +...
                    CHYE.*(prev.Ez((1:end-1)+1,:) - prev.Ez(1:end-1,:));
                    
    results.Ez(2:end-1,2:end-1) =    CEZE(2:end-1,2:end-1).*     prev.Ez(2:end-1,2:end-1) +...
                                        CEZH(2:end-1,2:end-1).*(    (results.Hy(2:end,2:end-1)   -   results.Hy((2:end)-1,2:end-1))-...
                                                                    (results.Hx(2:end-1,2:end)   -   results.Hx(2:end-1,(2:end)-1)));
%Update sources
    for s=1:numel(Source)
       sourceXloc = round(Source(s).X_loc(i)/conf.delta)+1;
       sourceYloc = round(Source(s).Y_loc(i)/conf.delta)+1;
       sourceValue = Source(s).value(i);
        % Check if source increasing
       if sourceMax
           if sourceValue >= prevsrc 
                prevsrc = sourceValue;
           else %Value decreases, meaning previous was max
               sourceMax = 0;
               src_i = i-1;
           end
       end
       results.Ez(sourceXloc,sourceYloc) = sourceValue;
    end
    

%Save current fields as old fields for next iteration
    prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
    
    E_fields(:,:,i)=prev.Ez;
    Hx_fields(:,:,i)=prev.Hx;
    Hy_fields(:,:,i)=prev.Hy;
end
%Get line
E_line = squeeze(E_fields(meter2index(sX,conf)+1,meter2index(sY,conf)+1:end,src_i:end));
Hx_line = squeeze(Hx_fields(meter2index(sX,conf)+1,meter2index(sY,conf)+1:end,src_i-1:end));
Hy_line = squeeze(Hy_fields(meter2index(sX,conf)+1,meter2index(sY,conf)+1:end,src_i-1:end));
% for i=1:size(Hx_line,1)
%    E_temp(i)= E_line(i,i); 
%    Hx_temp(i)= Hx_line(i,i); 
%    Hy_temp(i)= Hy_line(i,i); 
% end
% for i=1:size(Hx_line,1)
%    E_temp(i)= max(E_line(i,:)); 
%    Hx_temp(i)= Hx_line(i,i); 
%    Hy_temp(i)= Hy_line(i,i); 
% end

for i=1:size(Hx_line,1)
   [E_temp(i) j]= max(E_line(i,:)); 
   Hx_temp(i)= Hx_line(i,j); 
   Hy_temp(i)= Hy_line(i,j); 
end
%% Calculation of energy
P = abs(sqrt(Hx_temp'.^2+Hy_temp'.^2)).*abs(E_temp');%./(4*pi*1e-7);
% In theory you should divide by mu but here it is not necessary because we
% are only interested in the proportionality to 1/r^2.

% %% Free videofile
% close(gcf)
% close(v)
% 

%% Calculating ratios
ratio = E_temp(2:end)./ E_temp(1:end-1);
r = 1:numel(P)+1;
mainRatio = E_temp(1:end)./ E_temp(1); %Compare with starting value;
%% Plotting
figure
yyaxis left
plot(E_temp)
hold on
plot(E_temp(1)./r,'-o')
yyaxis right
plot(Hx_temp)
plot(abs(Hy_temp))
plot(Hx_temp(1)./r,'*')
plot(abs(Hy_temp(1))./r,'o')
legend('Ez','1/r','Hx','Hy','1/r','1/r')
figure
plot(P)
hold on
plot(P(1)./r./r,'*')
% plot(P(1)./r,'*')
legend('power ifo distance','P(1)/r^2')
hold off
xlabel('Distance from source (n*\delta)')
ylabel('P = E.B/\mu_0')
title('Evolution of power in function of distance to source.')
% figure
% plot(mainRatio)
%% Output
disp(['Ratio of E-field compared to itself 1 time instance ago and on previous location:' num2str(ratio)])
disp(['Ratio of E-field compared to original source value:' num2str(mainRatio)])
%% Draw and export to movie 
rmpath(genpath(pwd))