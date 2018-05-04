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
conf.nrOfFrames     = 25;
conf.Resolution_X   = 300;
conf.Resolution_Y   = 300;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

%% Initialising the different sources 
Source = struct;

f = 900e6;
sX=conf.x_length/2;
sY=conf.y_length/2;
Source = addSource( Source,conf,sX,sY,f,sin(2*pi*f*T) );

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
       results.Ez(sourceXloc,sourceYloc) = sourceValue;
    end
    
% %Plot calculated fields
% 
% % Prepare for plotting
% 
% 
%     [M,N] = size(results.(conf.ToPrint));
%     [X,Y] = meshgrid(linspace(0,conf.x_length,N), linspace(0,conf.y_length,M));
% 
%     ToPrintq=results.(conf.ToPrint);
%     temp = ToPrintq(:,:,20:end);
%     absMaxToPrint = max(ToPrintq(:));
% 
% 
% % Print
%     disp(['Frame: ',num2str(i),' / ',num2str(conf.nrOfFrames)])
%     surf(X,Y,ToPrintq,...
%             'LineStyle','none',...
%             'FaceColor','interp');
%     hold on 
%     surf(   Xq2,...
%             Yq2,...
%             ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint,...
%             'FaceAlpha','interp',...
%             'AlphaDataMapping','none',...
%             'AlphaData',EPSrelalpha(:,:,1),...
%             'LineStyle','none',...
%             'FaceColor','red');
%     surf(   Xq2,...
%             Yq2,...
%             ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint+0.1,...
%             'FaceAlpha','interp',...
%             'AlphaDataMapping','none',...
%             'AlphaData',MUrelalpha(:,:,1),...
%             'LineStyle','none',...
%             'FaceColor','blue');
%     hold off
%     colorbar;
%     caxis([-0.5,0.5])
%     view(2)
%     frame = getframe;
%     writeVideo(v,frame);

%Save current fields as old fields for next iteration
    prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
    
  %Track 1st nonzero source value and its evolution over 1 line until reflection  
    if i>=2 && i < sz(1)/2 && i< sz(2)/2
        E_fields(:,:,i-1)=prev.Ez; %Used to see the evolution manually
        E_temp=[E_temp prev.Ez(meter2index(sX,conf)+i-1,meter2index(sY,conf)+1)]; %Save E_field evolution on 1 row. 
    end
    % H field follows a time instance later
    if i>=3 && i < sz(1)/2+1 && i< sz(2)/2+1
        Hx_temp=[Hx_temp prev.Hx(meter2index(sX,conf)+i-2,meter2index(sY,conf)+1)];
        Hy_temp=[Hy_temp prev.Hy(meter2index(sX,conf)+i-2,meter2index(sY,conf)+1)];    
    end
    end
end

%% Calculation of energy
H_total = [Hx_temp' Hy_temp' zeros(size(Hx_temp,2),1)];
E_total = [zeros(size(E_temp,2),1) zeros(size(E_temp,2),1) E_temp'];
P = abs(sqrt(H_total(:,1).^2+H_total(:,2).^2)).*abs(E_total(:,3));%./(4*pi*1e-7);
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
plot(ratio)
figure
plot(P)
hold on
plot(P(1)./r./r,'*')
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