clc
clear
close all
addpath(genpath(pwd))

% NOTE:
% For comments and info on the general working of the code please check the
% file 'main_commented.m'. The structure is a bit different but the principles
% are the same. Not all general comments are present in this file.


% NOTE:
% This script is for plotting the KE-simulation. The calculations for
% different v are done in 'Knife_Edge.m'
%% Initialising the fields for simulation 
conf.fmax           = 900e6;
conf.x_length       = 15;
conf.y_length       = 15;
conf.nrOfFrames     = 500;
conf.Resolution_X   = 500;
conf.Resolution_Y   = 500;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

%Defining locations and heights of TX,RX
factor = 3;
sloc=2+factor/2;
Recloc= 8+factor;
source_h=2.5;
rec_h=1.5;

%% Initialising the different sources 
Source = struct;

f = 900e6;
Source = addSource( Source,conf,source_h,sloc,f,exp(-((T-30*conf.deltat)./10./conf.deltat).^2));%sin(2*pi*f*T));
%% Filling the field with objects
buildx=3+factor;
% buildx = 10,5;

%Knife-edge
field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,buildx,2.5,0.1,5,conf);
% Note you can make it thinner but than you have to increase EPSrel

%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];

%Knife-edge
field(1).Sig=draw_rectangle(field(1).Sig,10e10,buildx,2.5,0.1,5,conf);


%% Prep video
v = VideoWriter('KE-model_10.5','MPEG-4');
v.Quality = 100; 
open(v);
figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])

EpsRel = field(1).EpsRel();
MuRel = field(1).MuRel();
SIG=field(1).Sig();

[M,N,T] = size(EpsRel);
[X2,Y2]         = meshgrid(      linspace(0,conf.x_length,M),...
                                    linspace(0,conf.y_length,N)...
                                    );


[Xq2,Yq2]       = meshgrid(     linspace(0,conf.x_length,conf.Resolution_X),...
                       linspace(0,conf.y_length,conf.Resolution_Y)); 

EPSrelalpha     = (interp2(X2,Y2,EpsRel,Xq2,Yq2) -1) /   (max(EpsRel(:))-1)/2.5;
MUrelalpha      = (interp2(X2,Y2,MuRel,Xq2,Yq2)  -1) /   (max(MuRel(:))-1)/2.5;

SIGrelalpha      = (interp2(X2,Y2,SIG,Xq2,Yq2)  -1) /   (max(SIG(:))-1)/2.5;
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
    
%Plot calculated fields

% Prepare for plotting


    [M,N] = size(results.(conf.ToPrint));
    [X,Y] = meshgrid(linspace(0,conf.x_length,N),...
        linspace(0,conf.y_length,M));

    ToPrintq=results.(conf.ToPrint);
    absMaxToPrint = max(ToPrintq(:));


% Print
    disp(['Frame: ',num2str(i),' / ',num2str(conf.nrOfFrames)])
    surf(X,Y,ToPrintq,...
            'LineStyle','none',...
            'FaceColor','interp');
    hold on 
    surf(   Xq2,...
            Yq2,...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',EPSrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','red');
    surf(   Xq2,...
            Yq2,...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint+0.1,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',MUrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','blue');
        surf(   Xq2,...
            Yq2,...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint+0.1,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',SIGrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','red');
    hold off
    colorbar;
    caxis([-0.5,0.5])
    view(2)
    frame = getframe;
    writeVideo(v,frame);

%Save current fields as old fields for next iteration
    prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
end
%% Free videofile
close(gcf)
close(v)

%% Free path
rmpath(genpath(pwd))