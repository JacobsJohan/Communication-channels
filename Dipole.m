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
conf.nrOfFrames     = 50;
conf.Resolution_X   = 300;
conf.Resolution_Y   = 300;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

%% Initialising the different sources 
Source = struct;

f = 900e6; %freq of source
sourceX = 0.95; %Horizontal
sourceY = 0.95;   %Vertical
% Source = addSource( Source,conf,sourceY,sourceX,f,(1-2*(pi*f*T).^2).*exp(-(pi*f*T).^2)); %Gaussian pulse
Source = addSource( Source,conf,sourceY,sourceX,f,zeros(1,numel(T)));
fi=500e06;
i=sin(2*pi*fi*T);
%% Filling the field with objects

x1=0.5;
y1= 0.25;
y2=0.25+0.01+1/12;
thick=0.01;
height=1/6;
field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,x1,y2,thick,height,conf);
field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,x1,y1,thick,height,conf);
%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];
field(1).Sig = draw_rectangle(field(1).Sig,10^10,x1,y2,thick,height,conf);
field(1).Sig = draw_rectangle(field(1).Sig,10^10,x1,y1,thick,height,conf);
%% Prepare video
v = VideoWriter('Dipole','MPEG-4');
v.Quality = 100; 
open(v);
figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])

%Interpolate spatial properties for plotting
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

%Calculate constants used in FDTD algorithm
CHXH = (1-SIGm(:,end-1).*deltat/2./MUXrel)./(1+SIGm(:,end-1).*deltat/2./MUXrel);
CHXE = (Sc/Z0)./MUXrel * 1./(1+SIGm(:,end-1).*deltat/2./MUXrel);
CHYH = (1-SIGm(end-1,:).*deltat/2./MUYrel)./(1+SIGm(end-1,:).*deltat/2./MUYrel);
CHYE = (Sc/Z0)./MUYrel * 1./(1+SIGm(end-1,:).*deltat/2./MUYrel);
CEZE = (1-deltat./2*SIG./EPSrel) ./ (1+deltat/2*SIG./EPSrel);
CEZH = 1./(1+deltat/2*SIG./EPSrel) * Z0*Sc./EPSrel;

%Initializing results
results=struct;
results.Ez=field(1).Ez;%All 3 are initialized with right dimensions
... and all zero
results.Hx=field(1).Hx;
results.Hy=field(1).Hy;

%Initializing memoryMatrices to store results from previous time-instance
prev=struct;
prev.Ez=field(1).Ez;
prev.Hx=field(1).Hx;
prev.Hy=field(1).Hy;

% Prepare for plotting
[M,N] = size(results.(conf.ToPrint));
[X,Y] = meshgrid(linspace(0,conf.x_length,N), linspace(0,conf.y_length,M));


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
    ToPrintq=results.(conf.ToPrint);%Choose field to plot
    absMaxToPrint = max(ToPrintq(:));
% Print
    disp(['Frame: ',num2str(i),' / ',num2str(conf.nrOfFrames)])
    surf(X,Y,ToPrintq,...
            'LineStyle','none',...
            'FaceColor','interp');
    hold on 
    %Display relative permittivity
    surf(   Xq2,...
            Yq2,...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',EPSrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','red');
    %Display relative permeability
    surf(   Xq2,...
            Yq2,...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint+0.1,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',MUrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','blue');
    hold off
    colorbar;
    if conf.ToPrint=='Ez'
        caxis([-0.5,0.5])
    else
        caxis([-0.001,0.001])
    end
    view(2)
    frame = getframe;
    writeVideo(v,frame);

    %Save current fields as old fields for next iteration
    prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
end
%% Free videofile
close(gcf)
close(v)

%% Release path
disp(['Note: if you look at the Hx and Hy field you see that they are not isotropic'])
rmpath(genpath(pwd))