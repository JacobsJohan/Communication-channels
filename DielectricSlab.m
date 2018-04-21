clear
close all
clc
addpath(genpath(pwd))

% NOTE:
% For comments and info on the general working of the code please check teh
% file main_commented. The structure is a bit different but the principles
% are the same. The comment mentioned here are only those specifically for
% this script.

%% Dielectric slab
% The dielectic slab consists of 3 layers of dielectrics, each a different
% dielectric constant and therefore a different refraction index. To get
% internal reflection, the dielectric index of the middle layer must be
% (slightly) larger than the other 2. Here the assumption is made that all
% relative permeablilities are the same.

%% Initialising the fields for simulation 
conf.fmax           = 2e09;
conf.x_length       = 1;
conf.y_length       = 1;
conf.nrOfFrames     = 150;
conf.Resolution_X   = 300;
conf.Resolution_Y   = 300;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

%% Initialising the different sources 
Source = struct;

f = 2e09;

%0° angle
% xloc = 0.1;
% yloc = 0.5; 
% internal=0;
 
% %some angle
xloc = 0.1;
yloc = 0.3; 
internal=0;
 
% %other angle
% xloc = 0.4;
% yloc = 0.5; 
% internal=1; %We place the source at input to get high input angles. Therefore we clear shielding
Source = addSource( Source,conf,yloc,xloc,f,sin(2*pi*f*T) );

%% Filling the field with objects
middleX=0.7;
slabL=0.7;
slabT=0.1;

% NOTE:
% The dimensions of the slab are much to large to be realistic. However due
% to the goal and the nature of the simulation, it is okay. If we would
% make it small, it would not be properly visible or the reflections would
% bother the main goal of the simulation too much. The main properties 
% still hold so we can use this as a model.

%To only study the influence of the right angels, we will shield the side
%of the waveguide

if ~internal
field(1).EpsRel = draw_rectangle(field(1).EpsRel,200,middleX,0.655,slabL,0.01,conf); %Isolating shield around wire
field(1).EpsRel = draw_rectangle(field(1).EpsRel,200,0.7-0.35-0.005,0.6,0.01,0.11,conf); %Isolate upper input
field(1).EpsRel = draw_rectangle(field(1).EpsRel,200,0.7-0.35-0.005,0.4,0.01,0.11,conf); %Isolate bottom input
field(1).EpsRel = draw_rectangle(field(1).EpsRel,200,middleX,0.345,slabL,0.01,conf); %Isolating shield around wire
end

field(1).EpsRel = draw_rectangle(field(1).EpsRel,2,middleX,0.6,slabL,slabT,conf); %Upper
field(1).EpsRel = draw_rectangle(field(1).EpsRel,5,middleX,0.5,slabL,slabT,conf); %Middle
field(1).EpsRel = draw_rectangle(field(1).EpsRel,2,middleX,0.4,slabL,slabT,conf); %Bottom


%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];

%% Prep video
v = VideoWriter('DSangle','MPEG-4');
v.Quality = 100; 
v.FrameRate = 5;
open(v);
figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])

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

for i=1:conf.nrOfFrames-1
%     tempfield = FDTDMaxwellCore2(tempfield,field,conf,Source );

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
    [X,Y] = meshgrid(linspace(0,conf.x_length,N), linspace(0,conf.y_length,M));

    ToPrintq=results.(conf.ToPrint);
    temp = ToPrintq(:,:,20:end);
    minToPrint = min(temp(:));
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

%% Draw and export to movie 
rmpath(genpath(pwd))
