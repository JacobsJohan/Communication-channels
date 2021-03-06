clc
clear
close all
addpath(genpath(pwd))

%% Initialising the fields for simulation 
conf.fmax           = 900e6;
conf.x_length       = 3;
conf.y_length       = 3;
conf.nrOfFrames     = 50;
conf.Resolution_X   = 300;
conf.Resolution_Y   = 300;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

%% Initialising the different sources 
Source = struct;

f = 900e6;
Source = addSource( Source,conf,1.5,1.5,f,sin(2*pi*f*T) );
%% Simulating losses
field(1).SigM(:) = 0;
field(1).Sig(:) = 0;

%% Speedtest point
x=2.5;
y=1.5;
xind = meter2index(2.5,conf);
yind = meter2index(1.5,conf);

checking=1;%Boolean 

%% Prep video
v = VideoWriter('Output','MPEG-4');
v.Quality = 100; 
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


    [M,N] = size(results.Ez);
    [X,Y]         = meshgrid(     linspace(0,conf.x_length,M),...
                                    linspace(0,conf.y_length,N)...
                                    );

    ToPrintq=results.Ez;
    temp = ToPrintq(:,:,20:end);
    % maxToPrint = max(temp(:));
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
%     text(X2(end,end)/10,Y2(end,end)/10,absMaxToPrint+0.2,['time = ',num2str(Zq(1,1,i)),'s']);
    hold off
    colorbar;
%     zlim([minToPrint,absMaxToPrint+0.2]);
    caxis([-0.5,0.5])
    view(2)
    frame = getframe;
    writeVideo(v,frame);

%Save current fields as old fields for net iteration
    prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
    
    %Check if nonzero
    if checking && abs(prev.Ez(yind,xind))>1e-10
        iteration=i;
        checking=0;
    end
end
%% Free videofile
close(gcf)
close(v)

%% Speed validation
% How will we validate the speed?
% We can take a random point in space and check after which nr of frames it
% becomes different from zero. From there on we can calculate back the
% speed because we know the time step and spatial resolution

% Let's center our source and take a point at 1m to the right. This yields
% x=2.5,y=1.5 Indices etc are calculated above

% time=frame-2 because in the first frame, everything is being initialized.
% In the second frame the source starts transmitting, but there is no 
% propagation yet. This means we have 2 all zero frames. Since we check the
% first value that is nonzero, we have to drop 2 frames
% The propagation distance is the difference in indices between source and
% point 
speed = abs(meter2index(1.5,conf)-xind)/(iteration-2);%indexdifference/time

% Remember conf.deltat = conf.delta/3e8/sqrt(2);  

%% Draw and export to movie 
% plotAndToVid2('Output/simulation2',field,conf,0.5,-0.5)
rmpath(genpath(pwd))