clc
clear
close all
addpath(genpath(pwd))
% NOTE:
% For comments and info on the general working of the code please check the
% file 'main_commented.m'. The structure is a bit different but the principles
% are the same. Not all general comments are present in this file.

%Only the speed of the E-field is checked.
%% Initialising the fields for simulation 
conf.fmax           = 900e6;
conf.x_length       = 3;
conf.y_length       = 3;
conf.nrOfFrames     = 100;
conf.Resolution_X   = 300;
conf.Resolution_Y   = 300;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

%% Initialising the different sources 
Source = struct;
f = 900e6;
Source = addSource( Source,conf,1.5,1.5,f,sin(2*pi*f*T) );
%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];


%% Changing epsrel for entire field
% n = sqrt(epsrel*murel);
% We want n of 1.5, n=c/v, murel = 1;
% Select here a wanted n
n=1;
n=1.2;
epsrel = n^2;
field(1).EpsRel(:) = epsrel;

%% Setting treshold for detection
%  Due to the nature of the simulation, the field is updated every
%  iteration, which is a fixed 'speed'. Depending on the value of epsrel
%  the values in a point where there should not yet be a wave (if n~=1) is
%  already nonzero (or larger than 1e-15 due to Matlab accuracy). Therefore
%  one should alter the treshold. In this kind of simulations it is only
%  encouraged to check the speed of light. One could say that if the field
%  is small enough, we can neglect it so 'it is not yet there', but picking
%  a treshold in that case is hard. So to be mathematically correct and
%  to avoid suspicious assumptions, one should only validate the speed if
%  n=1. We checked other speeds with a treshold of 1e-7 and 1e-8 which gave good
%  results in a broad range of realistic n~=1 values.

if n == 1
    treshold = 1e-15;
elseif n <1.3
    treshold = 1e-8;
else
    treshold = 1e-7;
end

% Circumvension: look for maxima 

%%  Speedtest point
x = 2.5;
x2 = 2;
y = 1.5;
xind = meter2index(x,conf);
yind = meter2index(y,conf)+1;
x2ind = meter2index(x2,conf);

checking=1;%Boolean 
checkMax = 1;
sourceMax = 1;

prevval = -50;
prevsrc = -50;


%% Prep video
% v = VideoWriter('SpeedVal','MPEG-4');
% v.Quality = 100; 
% open(v);
% figure()
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])
% 
% EpsRel = field(1).EpsRel();
% MuRel = field(1).MuRel();
% 
% [M,N,T] = size(EpsRel);
% [X2,Y2]         = meshgrid(      linspace(0,conf.x_length,M),...
%                                     linspace(0,conf.y_length,N)...
%                                     );
% 
% 
% [Xq2,Yq2]       = meshgrid(     linspace(0,conf.x_length,conf.Resolution_X),...
%                        linspace(0,conf.y_length,conf.Resolution_Y)); 
% 
% EPSrelalpha     = (interp2(X2,Y2,EpsRel,Xq2,Yq2) -1) /   (max(EpsRel(:))-1)/2.5;
% MUrelalpha      = (interp2(X2,Y2,MuRel,Xq2,Yq2)  -1) /   (max(MuRel(:))-1)/2.5;


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


%% Loop 
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
               src_i = i;
           end
       end
       results.Ez(sourceXloc,sourceYloc) = sourceValue;
    end
    
%Plot calculated fields

% % Prepare for plotting
% 
% 
%     [M,N] = size(results.(conf.ToPrint));
%     [X,Y] = meshgrid(linspace(0,conf.x_length,N),...
%         linspace(0,conf.y_length,M));
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
    
    %Check if nonzero
    if checking && abs(prev.Ez(yind,xind))>treshold
        iteration=i; %This is the frame where the E-field in the point becomes nonzero.
        checking=0; %Frame found so checking can stop
    end
    save(1,i)=prev.Ez(yind,x2ind);
    save(2,i)=prev.Ez(yind,yind);
    if checkMax
        
       if  prev.Ez(yind,x2ind)>= prevval
           prevval = prev.Ez(yind,x2ind);  
       else
           rec_i = i-1;
           checkMax = 0;
       end
    end
end
% %% Free videofile
% close(gcf)
% close(v)

%% Speed validation
% How will we validate the speed?
% We can take a random point in space and check after which nr of frames it
% becomes different from zero. From there on we can calculate back the
% speed because we know the time step and spatial resolution

% Let's center our source and take a point at 1m to the right. This yields
% x=2.5,y=1.5 Indices etc are calculated above.

%time=frame-2 because first frame is all zero and
%because in the first frame the source value is placed after calculations
%of the fields and the first sourcevalue is 0. This means we have 2 all zero frames. Since we check for the
%first source value that is nonzero, we have to drop 2 frames.
%The propagation distance is the difference in indices between source and
%point 

speed = abs(meter2index(1.5,conf)+1-xind)/(iteration-2);%indexdifference/time
speed_max = abs(meter2index(1.5,conf)+1-x2ind)/(rec_i-src_i);

disp(['The ratio of the simulation speed and the speed of light is: ' num2str(speed)])
if n==1
    if speed == 1
        disp('Speed in simulation is speed of light')
    else
        disp(['Error in simulation'])
    end
else
    disp(['Speed is v = ' num2str(c*speed) ' m/s'])
    disp(['v/c ratio = ' num2str(speed) ', 1/n = ' num2str(1/n)]);
end
disp(['The ratio of the simulation speed and the speed of light looking at propagation of max is: ' num2str(speed_max)])
%% Free path
rmpath(genpath(pwd))