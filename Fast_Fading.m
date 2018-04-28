clc
clear
close all
addpath(genpath(pwd))

%% Initialising the fields for simulation 
conf.fmax           = 2e09;
conf.nrOfFrames     = 1000;
conf.Resolution_X   = 500;
conf.Resolution_Y   = 500;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit_fastfading( conf );

%% Initialising the different sources 
Source = struct;

f = 2e09; %freq of source
sourceX = conf.x_length/2; %Horizontal
sourceY = conf.y_length/2  ;   %Vertical
Source = addSource( Source,conf,sourceY,sourceX,f,sin(2*pi*f*T) );

%% Selecting measurement points

% %Check opgave om beter te hebbe
% 
% % xs=[1/4 1/4 3/4 3/4]*conf.x_length;
% % ys=[1/4 3/4 1/4 3/4]*conf.y_length;
% 
% xs=[0.25 0.75]*conf.x_length;
% ys=[1/4 3/4 ]*conf.y_length;
% 
% indx = meter2index(xs,conf);
% indy = meter2index(ys,conf);

%% Define local areas
boxSize = 30*conf.delta;
y11 = sourceY-50*conf.delta;
y12 = sourceY-50*conf.delta-boxSize;
x11 = sourceX-boxSize/2;
x12 = sourceX+boxSize/2;
box1=[x11 x12 y12 y11];

y21 = 50*conf.delta+boxSize/2;
y22 = 50*conf.delta-boxSize/2;
x21 = sourceX-boxSize/2;
x22 = sourceX+boxSize/2;
box2=[x21 x22 y22 y21];

y31 = 50*conf.delta+boxSize/2;
y32 = 50*conf.delta-boxSize/2;
x31 = 50*conf.delta-boxSize/2;
x32 = 50*conf.delta+boxSize/2;
box3=[x31 x32 y32 y31];

y41 = 50*conf.delta+boxSize/2;
y42 = 50*conf.delta-boxSize/2;
x41 = 375*conf.delta-boxSize/2;
x42 = 375*conf.delta+boxSize/2;
box4=[x41 x42 y42 y41];

%converting to indices
boxind = meter2index([box1;box2;box3;box4],conf);
%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];

%% Prepare video
v = VideoWriter('FF','MPEG-4');
v.FrameRate = 8;
v.Quality = 100; 
open(v);
figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])

%Interpolate spatial properties for plotting
EpsRel = field(1).EpsRel();
MuRel = field(1).MuRel();

[M,N,T] = size(EpsRel);
[X2,Y2]         = meshgrid(linspace(0,conf.x_length,M),...
                                    linspace(0,conf.y_length,N)...
                                    );


[Xq2,Yq2]       = meshgrid(linspace(0,conf.x_length,conf.Resolution_X),...
                       linspace(0,conf.y_length,conf.Resolution_Y)); 

EPSrelalpha = (interp2(X2,Y2,EpsRel,Xq2,Yq2)-1)/(max(EpsRel(:))-1)/2.5;
MUrelalpha  = (interp2(X2,Y2,MuRel,Xq2,Yq2) -1) /(max(MuRel(:))-1)/2.5;


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
CHXH = (1-SIGm(:,end-1).*deltat/2./MUXrel)./...
                                     (1+SIGm(:,end-1).*deltat/2./MUXrel);
CHXE = (Sc/Z0)./MUXrel * 1./(1+SIGm(:,end-1).*deltat/2./MUXrel);
CHYH = (1-SIGm(end-1,:).*deltat/2./MUYrel)./...
                                     (1+SIGm(end-1,:).*deltat/2./MUYrel);
CHYE = (Sc/Z0)./MUYrel * 1./(1+SIGm(end-1,:).*deltat/2./MUYrel);
CEZE = (1-deltat./2*SIG./EPSrel) ./ (1+deltat/2*SIG./EPSrel);
CEZH = 1./(1+deltat/2*SIG./EPSrel) * Z0*Sc./EPSrel;

%Initializing results
results=struct;
%All 3 fields are initialized with right dimensions and all zero
results.Ez=field(1).Ez;
results.Hx=field(1).Hx;
results.Hy=field(1).Hy;

%Initializing memoryMatrices to store results from previous time-instance
prev=struct;
prev.Ez=field(1).Ez;
prev.Hx=field(1).Hx;
prev.Hy=field(1).Hy;

% Prepare for plotting
[M,N] = size(results.(conf.ToPrint));
[X,Y] = meshgrid(linspace(0,conf.x_length,N),...
                                            linspace(0,conf.y_length,M));

% Initialize temp saveFile
% E_save=ones(4,conf.nrOfFrames-1);

for i=1:conf.nrOfFrames-1
%%Calculate new fields
    disp([num2str(i),' / ',num2str(conf.nrOfFrames)])
    results.Hx =    CHXH.*prev.Hx -...
                    CHXE.*(prev.Ez(:,(1:end-1)+1) - prev.Ez(:,1:end-1));
    results.Hy =    CHYH.*prev.Hy +...
                    CHYE.*(prev.Ez((1:end-1)+1,:) - prev.Ez(1:end-1,:));
                    
    results.Ez(2:end-1,2:end-1) = CEZE(2:end-1,2:end-1).*...
        prev.Ez(2:end-1,2:end-1) +CEZH(2:end-1,2:end-1).*...
       ((results.Hy(2:end,2:end-1)-results.Hy((2:end)-1,2:end-1))...
             -(results.Hx(2:end-1,2:end)-results.Hx(2:end-1,(2:end)-1)));
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
    surf(   Xq2,Yq2,...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',EPSrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','red');
    %Display relative permeability
    surf(   Xq2,Yq2,...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint+0.1,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',MUrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','blue');
    %Other properties of the space can also be plotted in similar way as
    %EPSrel and MUrel.
    hold off
    colorbar;
    caxis([-0.5,0.5]) %Set up color bar
    view(2)     %View from top
    frame = getframe;
    writeVideo(v,frame);

    %Save current fields as old fields for next iteration
    prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
    
    % Save values in measurement points
%     E_save(:,i) = reshape(prev.Ez(indx,indy),numel(xs)*numel(ys),1);

    %Save E-field in boxes
    E_box1(:,:,i) = prev.Ez(boxind(1,3):boxind(1,4),boxind(1,1):boxind(1,2));
    E_box2(:,:,i) = prev.Ez(boxind(2,3):boxind(2,4),boxind(2,1):boxind(2,2));
    E_box3(:,:,i) = prev.Ez(boxind(3,3):boxind(3,4),boxind(3,1):boxind(3,2));
    E_box4(:,:,i) = prev.Ez(boxind(4,3):boxind(4,4),boxind(4,1):boxind(4,2));

end
%% Free videofile
close(gcf)
close(v)

[M,N] = size(E_box1(:,:,1));
[X,Y] = meshgrid(linspace(0,conf.x_length,N),...
                                            linspace(0,conf.y_length,M));

%% Make videos
 
% NOTE: give time interval as parameter(based on existing videos)
plotBox('Box1',E_box1,conf,);
plotBox('Box2',E_box2,conf,);
plotBox('Box3',E_box3,conf,);
plotBox('Box4',E_box4,conf,);

%% Plot graphs
center1=squeeze(E_box1(round(size(E_box1,1)/2),round(size(E_box1,2)/2),:));
center2=squeeze(E_box2(round(size(E_box2,1)/2),round(size(E_box2,2)/2),:));
center3=squeeze(E_box3(round(size(E_box3,1)/2),round(size(E_box3,2)/2),:));
center4=squeeze(E_box4(round(size(E_box4,1)/2),round(size(E_box4,2)/2),:));

figure
subplot(2,2,1)
plot(center1)
subplot(2,2,2)
plot(center2)
subplot(2,2,3)
plot(center3)
subplot(2,2,4)
plot(center4)
%% Release path
rmpath(genpath(pwd))