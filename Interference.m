clc
clear
close all
addpath(genpath(pwd))


%%Load data from one source
load('OneSource')
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
% Base source
Source = addSource( Source,conf,1,1.5,f,sin(2*pi*f*T) );
% Source in phase
Source = addSource( Source,conf,2,1.5,f,sin(2*pi*f*T) );
% Source i anti-phase
% Source = addSource( Source,conf,2,1.5,f,sin(2*pi*f*T+pi));
%% Simulating losses
field(1).SigM(:) = 300;
field(1).Sig(:) = 300;

%% Prep video
% v = VideoWriter('Output_Interference','MPEG-4');
% v.Quality = 100; 
% open(v);
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

E_temp=[]; %Initialize matrix to keep values at middel position between the 2 sources
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
    [X,Y] = meshgrid(linspace(0,conf.x_length,N),...
        linspace(0,conf.y_length,M));

    ToPrintq=results.(conf.ToPrint);
    temp = ToPrintq(:,:,20:end);
    % maxToPrint = max(temp(:));
    minToPrint = min(temp(:));
    absMaxToPrint = max(ToPrintq(:));


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
% %     text(X2(end,end)/10,Y2(end,end)/10,absMaxToPrint+0.2,['time = ',num2str(Zq(1,1,i)),'s']);
%     hold off
%     colorbar;
% %     zlim([minToPrint,absMaxToPrint+0.2]);
%     caxis([-0.5,0.5])
%     view(2)
%     frame = getframe;
%     writeVideo(v,frame);

%Save current fields as old fields for net iteration
    prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
    
%    Saving current field in middle between 2 sources.
   E_temp=[E_temp prev.Ez(meter2index(1.5,conf)+1,meter2index(1.5,conf))];
   
end
%% Free videofile
% close(gcf)
% close(v)

figure
plot(E_temp)
title('E-field between 2 sources in phase')
xlabel('Frame number')
ylabel('E_z-field (V/m)')
figure
plot(E_temp./2)
hold on
plot(OneSource,'*')
legend('Halved interfered E_z-field','E_z -field for one source')
title({'E-field halved between 2 sources in phase' , 'compared to field for 1 source'})
xlabel('Frame number')
ylabel('E_z-field (V/m)')

figure
plot(E_temp(1:65))
hold on
plot(OneSource(1:65))
plot(E_temp(1:65)./2,'bl*')
legend('E_z-field for 2 interfering sources','E_z -field for one source','For 2 sources, halved')
title({'E-field in middle between 2 sources in phase,' , 'compared to field for 1 source'})
xlabel('Frame number')
ylabel('E_z-field (V/m)')


%% Draw and export to movie 
% plotAndToVid2('Output/simulation2',field,conf,0.5,-0.5)
rmpath(genpath(pwd))