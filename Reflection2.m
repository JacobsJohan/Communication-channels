clc
clear
close all
addpath(genpath(pwd))

%PROBLEM: De amplitude wijzigt, er zitten erges losses in denkik... Ale als
%ge de amplitude van 2 punte naast elkaar neemt dan is de ene ni de waarde
%van het andere vanuit de vorige frame wa ik wel zou verwachten... Iemand
%een oplossing?

%% Normal total reflection
%
%For normal total reflection you just have to check if the amplitude after
%moment of reflection is the sum of the incident amplitude and the one
%reflected(incident one frame before or 2?)
%
%For normal non-reflection you just have to check if the amplitude after
%moment of reflection is the sum of the incident amplitude and the one
%reflected(incident one frame before or 2?) and the amplitude right behind
%the object (very thin)


%% Non normal reflection
%
%For non normal total reflection you just have to check if the amplitude after
%moment of reflection is the sum of the incident amplitude and the one
%reflected(incident one frame before or 2?) at th epoint corresponding to
%the right angle
%
%For non normal non-reflection you just have to check if the amplitude after
%moment of reflection is the sum of the incident amplitude and the one
%reflected(incident one frame before or 2?) corresponding to the right
%angle and the amplitude right behind the object (very thin)(in the
%direction of the incident wave)

%% Initialising the fields for simulation 
conf.fmax           = 900e6;
conf.x_length       = 1;
conf.y_length       = 1;
conf.nrOfFrames     = 30;
conf.Resolution_X   = 500;
conf.Resolution_Y   = 500;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf);

%% Initialising the different sources 
Source = struct;

f = 900e6;
sourceY = 0.5; %For perpendicular
% sourceY = 0.2;      %For angle
Source = addSource( Source,conf,sourceY,0.5,f,sin(2*pi*f*T));
% indY = meter2index(sourceY,conf);
indY = meter2index(0.5,conf);
indX= meter2index(0.65,conf)-1; %0.65=0.7-0.1/2

for loop=1:2
    %% Filling the field with objects

    %Full reflection normal:
    if loop==2
     field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,0.7,sourceY,0.1,3,conf);
    end
    %Partial reflection normal
    %  field(1).EpsRel = draw_rectangle(field(1).EpsRel,50,1.5,1.5,0.01,3,conf);

    %Full reflection non normal
    %  field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,1.5,2,0.1,3,conf);

    %Partial reflection non normal
    %  field(1).EpsRel = draw_rectangle(field(1).EpsRel,50,1.5,2,0.01,3,conf);


    % Simulating losses
    field(1).SigM(:) = 300;
    field(1).Sig(:) = 300;

    %% Prep video
    EpsRel = field(1).EpsRel();
    MuRel = field(1).MuRel();



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

    magnitudebefore=0;
    magnitudeafter=0;
    
    %Booleans
    checking=1;
    next=0;
    nextCheck=0;
    nextnextCheck=0;
    nextnextnextCheck=0;
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
    %Perpendicular
%         if i==8
%            prev.Ez
%            pre(loop)=prev.Ez(indY,indX);
%            preH(loop)=prev.Hx(indY,indX);
%         elseif i==9
%            post(loop)=prev.Ez(indY,indX);
%            postH(loop)=prev.Hx(indY,indX);
%         end
% 
% %Angle
%         if i==15
%             prev.Ez
%            pre(loop)=prev.Ez(indY,indX);
%            preH(loop)=prev.Hx(indY,indX);
%            pre2(loop)=prev.Ez(indY+1,indX);
%         elseif i==16
%            post(loop)=prev.Ez(indY,indX);
%            postH(loop)=prev.Hx(indY,indX);
%            pre3(loop)=prev.Ez(indY+2,indX);
%         elseif i==17
%            post2(loop)=prev.Ez(indY+1,indX);
%         elseif i==18
%            post3(loop)=prev.Ez(indY+2,indX);
%         end


% Without indices:
%     Perpendicular
%         if checking && abs(prev.Ez(indY,indX))>0.001
%            next=1;
%            checking=0;
%         elseif next %Due to the index where we check we have to wait 1 extra cycle
%            pre(loop)=prev.Ez(indY,indX);
%            preH(loop)=prev.Hx(indY,indX);
%            next=0;
%            nextCheck=1;
%         elseif nextCheck
%            post(loop)=prev.Ez(indY,indX);
%            postH(loop)=prev.Hx(indY,indX);
%            nextCheck=0;
%         end
% 

    %Angle
        if checking && abs(prev.Ez(indY,indX))>0.001
            checking=0;
            next=1;
        elseif next
           next=0
           nextCheck=1;
           pre(loop)=prev.Ez(indY,indX);
           preH(loop)=prev.Hx(indY,indX);
           pre2(loop)=prev.Ez(indY+1,indX);
        elseif nextCheck
           post(loop)=prev.Ez(indY,indX);
           postH(loop)=prev.Hx(indY,indX);
           pre3(loop)=prev.Ez(indY+2,indX);
           nextCheck=0;
           nextnextCheck=1;
        elseif nextnextCheck
           post2(loop)=prev.Ez(indY+1,indX);
           nextnextCheck=0;
           nextnextnextCheck=1;
        elseif nextnextnextCheck
           post3(loop)=prev.Ez(indY+2,indX);
           nextnextnextCheck=0;
        end




        %Save current fields as old fields for net iteration
            prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
    end
end
rmpath(genpath(pwd))