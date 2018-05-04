clc
clear
close all
addpath(genpath(pwd))

% NOTE:
% For comments and info on the general working of the code please check the
% file 'main_commented.m'. The structure is a bit different but the principles
% are the same. Not all general comments are present in this file.

%% Plan/Reasoning
%
% For normal total reflection you just have to check if the amplitude at
% the reflective object after reflection is different from the amplitude 
% without reflective object. Also one can check the value of the E-field in
% the reflective object.
%
% For non-normal reflection the amplitudes around the point of reflection
% are observed and compared to simulations without reflective object.
% Again the values inside the object of reflection can be checked.

% To generate an accompanying video, please refer to 'Reflection_Visual'


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
value=zeros(1,numel(T));
value(1)=1;
Source = addSource( Source,conf,sourceY,0.5,f,value);%sin(2*pi*f*T)
indY = meter2index(0.5,conf); %Check is always performed at this Y index.
indX= meter2index(0.65,conf)-1; %0.65=0.7-0.1/2

%% Initializing saves
pre=[];
pre2=[];
pre3=[];

post=[];
post2=[];
post3=[];

inside=[];
E_temp2=[];
for loop=1:2
    %% Filling the field with objects

    %Full reflection:
    if loop==2
     field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,0.7,sourceY,0.1,3,conf);
    end

    % Simulating losses
    field(1).SigM(:) = 0;       % No magnetic conductivity in the air
    field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];

    % For more conveniet references
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
    
    %Booleans to steer the checking
    checking=1;
    next=0;
    nextCheck=0;
    nextnextCheck=0;
    nextnextnextCheck=0;
    
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
        prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
% Comment and uncomment the appropriate part for which calculation you
% want.
%     Perpendicular
        if checking && abs(prev.Ez(indY,indX))>0.001
           next=1;
           checking=0;
           pre(loop)=prev.Ez(indY,indX);
           preH(loop)=prev.Hx(indY,indX);
           save1(:,:,loop)= prev.Ez;
        elseif next %Due to the index where we check we have to wait 1 extra cycle
           next=0;
           nextCheck=1;
           pre2(loop)=prev.Ez(indY,indX-1);
           post2(loop)=prev.Ez(indY,indX);
           inside(1,loop) = prev.Ez(indY,indX+1);
           save2(:,:,loop)= prev.Ez;
        elseif nextCheck
           post(loop)=prev.Ez(indY,indX);
           postH(loop)=prev.Hx(indY,indX);
           nextCheck=0;
           inside(2,loop)= prev.Ez(indY,indX+1);
           save3(:,:,loop)= prev.Ez;
        end
        

%     %Angle
%         if checking && abs(prev.Ez(indY,indX))>0.001
%             checking=0;
%             next=1;
%         elseif next
%            next=0
%            nextCheck=1;
%            pre(loop)=prev.Ez(indY,indX);
%            preH(loop)=prev.Hx(indY,indX);
%            pre2(loop)=prev.Ez(indY+1,indX);
%            inside(1,loop) = prev.Ez(indY,indX+1);
%         elseif nextCheck
%            post(loop)=prev.Ez(indY,indX);
%            postH(loop)=prev.Hx(indY,indX);
%            pre3(loop)=prev.Ez(indY+2,indX);
%            nextCheck=0;
%            nextnextCheck=1;
%            inside(2,loop)= prev.Ez(indY,indX+1);
%         elseif nextnextCheck
%            post2(loop)=prev.Ez(indY+1,indX);
%            nextnextCheck=0;
%            nextnextnextCheck=1;
%            inside(3,loop)= prev.Ez(indY,indX+1);
%         elseif nextnextnextCheck
%            post3(loop)=prev.Ez(indY+2,indX);
%            nextnextnextCheck=0;
%            inside(4,loop)= prev.Ez(indY,indX+1);
%         end

        %Save current fields as old fields for next iteration
    end
    prev.Ez
end

%% Reflection coeff (for normal incidence only!)
%n=sqrt(epsrel*murel)
n1 = sqrt(1*1);
n2 = sqrt(1*50000);
R = ((n1-n2)/(n1+n2))^2;
T = 1-R;

%% Show results
disp(['Values before reflection for both simulations at pixel of approach: ' num2str(pre)])
disp(['Values before reflection for both simulations at under angle: ' num2str(pre2) '. (Only valid for approach under angel)'])
disp(['Values before reflection for both simulations at under other angle: ' num2str(pre3) '. (Only valid for approach under angel)'])

disp(['Values after reflection for both simulations at pixel of approach: ' num2str(post)])
disp(['Values after reflection for both simulations at under angle: ' num2str(post2) '. (Only valid for approach under angel)'])
disp(['Values after reflection for both simulations at under other angle: ' num2str(post3) '. (Only valid for approach under angel)'])

disp(['E-field at index of object when it is not there: ' num2str(inside(1,:))])
disp(['E-field inside object at point of relection: ' num2str(inside(2,:))])

disp(['Reflection coefficient = ' num2str(R) ' Transmission coeff = ' num2str(T)])
disp(['Reflected = ' num2str((post(1)-post(2))*2/R) ' Transmitted = ' num2str(inside(2,2))])
disp(['Ratio transmitted compared to ingoing = '  num2str(inside(2,2)/pre(1))])
disp(['Ratio reflected compared to ingoing = '  num2str((post(1)-post(2))/pre(1)*2)])
%% Free path
rmpath(genpath(pwd))