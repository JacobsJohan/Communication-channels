clc
clear
close all
addpath(genpath(pwd))

srcLoc=2;
recLoc= 8;
srcH=2.5;
recH=1.5;

conf.fmax           = 900e6;
conf.x_length       = 10;
conf.y_length       = 10;
conf.nrOfFrames     = 400;
conf.Resolution_X   = 500;
conf.Resolution_Y   = 500;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

% Get indices
Sindx=meter2index(srcLoc,conf);
Rindx=meter2index(recLoc,conf);
Sindy=meter2index(srcH,conf);
Rindy=meter2index(recH,conf);

c       = 3e8;              %Speed of light
mu0     = pi*4e-7;          %Permeabillity of free space 
eps0    = 1/mu0/c/c;        %Permitivity of free space 
Z0      = sqrt(mu0/eps0);   %Free space impedance 

lambda=c/conf.fmax;

KE=2.5:0.1:7.5;
%% for loop over location of 
for j=1:numel(KE)
    disp(['Running KE ' num2str(j) '/' num2str(numel(KE))])
%% Initialising the fields for simulation 
    conf.fmax           = 900e6;
    conf.x_length       = 10;
    conf.y_length       = 10;
    conf.nrOfFrames     = 500;
    conf.Resolution_X   = 500;
    conf.Resolution_Y   = 500;
    conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
    [ field,conf,T ] = FDTDInit( conf );


    deltat  = conf.deltat;
    delta   = conf.delta;
    Sc      = c*deltat/delta;  	%Courant number

%% Initialising the different sources 
    Source = struct;

    f = 900e6;
    Source = addSource( Source,conf,srcH,srcLoc,f,sin(2*pi*f*T));

%% Simulating losses
    field(1).SigM(:) = 0;
    field(1).Sig(:) = 8e-15;

%% Filling the field with objects
    %Knife-edge
    field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,KE(j),2.5,0.1,5,conf);%0.1-0.5m is a good thickness for the metal plate. If you want thinner you need higher epsrel I think.

%% Simulate field

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
    
    

    %Save current fields as old fields for net iteration
        prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
    %Update values for reciever
        recval(j,i)=prev.Ez(Rindy,Rindx);
        soval(j,i)=prev.Ez(Sindy,Sindx);
    end
    d1=KE(j)-srcLoc;
    d2=recLoc-KE(j);
    alpha=atan((srcH-recH)/(d1+d2));
    h=5-recH-tan(alpha)*d2;
    v(j)=h*sqrt(2/lambda*(1/d1+1/d2));
    
end

%Sort results
rec = max(abs(recval),[],2);
src = max(abs(soval),[],2);
ratio=rec./src;
v1=[v' rec];
v2=[v' ratio];

[~,idx] = sort(v1(:,1));
v1 = v1(idx,:);
[~,idx] = sort(v2(:,1));
v2 = v2(idx,:);

%Theory
vth=5:0.1:13;

Fv=1./(2*pi^2*vth.^2);%Already squared and absd
Le=10*log10(Fv);

%Plot results
figure
plot(v2(:,1),10*log10(v2(:,2)))
title('Ratio of max E-fields in reciever and source i.f.o. \nu')
xlabel('\nu')
ylabel('max(E_{rec})/max(E_{src}) (dB)')
hold on 
plot(vth,Le)
legend('Simulation result','Theory')

rmpath(genpath(pwd))