clc
clear
close all
addpath(genpath(pwd))

% NOTE:
% For comments and info on the general working of the code please check the
% file 'main_commented.m'. The structure is a bit different but the principles
% are the same. Not all general comments are present in this file.

%% Initialisations

%Locations and heights of TX and RX
sloc=2;
Recloc= 8;
source_h=2.5;
rec_h=1.5;


conf.fmax           = 900e6;
conf.x_length       = 10;
conf.y_length       = 10;
conf.nrOfFrames     = 400;
conf.Resolution_X   = 500;
conf.Resolution_Y   = 500;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  
[ field,conf,T ] = FDTDInit( conf );

%% Get indices
Sindx=meter2index(sloc,conf);
Rindx=meter2index(Recloc,conf);
Sindy=meter2index(source_h,conf);
Rindy=meter2index(rec_h,conf);

%% Constants
c       = 3e8;              %Speed of light
mu0     = pi*4e-7;          %Permeabillity of free space 
eps0    = 1/mu0/c/c;        %Permitivity of free space 
Z0      = sqrt(mu0/eps0);   %Free space impedance 

lambda=c/conf.fmax;

%% Range of position building
KE=2.5:0.1:7.5;
%% for loop over location of building
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

%% Initialising the sources
Source = struct;
f = 900e6;
Source = addSource( Source,conf,source_h,sloc,f,sin(2*pi*f*T));

%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];

%% Filling the field with objects
    %Knife-edge
    %Both together deliver good results (just don't do sigma on its own)
    field(1).Sig=draw_rectangle(field(1).Sig,10e10,KE(j),2.5,0.1,5,conf);
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
    
    

    %Save current fields as old fields for next iteration
        prev.Ez=results.Ez;prev.Hx=results.Hx;prev.Hy=results.Hy;
    %Update values for reciever
        recval(j,i)=prev.Ez(Rindy,Rindx);
        soval(j,i)=prev.Ez(Sindy,Sindx);
        if i>=2
        Hx_src(j,i-1)=prev.Hx(Sindy,Sindx);
        Hy_src(j,i-1)=prev.Hy(Sindy,Sindx);
        Hx_rec(j,i-1)=prev.Hx(Rindy,Rindx);
        Hy_rec(j,i-1)=prev.Hy(Rindy,Rindx);
        end
    end
    d1=KE(j)-sloc;
    d2=Recloc-KE(j);
    alpha=atan((source_h-rec_h)/(d1+d2));
    h=5-rec_h-tan(alpha)*d2;
    v(j)=h*sqrt(2/lambda*(1/d1+1/d2));
    
end
H_src = [Hx_src Hy_src zeros(size(Hx_src,1),1)];
E_src = [zeros(size(soval,1),1) zeros(size(soval,1),1) soval];
H_rec = [Hx_rec Hy_rec zeros(size(Hx_rec,1),1)];
E_rec = [zeros(size(recval,1),1) zeros(size(recval,1),1) recval];
P_src = abs(sqrt(H_src(:,1).^2+H_src(:,2).^2)).*abs(E_src(:,3));
P_rec = abs(sqrt(H_rec(:,1).^2+H_rec(:,2).^2)).*abs(E_rec(:,3));

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

Fv=1./(2*pi^2*vth.^2); %Includes already squaring and abs
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
hold on
plot(P_rec./P_src)

%% Free path
rmpath(genpath(pwd))