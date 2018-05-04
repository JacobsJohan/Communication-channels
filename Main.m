clc
clear
close all
addpath(genpath(pwd))

%Better to use main_commented

%% Initialising the fields for simulation 
conf.fmax           = 900e6;
conf.x_length       = 3;
conf.y_length       = 3;
conf.nrOfFrames     = 200;
conf.Resolution_X   = 200;
conf.Resolution_Y   = 200;
conf.Resolution_T   = 300;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  

% delta = 3e8./conf.fmax/10/5; % For SAR study with head higher resolution is needed
[ field,conf,T ] = FDTDInit( conf );


%% Initialising the different sources 
Source = struct;
% use:  addSource( Sources,conf,X_loc,Y_loc,frequency,value ) 
%       X_loc/Y_loc/Value can either be a scalar or vector of length
%       conf.numberOfFrames

f = 900e6;
% Source = addSource( Source,conf,0.5,0.5,f,sin(2*pi*f*T) );

% Source = addSource( Source,conf,1+1*cos(2*pi*f/16*T),1+1*sin(2*pi*f/16*T),f,sin(2*pi*f*T) );
% Source = addSource( Source,conf,2,2,f,sin(2*pi*f*T) );%Voor 10x10
% Source = addSource( Source,conf,0.33,0.33,f,sin(2*pi*f*T) );%voor 0.66x0.66


Source = addSource( Source,conf,1.5,1,f,(1-2*(pi*f*T).^2).*exp(-(pi*f*T).^2) );

%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % Conductivity of air is [3e-15, 8e-15];


%% Filling the field with objects
% field(1).EpsRel = draw_rectangle(field(1).EpsRel,4,2,2,1,1,conf);
% field(1).EpsRel = draw_ellipse(field(1).EpsRel,8,2,2,0.5,0.25,conf);

% field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,1.5,1,3,.2,conf);
% field(1).EpsRel = draw_rectangle(field(1).EpsRel,50000,1.5,2,3,.2,conf);

% field(1).EpsRel = simHead(field(1).EpsRel, field(1).Sig, conf.delta, 0.6, 0.5, conf);

%% SAR study
% SAR = sigma*E^2/rho
% E = rms values of electric field amplitude
% sigma = electric conductivity
% rho = tissue mass density

%% Simulate field
[field] = FDTDMaxwellCore( field,conf,Source );

%% Draw and export to movie 
plotAndToVid('Output/simulation',field,conf,0.5,-0.5)
rmpath(genpath(pwd))