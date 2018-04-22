clc
clear
close all
addpath(genpath(pwd))
% NOTE:
% For comments and info on the general working of the code please check the
% file 'main_commented.m'. The structure is a bit different but the principles
% are the same. Not all general comments are present in this file.

%% Initialising the configuration and fields for simulation 
conf.fmax           = 900e6;
conf.x_length       = 10;
conf.y_length       = 10;
conf.nrOfFrames     = 500;
conf.Resolution_X   = 300;
conf.Resolution_Y   = 300;
conf.Resolution_T   = 500;
conf.ToPrint        = 'Ez';  %Needs to be field of the structure 'Field'  

% delta = 3e8./conf.fmax/10/5; % For SAR study with head higher resolution is needed
[ field,conf,T ] = FDTDInit( conf );


%% Initialising the different sources 
Source = struct;
% use:  addSource( Sources,conf,X_loc,Y_loc,frequency,value ) 
%       X_loc/Y_loc/Value can either be a scalar or vector of length
%       conf.nrOfFrames

f = 900e6;
Source = addSource( Source,conf,5,5,f,sin(2*pi*f*T) );

%% Simulating losses
field(1).SigM(:) = 0; % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15; % conductivity of air is [3e-15, 8e-15];

%% Filling the field with objects
% draw_rectangle(A,value,center_X,center_Y,width,height,conf)
% draw_ellipse(A,value,center_X,center_Y,radius_X,radius_Y,conf)

% Example: Place some trees with relative permittivity of wood between 2 and 6.
%          Place some buildings with relative permittivity of 4.5
%          (concrete). It can be a bit higher if you consider glass etc.

field(1).EpsRel = draw_rectangle(field(1).EpsRel, 4.5, 7, 5, 1, 6, conf);
field(1).EpsRel = draw_rectangle(field(1).EpsRel, 4.5, 4, 2, 4, 1, conf);

field(1).EpsRel = draw_ellipse(field(1).EpsRel, 6, 4, 8, 1, 1, conf);
field(1).EpsRel = draw_ellipse(field(1).EpsRel, 6, 2, 4, 1, 1, conf);

% The relative permeabilities for wood and concrete are very close to 1

% The electric conductivity for wood = 1e-15, for concrete = 0
% These 4 lines of code don't make a difference since the conductivity for
% air is also ~0
field(1).Sig = draw_rectangle(field(1).Sig, 0, 7, 5, 1, 6, conf);
field(1).Sig = draw_rectangle(field(1).Sig, 0, 4, 2, 4, 1, conf);

field(1).Sig = draw_ellipse(field(1).Sig, 1e-15, 4, 8, 1, 1, conf);
field(1).Sig = draw_ellipse(field(1).Sig, 1e-15, 2, 4, 1, 1, conf);

%% Simulate field
[field] = FDTDMaxwellCore( field,conf,Source );

%% Calculate average electric field amplitude at a constant distance from the source
% In this case at 4.5 m from the source (circle with radius 9 m)
Locations = draw_ellipse(field(1).Ez,10,5,5,4.5,4.5,conf);
Locations = draw_ellipse(Locations,1,5,5,4.46,4.46,conf);

% Define starting values
j = 150;
i = 16;
AvgField = 0;
% Compute average field in first quarter circle
while i<151
    if (Locations(i,j) == 10)
        temp = 0;
        for k=1:conf.nrOfFrames
            temp = temp + field(k).Ez(i,j);
        end
        AvgField = [AvgField temp/conf.nrOfFrames]; 
    end
    % Check where the next contiguous value will be
    if (Locations(i,j-1) == 10)
        j = j-1;
    elseif (Locations(i+1,j) == 10)
        i = i+1;
    elseif (Locations(i+1,j-1) == 10)
        i = i+1;
        j = j-1;
    end
end
% Compute average field in 2nd quarter circle
while j<151
    if (Locations(i,j) == 10)
        temp = 0;
        for k=1:conf.nrOfFrames
            temp = temp + field(k).Ez(i,j);
        end
        AvgField = [AvgField temp/conf.nrOfFrames]; 
    end
    % Check where the next contiguous value will be
    if (Locations(i+1,j) == 10)
        i = i+1;
    elseif (Locations(i,j+1) == 10)
        j = j+1;
    elseif (Locations(i+1,j+1) == 10)
        i = i+1;
        j = j+1;
    end
end
% Compute average field in 3rd quarter circle
% Evaluation of the resulting plot gives a minimum value for the AvgField
% at index 532. The commented code is to find the correponding indices i
% and j.
while i>150
    if (Locations(i,j) == 10)        
%         if (length(AvgField) == 533)
%             i
%             j
%         end
        temp = 0;
        for k=1:conf.nrOfFrames
            temp = temp + field(k).Ez(i,j);
        end
        AvgField = [AvgField temp/conf.nrOfFrames]; 
    end
    % Check where the next contiguous value will be
    if (Locations(i,j+1) == 10)
        j = j+1;
    elseif (Locations(i-1,j) == 10)
        i = i-1;
    elseif (Locations(i-1,j+1) == 10)
        i = i-1;
        j = j+1;
    end
end
% Compute average field in 4th quarter circle
while j>151
    if (Locations(i,j) == 10)
        temp = 0;
        for k=1:conf.nrOfFrames
            temp = temp + field(k).Ez(i,j);
        end
        AvgField = [AvgField temp/conf.nrOfFrames]; 
    end
    % Check where the next contiguous value will be
    if (Locations(i-1,j) == 10)
        i = i-1;
    elseif (Locations(i,j-1) == 10)
        j = j-1;
    elseif (Locations(i-1,j-1) == 10)
        i = i-1;
        j = j-1;
    end
end

AvgField = AvgField(2:end);


%% Plots
figure
plot(AvgField)
xlim([1,length(AvgField)])
title('Amplitude of electrical field at constant distance from the source')
ylabel('E [V/m]')

%% First draw circle in which we compute electric field amplitude, then draw buildings and trees again (they get overwritten)
Locations2 = draw_ellipse(field(1).EpsRel,10,5,5,4.5,4.5,conf);
field(1).EpsRel = draw_ellipse(Locations2,1,5,5,4.46,4.46,conf);

field(1).EpsRel = draw_rectangle(field(1).EpsRel, 4.5, 7, 5, 1, 6, conf);
field(1).EpsRel = draw_rectangle(field(1).EpsRel, 4.5, 4, 2, 4, 1, conf);

field(1).EpsRel = draw_ellipse(field(1).EpsRel, 6, 4, 8, 1, 1, conf);
field(1).EpsRel = draw_ellipse(field(1).EpsRel, 6, 2, 4, 1, 1, conf);

% Minimum value in shadowing plot is at i = 282, j = 174
% We will show this location in the figure
field(1).EpsRel = draw_ellipse(field(1).EpsRel, 10, 10*282*2/599, 10*174*2/599, .2, .2, conf);
% Draw source in middle as small dot
field(1).EpsRel = draw_ellipse(field(1).EpsRel, 10, 5, 5, .1, .1, conf);

figure
imagesc(field(1).EpsRel)
title('Shadowing at a constant distance from the source (in center)')
axis off
% With the current settings, the minimum value for the electric field
% amplitude can be found at the yellow spot. Note that this is not a tree

rmpath(genpath(pwd))