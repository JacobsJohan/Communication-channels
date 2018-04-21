clc
clear
close all
addpath(genpath(pwd))
%% Description
% This code is used to model the effect of shadowing. Some trees and
% buildings are placed inside a 20m x 20m frame and the average electric 
% field amplitude at a constant distance (4.5 meters in this case) from the
% central source is measured.

%% Initialising the configuration and fields for simulation 
conf.fmax           = 900e6;
conf.x_length       = 20;
conf.y_length       = 20;
conf.nrOfFrames     = 500;
conf.Resolution_X   = 300;
conf.Resolution_Y   = 300;
conf.Resolution_T   = 500;
conf.ToPrint        = 'Ez';

[ field,conf,T ] = FDTDInit( conf );


%% Initialising the different sources 
Source = struct;
% use:  addSource( Sources,conf,X_loc,Y_loc,frequency,value ) 
%       X_loc/Y_loc/Value can either be a scalar or vector of length
%       conf.nrOfFrames

f = 900e6;
Source = addSource( Source,conf,10,10,f,sin(2*pi*f*T) );

%% Simulating losses
field(1).SigM(:) = 0;       % No magnetic conductivity in the air
field(1).Sig(:) = 8e-15;    % conductivity of air is [3e-15, 8e-15];

%% Filling the field with objects
% draw_rectangle(A,value,center_X,center_Y,width,height,conf)
% draw_ellipse(A,value,center_X,center_Y,radius_X,radius_Y,conf)

% Example: Place some trees with relative permittivity of wood between 2 and 6.
%          Place some buildings with relative permittivity of 4.5 (concrete).

field(1).EpsRel = draw_rectangle(field(1).EpsRel, 4.5, 14, 10, 1, 6, conf);
field(1).EpsRel = draw_rectangle(field(1).EpsRel, 4.5, 10, 6, 6, 1, conf);

field(1).EpsRel = draw_ellipse(field(1).EpsRel, 2, 8, 16, 1, 1, conf);
field(1).EpsRel = draw_ellipse(field(1).EpsRel, 2, 4, 10, 1, 1, conf);

% The relative permeabilities for wood and concrete are very close to 1

% The electric conductivity for wood = 1e-15, for concrete = 0
% For the sake of completeness, these losses are added, but they don't make
% any difference since the conductivity of air is also ~0

field(1).Sig = draw_rectangle(field(1).Sig, 0, 7, 5, 1, 6, conf);
field(1).Sig = draw_rectangle(field(1).Sig, 0, 4, 2, 4, 1, conf);

field(1).Sig = draw_ellipse(field(1).Sig, 1e-15, 4, 8, 1, 1, conf);
field(1).Sig = draw_ellipse(field(1).Sig, 1e-15, 2, 4, 1, 1, conf);

%% Simulate field
[field] = FDTDMaxwellCore( field,conf,Source );

%% Calculate average electric field amplitude at a constant distance from the source
% In this case at 4.5 m from the source (circle with radius 9 m)
% Locations is a matrix that represents the circle in which we will study
% the shadowing.

Locations = draw_ellipse(field(1).Ez,10,10,10,4.5,4.5,conf);
Locations = draw_ellipse(Locations,1,10,10,4.46,4.46,conf);
% Locations = draw_ellipse(Locations,1,10,10,4.3,4.3,conf);     % For thicker circle

[row, col] = find(Locations == 10);

temp = reshape([field.Ez],length(Locations),length(Locations),[]);

% Calculate the mean value of the electric field over all time frames
tempMean = mean(temp,3);        

% AvgField will be a matrix containing the average electric field values at
% all Locations we are interested in. The other values of the matrix are 0.
AvgField=zeros(size(field(1).Ez));

% CAN BE DELETED
% AvgFieldVals will be a vector containing all electric field values we are
% interested in, such that we can make a shadowing plot. Note that due to
% the way the function 'find' works, the values in AvgFieldVals will not
% represent contiguous values of the circle of interest. This will be dealt
% with in AvgFieldVals2
AvgFieldVals = 0;

% Save minimal value of average electric field
minField = 100;
for i=1:length(row)
    AvgField(col(i),row(i)) = tempMean(col(i),row(i));
    AvgFieldVals = [AvgFieldVals, tempMean(col(i),row(i))];
    
    if (tempMean(col(i),row(i)) < minField)
        minField = tempMean(col(i),row(i));
    end
end
AvgFieldVals = AvgFieldVals(2:end);                 % Remove initial 0
[minRow, minCol] = find(AvgField == minField);

% Add minVal to AvgField to plot later
% AvgField = draw_ellipse(AvgField, 3e-3, minCol*conf.delta, minRow*conf.delta, 0.5, 0.5, conf);
tempMean = draw_ellipse(tempMean, 3e-3, minCol*conf.delta, minRow*conf.delta, 0.5, 0.5, conf);

%% Make vector with contiguous average electric field values.
% Instead of using an advanced pathfinding algorithm to find electric field
% values of contiguous locations on the circle of interest, a more simple
% approach will be used. We will run over the circle one quadrant at a time
% in the following way: quadrant II - I - IV - III

% Define starting values
j = 300;
i = 166;
AvgFieldVals2 = 0;
% Compute average field in first quarter circle
while i<301
    if (AvgField(i,j) ~= 0)
        AvgFieldVals2 = [AvgFieldVals2, AvgField(i,j)];
    end
    % Check where the next contiguous value will be
    if (AvgField(i,j-1) ~= 0)
        j = j-1;
    elseif (AvgField(i+1,j) ~= 0)
        i = i+1;
    elseif (AvgField(i+1,j-1) ~= 0)
        i = i+1;
        j = j-1;
    end
end
% Compute average field in 2nd quarter circle
while j<301
    if (AvgField(i,j) ~= 0)
        AvgFieldVals2 = [AvgFieldVals2, AvgField(i,j)];
    end
    % Check where the next contiguous value will be
    if (AvgField(i+1,j) ~= 0)
        i = i+1;
    elseif (AvgField(i,j+1) ~= 0)
        j = j+1;
    elseif (AvgField(i+1,j+1) ~= 0)
        i = i+1;
        j = j+1;
    end
end
% Compute average field in 3rd quarter circle
while i>300
    if (AvgField(i,j) ~= 0)
        AvgFieldVals2 = [AvgFieldVals2, AvgField(i,j)];
    end
    % Check where the next contiguous value will be
    if (AvgField(i,j+1) ~= 0)
        j = j+1;
    elseif (AvgField(i-1,j) ~= 0)
        i = i-1;
    elseif (AvgField(i-1,j+1) ~= 0)
        i = i-1;
        j = j+1;
    end
end
% Compute average field in 4th quarter circle
while j>301
    if (AvgField(i,j) ~= 0)
        AvgFieldVals2 = [AvgFieldVals2, AvgField(i,j)];
    end
    % Check where the next contiguous value will be
    if (AvgField(i-1,j) ~= 0)
        i = i-1;
    elseif (AvgField(i,j-1) ~= 0)
        j = j-1;
    elseif (AvgField(i-1,j-1) ~= 0)
        i = i-1;
        j = j-1;
    end
end

AvgFieldVals2 = AvgFieldVals2(2:end); % Remove initial zero

minIndex = find(AvgFieldVals2 == minField);

%% Plots
% The first plot shows the locations of the buildings and trees we placed.
figure
imagesc(field(1).EpsRel)
title('Figure showing locations of placed objects')
axis off

% The next plot is a shadowing plot of non contiguous electric field
% values, hence it is probably not as intersting
% figure
% plot(AvgFieldVals)
% xlim([1,length(AvgFieldVals)])
% title('Amplitude of electrical field at constant distance from the source')
% ylabel('E [V/m]')

% The next plot is a shadowing plot for contiguous electric field values.
% The minimum electric field value is shown by a red star.

figure
plot(AvgFieldVals2)
hold on
plot(minIndex,minField,'*')
xlim([1,length(AvgFieldVals2)])
title('Amplitude of electrical field at constant distance from the source')
ylabel('E [V/m]')


% The next plot shows an image of the locations where we study the
% shadowing. The big dot show the minimal value of the shadowing which
% corresponds to the red star in the previous plot.
figure
imagesc(AvgField)
colorbar
title('Amplitude of electrical field at constant distance from the source')
axis off

% The next plot shows an image of the average electric field over the
% entire area. The location of the buildings can clearly be distinguished.
% The yellow dot at the lower left side of the image corresponds to the
% same minimal location as before.
figure
imagesc(tempMean)
colorbar
axis off

%%
rmpath(genpath(pwd))