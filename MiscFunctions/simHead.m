function [ EPSrel, SIG ] = simHead( EPSrel, SIG, delta, xLoc, yLoc, conf )
% Get relative permitivity of head
head = imread('headBlue.jpg');
head = head(:,:,1);
% Relative permittivities and conductivities at 900 MHz
for i=1:length(head(:,1))
    for j=1:length(head(1,:))
        if (head(i,j) >= 200)
            eps(i,j) = 1; %Air
            sig(i,j) = 0;
        elseif (head(i,j) >= 170)
            eps(i,j) = 41; %Skin
            sig(i,j) = 0.9;
        elseif (head(i,j) >= 150)
            eps(i,j) = 5; %Fat
            sig(i,j) = 0.05;
        elseif (head(i,j) >= 130)
%             eps(i,j) = 21; %Skull: Cancellous bones
            eps(i,j) = 12; %Skull (bone cortical)
            sig(i,j) = 0.16;
        elseif (head(i,j) >= 110)
            eps(i,j) = 69; %Cerebrospinal fluid
            sig(i,j) = 2.5;
        elseif (head(i,j) >= 90)
            eps(i,j) = 53; %Gray matter
            sig(i,j) = 1;
        elseif (head(i,j) >= 45)
            eps(i,j) = 39; %White matter
            sig(i,j) = 0.62;
        else
            eps(i,j) = 1;
            sig(i,j) = 0;
        end
    end
end

%% Rescale (Average human head circumference is 55 cm)
% Choose 20 cm in y direction, 16 in x direction
% delta = lambda/10 = 0.03m for 900 MHz
headWidth = 0.16;
headLength = 0.20;
% The head from the figure is 357x276 pixels. Assume this maps for a
% 0.2x0.16 m head.
[X, Y] = meshgrid(linspace(0,headWidth,length(head(1,:))), ...
                  linspace(0,headLength,length(head(:,1))));
% This should be transformed to a xsize by ysize head
xsize = round(headWidth/delta)*2; % x2 for eps and mu?
ysize = round(headLength/delta)*2;

[Xq, Yq] = meshgrid(linspace(0,headWidth,xsize), ...
                    linspace(0,headLength,ysize));

eps_small = interp2(X,Y,eps,Xq,Yq);
sig_small = interp2(X,Y,sig,Xq,Yq);

xIndex = meter2index(xLoc-headWidth/2, conf)*2; % for EpsRel and MuRel multiply result by 2
yIndex = meter2index(yLoc-headLength/2, conf)*2; % -headLength/2 for centered coordinates


for i=1:ysize
    for j=1:xsize
        EPSrel(yIndex+i,xIndex+j) = eps_small(i,j);
        SIG(yIndex+i,xIndex+j) = sig_small(i,j);
    end
end

end

