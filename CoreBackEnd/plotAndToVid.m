function [] = plotAndToVid(filename,field,conf,plotMax,plotMin)
% Functions plots field with the characteristics given by conf and writes
% frames to video with name 'filename'.
%% Resize fields 

for i = 1:numel(field)
    V(:,:,i)        = field(i).(conf.ToPrint);
    EpsRel(:,:,i)   = field(i).EpsRel;
    MuRel(:,:,i)    = field(i).MuRel;
end

[M,N,T] = size(V);
[X,Y,Z]         = meshgrid(     linspace(0,conf.x_length,N),...
                                linspace(0,conf.y_length,M),...
                                linspace(0,conf.deltat*T,T));
[Xq,Yq,Zq]      = meshgrid(     linspace(0,conf.x_length,conf.Resolution_X),...
                                linspace(0,conf.y_length,conf.Resolution_Y),...
                                linspace(0,conf.deltat*T,conf.Resolution_T));                 
ToPrintq        = interp3(X,Y,Z,V,Xq,Yq,Zq);

[M,N,T] = size(EpsRel);
[X2,Y2,Z2]         = meshgrid(      linspace(0,conf.x_length,M),...
                                    linspace(0,conf.y_length,N),...
                                    linspace(0,conf.deltat*T,T));
[Xq2,Yq2,Zq2]       = meshgrid(     linspace(0,conf.x_length,conf.Resolution_X),...
                                    linspace(0,conf.y_length,conf.Resolution_Y),...
                                    linspace(0,conf.deltat*T,conf.Resolution_T)); 
% Interpolate spatial characteristics                 
EPSrelalpha     = (interp3(X2,Y2,Z2,EpsRel,Xq2,Yq2,Zq2) -1) /   (max(EpsRel(:))-1)/2.5;
MUrelalpha      = (interp3(X2,Y2,Z2,MuRel,Xq2,Yq2,Zq2)  -1) /   (max(MuRel(:))-1)/2.5;

%Select field to print
temp = ToPrintq(:,:,20:end);
maxToPrint = max(temp(:));
minToPrint = min(temp(:));
absMaxToPrint = max(ToPrintq(:));


%% Print
%Prepare video
v = VideoWriter(filename,'MPEG-4');
v.Quality = 100; 
open(v);
figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])

% Print frames and write them to video
for i = 1:conf.Resolution_T
    disp(['Frame: ',num2str(i),' / ',num2str(conf.Resolution_T)])
    surf(Xq(:,:,i),Yq(:,:,i),ToPrintq(:,:,i),...
            'LineStyle','none',...
            'FaceColor','interp');
    hold on 
    surf(   Xq2(:,:,i),...
            Yq2(:,:,i),...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',EPSrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','red');
    surf(   Xq2(:,:,i),...
            Yq2(:,:,i),...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint+0.1,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',MUrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','blue');
    text(X2(end,end)/10,Y2(end,end)/10,absMaxToPrint+0.2,['time = ',num2str(Zq(1,1,i)),'s']);
    hold off
    colorbar;
    zlim([minToPrint,absMaxToPrint+0.2]);
    caxis([plotMin,plotMax])
    view(2)
    frame = getframe;
    %img = print('-RGBImage');
    %frame = im2frame(img);
    writeVideo(v,frame);
end
% Free video
close(gcf)
close(v)

end

