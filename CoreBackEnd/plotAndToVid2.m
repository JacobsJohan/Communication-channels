function [] = plotAndToVid2(filename,field,conf,plotMax,plotMin)
% Functions plots field with the characteristics given by conf and writes
% frames to video with name 'filename'.

%Note this is the same as plotAndToVid but without interpolation and for 1
%frame. Not used in practice.
% Prepare for plotting

EpsRel = field(1).EpsRel();
MuRel = field(1).MuRel();

[M,N] = size(results.(conf.ToPrint));
[X,Y] = meshgrid(linspace(0,conf.x_length,N), linspace(0,conf.y_length,M));

ToPrintq=result.(conf.ToPrint);

[M,N,T] = size(EpsRel);
[X2,Y2]         = meshgrid(      linspace(0,conf.x_length,M),...
                                    linspace(0,conf.y_length,N)...
                                    );
                    
EPSrelalpha     = (EpsRel -1) /   (max(EpsRel(:))-1)/2.5;
MUrelalpha      = (MuRel -1) /   (max(MuRel(:))-1)/2.5;


temp = ToPrintq(:,:,20:end);
% maxToPrint = max(temp(:));
% minToPrint = min(temp(:));
absMaxToPrint = max(ToPrintq(:));


%% Print
% v = VideoWriter(filename,'MPEG-4');
v = VideoWriter('Output','MPEG-4');
v.Quality = 100; 
open(v);
figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])
for i = 1:conf.nrOfFrames
    disp(['Frame: ',num2str(i),' / ',num2str(conf.Resolution_T)])
    surf(X(:,:,i),Y(:,:,i),ToPrintq(:,:,i),...
            'LineStyle','none',...
            'FaceColor','interp');
    hold on 
    surf(   X2(:,:,i),...
            Y2(:,:,i),...
            ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint,...
            'FaceAlpha','interp',...
            'AlphaDataMapping','none',...
            'AlphaData',EPSrelalpha(:,:,1),...
            'LineStyle','none',...
            'FaceColor','red');
    surf(   X2(:,:,i),...
            Y2(:,:,i),...
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
    writeVideo(v,frame);
end
close(gcf)
close(v)

end

