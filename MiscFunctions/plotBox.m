function [ output_args ] = plotBox( name, box,conf)

v = VideoWriter([pwd '\Vids\FF\' name],'MPEG-4');
v.FrameRate = 8;
v.Quality = 100; 
open(v);
figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])
[M,N] = size(box(:,:,1));
[X,Y] = meshgrid(linspace(0,conf.x_length,N),...
                                            linspace(0,conf.y_length,M));
for i=1:conf.nrOfFrames-1
 
% Print
    disp(['Frame: ',num2str(i),' / ',num2str(conf.nrOfFrames)])
    surf(X,Y,box(:,:,i),...
            'LineStyle','none',...
            'FaceColor','interp');
    hold on
    colorbar;
    caxis([-0.5,0.5]) %Set up color bar
    view(2)     %View from top
    frame = getframe;
    writeVideo(v,frame);
end
close(gcf)
close(v)

end

