function [ output_args ] = plotBox( name, box,conf, frames)

v = VideoWriter([pwd '\Vids\FF\' name],'MPEG-4');
v.Quality = 100; 
open(v);
figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])
[M,N] = size(box(:,:,1));
[X,Y] = meshgrid(linspace(0,30,N),...
                                            linspace(0,30,M));
for i=frames(1):frames(2)
 
% Print
    disp(['Frame: ',num2str(i-frames(1)+1),' / ',num2str(frames(2)-frames(1))])
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

