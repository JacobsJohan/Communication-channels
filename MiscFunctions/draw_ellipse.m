function [A] = draw_ellipse(A,value,center_X,center_Y,radius_X,radius_Y,conf)
%draw_ellipse Draw an ellipse on a matrix in meters.
%   Use: draw_ellipse( matrix to draw on, value of ellipse, center x,
%   center y, radius x, radius y, configuration structure) Uses values in
%   'meters'.
[Xnum, Ynum] = size(A);
[X, Y]= meshgrid(linspace(0,conf.x_length,Xnum),linspace(0,conf.y_length,Ynum));
B = (((X-center_X).^2/(radius_X^2)+(Y-center_Y).^2/(radius_Y^2)<=1))*value;
A = A + B - A.*(B~=0);
end

