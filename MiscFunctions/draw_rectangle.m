function [A] = draw_rectangle(A,value,center_X,center_Y,width,height,conf)
%draw_rectangle Draw a rectangle on a matrix in meters.
%   Use: draw_rectangle( matrix to draw on, value of rectangle, center x,
%   center y, width, height, configuration structure) Uses values in
%   'meters'.
[Xnum, Ynum] = size(A);
[X, Y]= meshgrid(linspace(0,conf.x_length,Xnum),linspace(0,conf.y_length,Ynum));
B = (abs(X-center_X)<=width/2 & abs(Y-center_Y)<=height/2)*value;
A = A + B - A.*(B~=0);
end

