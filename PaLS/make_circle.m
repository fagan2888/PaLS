% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Makes a white circle on a black background

% Inputs:
%
% Width of image |x_length|
% Height of image |y_length|
% X center of circle in image |x_c|
% Y center of circle in image |y_c|
% Radius of circle |r|

% Outputs:
%
% Matrix |circlePixels|
function circlePixels = make_circle(x_length,y_length,x_c,y_c,r)

% Create a logical image of a circle with specified
% diameter, center, and image size.

% First create the image.
imageSizeX = x_length;
imageSizeY = y_length;
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

% Next create the circle in the image.
centerX = x_c;
centerY = y_c;
radius = r;
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;