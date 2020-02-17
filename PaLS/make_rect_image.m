% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Makes a white rectangle on a black image

% Inputs:
% 
% Image size |img_size| makes images img_size by img_size pixels
% X center |x_cent|
% Y center |y_cent|
% Height |height_rect|
% Width |width_rect|

% Outputs:
% 
% Black and white image |rect|
function rect = make_rect_image(img_size,x_cent,y_cent,height_rect,width_rect)

if mod(height_rect,2)==0 || mod(width_rect,2)==0
    disp("Error: height and width must be odd number of pixels");
else 
    rect = zeros(img_size);
    
    left = x_cent - (width_rect-1)/2;
    right = x_cent + (width_rect-1)/2;
    
    top = y_cent - (height_rect-1)/2;
    bottom = y_cent + (height_rect-1)/2;
    
    rect(top:bottom,left:right) = ones(height_rect,width_rect);
end


