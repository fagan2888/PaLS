% Joshua Enxing
% Tufts University
% Done under supervision of Misha Kilmer and Eric Miller

% Produces surface plot of phi, sum of all basis functions

% Inputs:
%
% Parameter vector |p|
% Meshgrid elements |X| and |Y|
% Epsilon value |eps| that is used with the smooth heaviside function

% Outputs:
%
% Produces a surface plot
function plot_phi(p,X,Y,eps)

Z = zeros(length(X),length(Y));

figure; hold on

for m=1:length(p)/4
    a = p(4*(m-1)+1);
    b = p(4*(m-1)+2);
    c = p(4*(m-1)+3);
    d = p(4*(m-1)+4);
    for i=1:length(X)
        for j=1:length(Y)
            Z(i,j) = phi(a,b,c,d,X(i,j),Y(i,j),eps);
        end
    end
    surf(X,Y,Z);
end

view(17,22);