function [perc_v1,perc_v2, s] = percent_distance(v1,v2,pts)
%%
%{
author: Evan J. Ramos
date:   23 July 2020

description: this script determines the coordinate location of a data point
by projecting its location perpendicularly on a line. The distance that the
point is between the two vertices that define the line will help determine
the percent similarity the points are to the vertices.

inputs:
        v1: 1 x 2 vector (i.e. [x-coord, y-coord]) of vertex 1
        v2: 1 x 2 vector (i.e. [x-coord, y-coord]) of vertex 2
       pts: N x 2 vector of coordinate locations for points b/w v1 and v2

outputs:
        perc_v1 = N x 1 vector of percent distance from v1
        perc_v2 = N x 1 vector of percent distance from v2
        s       = scalar representing variance of data from best fit line

        100%: closest; 0%: farthest
%}

%% Fit line between v1 and v2, get polynomial coefficients

%normalize distances
max_x = max([v1(1) v2(1)]);
min_x = min([v1(1) v2(1)]);
max_y = max([v1(2) v2(2)]);
min_y = min([v1(2) v2(2)]);
pts(:,1) = (pts(:,1)-min_x)/(max_x-min_x);
pts(:,2) = (pts(:,2)-min_y)/(max_y-min_y);
v1(1) = (v1(1)-min_x)/(max_x-min_x); v1(2) = (v1(2)-min_y)/(max_y-min_y);
v2(1) = (v2(1)-min_x)/(max_x-min_x); v2(2) = (v2(2)-min_y)/(max_y-min_y);

poly_coeffs = polyfit([v1(1) v2(1)],[v1(2) v2(2)],1);
%poly_coeffs(1): slope; poly_coeffs(2): y-intercept

%% Determine perpendicular lines for points, get polynomial coefficients

N = length(pts);
perp_slopes  = -1/poly_coeffs(1)*ones(N,1); %sample slope for each point
y_intercepts = pts(:,2) - perp_slopes.*pts(:,1); %different y-intercepts

%% Get new x and y coordinates

%equate fit line and perpendicular line, solve for x-coordinates
x_fit = (y_intercepts - poly_coeffs(2))./(poly_coeffs(1) - perp_slopes);

%determine y-coordinates by passing x-coordinates to equation for fit
y_fit = polyval(poly_coeffs,x_fit);

%% Compute variance

%{
This variance is computed assuming that the line that fits through the two
vector quantities represents the "best-fit" solution.
%}

s = sqrt(N/(N^2 - N)*sum((pts(:,1)-x_fit).^2 + (pts(:,2)-y_fit).^2));
%% Determine percent distances of new x, y coords with vertices

%sum in quadriture (i.e. pythagorean distance)
whole_dist = sqrt((v2(2) - v1(2))^2 + (v2(1) - v1(1))^2);
v1_dists   = sqrt((v1(2) - y_fit).^2 + (v1(1) - x_fit).^2);
v2_dists   = sqrt((v2(2) - y_fit).^2 + (v2(1) - x_fit).^2);

%values that exceed the whole_dist will be equal to whole_dist
v1_exc = v1_dists > whole_dist;
v2_exc = v2_dists > whole_dist;
v1_dists(v1_exc) = whole_dist;
v2_dists(v1_exc) = 0;
v2_dists(v2_exc) = whole_dist;
v1_dists(v2_exc) = 0;

%percentages
perc_v1 = v2_dists/whole_dist;
perc_v2 = v1_dists/whole_dist;

end