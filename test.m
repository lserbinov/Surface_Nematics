% [V, F] = subdivided_sphere(3);
% IO.show(F,V,ones(size(V, 1), 1));
% geo = Geometry(F, V);
% Area = geo.area;
% volume = geo.volume;
close all; clc;
load("data/willmore/geo1.mat")

IO.show(M, P, ones(size(P, 1), 1));

figure 
load("data/willmore/geo100.mat")

IO.show(M, P, ones(size(P, 1), 1));
