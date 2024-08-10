clear all;
clc;
close all;
[P,M,~, ~, ~, ~] = readOBJ("~/Dev/active_nematics_shell/assets/models/bunny_regular.obj");
P = [P(:, 1), P(:, 3), P(:, 2)];
% [P,M] = subdivided_sphere(3);
% P(:, 1) = P(:, 1) * 2;
geo = Geometry(M, P);
figure;
data = ones(geo.mesh.n_v, 1);
IO.show(M, P, data);
p = 3;
v0 = geo.volume;
dt = 0.001;
A = sparse(1:geo.mesh.n_v, 1:geo.mesh.n_v, geo.v_area , geo.mesh.n_v, geo.mesh.n_v);
for ii = 1:10
    geo = Geometry(M, P); 
    P = (A + dt * geo.lap) \ (A * (P - dt * p * (geo.volume - v0) * geo.v_normal));
end
figure
IO.show(M, P, data);


