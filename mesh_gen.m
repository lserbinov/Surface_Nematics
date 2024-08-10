close all; clear all;
[P,M,~,~,~,~] = readOBJ("./assets/genus6.obj");
% [P,M,~,~] = loop(P,M,1);
% [P,M] = subdivided_sphere(3);
% P(:, 1) = 2 * P(:, 1);
figure;
geo = Geometry(M, P);
IO.show(M, P, dot(geo.bending_force(1), geo.v_normal, 2));
hold on;
qvr(geo.f_center, geo.f_normal);

%%

P = IO.smoothing(M, P, 0.005, 1e-2); % sphere
% save("./assets/bob.mat", "M", "P");

% load("./assets/bob.mat");
% M = [M(:, 3), M(:, 2), M(:, 1)];
% save("./assets/bob.mat", "M", "P");
geo = Geometry(M, P);
figure;
IO.show(M, P, dot(geo.bending_force(1), geo.v_normal, 2));
hold on;
%qvr(geo.f_center, geo.f_normal);
save("./assets/genus6_smooth.mat", "M", "P");