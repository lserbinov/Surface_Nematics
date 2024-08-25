close all; 
clc;
clear all;
load("./assets/gut_smooth.mat")
geo = Geometry(M, P);
[P, ~] = rm_rigid(P, zeros(size(P)), geo.v_area);

x = linspace(-1, 1, 10);
y = linspace(-1, 1, 10);
z = linspace(-1, 1, 10);
[x, y, z] = meshgrid(x, y, z);
P = [x(:), y(:), z(:)];

% changed variables

u = strain(P) + rot_shear(P);
qvr(P, reshape(u, [], 3), 3);





function u = strain(P)
    P = carte2cylin(P);
    A = 1;
    nv = size(P, 1);
    u = [-A * P(:, 1), zeros(nv, 1), 2 * A * P(:, 3)];
    [u, ~] = cylin2carte(u, P);
end

function u = ABC(P)
    % take cartesian coordinate
    A = 1;
    B = 1;
    C = 1;
    scale = [1, 1, 1];
    u = [A * sin(scale(3) * P(:, 3)) + C * cos(scale(2) * P(:, 2)), ...
        B * sin(scale(1) * P(:, 1)) + A * cos(scale(3) * P(:, 3)), ...
        C * sin(scale(2) * P(:, 2)) + B * cos(scale(1) * P(:, 1))];
end

function u = shear(coord)
    % take cartesian coordinate
    coord = reshape(coord, [], 3);
    nv = size(coord, 1);
    u = [zeros(nv, 1); coord(:, 3); zeros(nv, 1)];
end

function u = rot_shear(P)
    P = carte2cylin(P);
    A = 1;
    nv = size(P, 1);
    u = [zeros(nv, 1), A * P(:, 1) .* P(:, 3), zeros(nv, 1)];
    [u, ~] = cylin2carte(u, P);
end

function [u_, P_] = cylin2carte(u, P)
    % [u_r, u_t, u_z], [P_r, P_t, P_z] -> [u_x, u_y, u_z], [P_x, P_y, P_z]
    u_ = [u(:, 1) .* cos(P(:, 2)) - u(:, 2) .* sin(P(:, 2)), ...
        u(:, 1) .* sin(P(:, 2)) + u(:, 2) .* cos(P(:, 2)), ...
        u(:, 3)];
    P_ = [P(:, 1) .* cos(P(:, 2)), ...
        P(:, 1) .* sin(P(:, 2)), ...
        P(:, 3)];
end

function [P_] = carte2cylin(P)
    P_ = [sqrt(P(:, 1).^2 + P(:, 2).^2), atan2(P(:, 2), P(:, 1)), P(:, 3) ];
end


