close all;
clc;
clear all;

% load("./assets/genus6_smooth.mat")
global P v_area
[P, M] = subdivided_sphere(4);
geo = Geometry(M, P);
v_area = geo.v_area;
[~, K, ~, div, KTK, DTD] = geo.evolving_operators();

rank = 10;
scale = 10;
mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
mass0 = blkdiag(mass0, mass0, mass0);
[V, D] = eigs(KTK, mass0, rank, "smallestabs");
figure;
b = bar(abs(diag(D)));
b.FaceColor = '#ACB8B0';
set(gca,'YScale','log')
Xlabel = xlabel("$i$", "Interpreter", "latex");
Ylabel = ylabel("$|\lambda_{i}|$", "Interpreter", "latex");
set(gca, 'FontName', 'Helvetica', 'FontSize', 12);
set([Xlabel, Ylabel], 'FontName', 'Helvetica', 'FontSize', 20);

figure; 
IO.show(M, P, ones(size(P, 1), 1));
hold on; 
qvr(P, reshape(real(V(:, 1)), [], 3), scale);
qvr(P, reshape(real(V(:, 2)), [], 3), scale);
qvr(P, reshape(real(V(:, 3)), [], 3), scale);
qvr(P, reshape(real(V(:, 4)), [], 3), scale);
qvr(P, reshape(real(V(:, 5)), [], 3), scale);
qvr(P, reshape(real(V(:, 6)), [], 3), scale);

figure;
IO.show(M, P, ones(size(P, 1), 1));
hold on; 

V1_rm = local_rm_rigid(reshape(real(V(:, 1)), [], 3));

qvr(P, local_rm_rigid(reshape(real(V(:, 1)), [], 3)), scale);
qvr(P, local_rm_rigid(reshape(real(V(:, 2)), [], 3)), scale);
qvr(P, local_rm_rigid(reshape(real(V(:, 3)), [], 3)), scale);
qvr(P, local_rm_rigid(reshape(real(V(:, 4)), [], 3)), scale);
qvr(P, local_rm_rigid(reshape(real(V(:, 5)), [], 3)), scale);
qvr(P, local_rm_rigid(reshape(real(V(:, 6)), [], 3)), scale);

figure;
IO.show(M, P, ones(size(P, 1), 1));
hold on; 
qvr(P, rm_translation(reshape(real(V(:, 1)), [], 3)), 3 * scale);
% qvr(P, rm_translation(reshape(real(V(:, 2)), [], 3)), 2 * scale);
% qvr(P, rm_translation(reshape(real(V(:, 3)), [], 3)), 2 * scale);
qvr(P, rm_rotation(reshape(real(V(:, 1)), [], 3)), 2 * scale);
% qvr(P, rm_rotation(reshape(real(V(:, 2)), [], 3)), 2 * scale);
% qvr(P, rm_rotation(reshape(real(V(:, 3)), [], 3)), 2 * scale);
qvr(P, rm_translation(reshape(real(V(:, 4)), [], 3)), 3 * scale);
qvr(P, rm_translation(reshape(real(V(:, 5)), [], 3)), 3 * scale);
% qvr(P, rm_translation(reshape(real(V(:, 6)), [], 3)), 2 * scale);
qvr(P, rm_rotation(reshape(real(V(:, 4)), [], 3)), 2 * scale);
% qvr(P, rm_rotation(reshape(real(V(:, 5)), [], 3)), 2 * scale);
qvr(P, rm_rotation(reshape(real(V(:, 6)), [], 3)), 2 * scale);

function velocity = local_rm_rigid(vel)
    global P v_area
    [~, velocity] = rm_rigid(P, vel, v_area);
    velocity = reshape(velocity, [], 3);
end


function velocity = rm_translation(velocity)
    global P v_area
    vertex = P;

    % format velocity 
    n_v = size(vertex, 1);
    velocity = reshape(velocity, n_v, 3);

    % remove translation
    mass = sum(v_area); % total mass
    u_L = sum(velocity .* v_area) / mass; % trans vel

    velocity = velocity - repmat(u_L, n_v, 1);
end


function velocity = rm_rotation(velocity)
    global P v_area
    vertex = P;

    % format velocity 
    n_v = size(vertex, 1);
    velocity = reshape(velocity, n_v, 3);

    % remove translation
    mass = sum(v_area); % total mass
    com = sum(vertex .* v_area) / mass; % center of mass

    % shift com and remove translation
    vertex = vertex - repmat(com, n_v, 1);
    ang_momentum = sum(cross(vertex, velocity, 2).* v_area, 1);
    
    % construct moment of inertia tensor
    Rsq_id = eye(3);
    Rsq_id = permute(Rsq_id(:, :, ones(n_v, 1)), [3, 1, 2]);
    Rsq_id = Rsq_id .* sum(vertex .* vertex, 2);
    VVT = vertex(:, :, ones(3, 1)) .* permute(vertex(:, :, ones(3, 1)), [1, 3, 2]);
    moment = sum((Rsq_id - VVT) .* v_area, 1);

    % remove rotation
    w = (squeeze(moment) \ ang_momentum')'; % angular velocity
    w = repmat(w, n_v, 1);
    velocity = velocity - cross(w, vertex, 2);
end









