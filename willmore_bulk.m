verbose = false;
%%% directory
dir = "./data/bulk/sphere/strainrot_alpha5e-1_kappa1e-2/"; 
[status, msg, msgID] = mkdir(dir); 

%% system init
%%% 0 for new simulation, otherwise continue from geo"start".mat
start = 0; 
if start == 0
    %%% geometry
    % load("./assets/gut_smooth.mat")
    % P = [P(:, 3), P(:, 2), -P(:, 1)];
    [P, M] = subdivided_sphere(4);
    geo = Geometry(M, P);
    [P, ~] = rm_rigid(P, zeros(size(P)), geo.v_area);
    geo = Geometry(M, P);
    %%% parameters
    p.dt = 1e-2; % time 
    p.kappa = 1e-2; % bending
    p.T = 3000; % total time (frames)
    p.lambda = 1e3; % area penalty
    p.alpha = 5e-1; p.amp = 1; p.type = "strainrot"; p.scale = 5;
    %%% optimizer
    o.h = 0.1; o.eta = 1e2; o.k = 0; o.metric = "l2"; 
    o.tol_f = 5e-4; o.tol_d = 5e-4; o.max_iter = 10000;
    %%% remesh
    r.edge_length = mean(geo.he_length); % target edge length
    r.n_iter = 50; % remeshing iterations
    r.smooth = 0.01; % smoothing factor
    %%% initialize
    velocity = zeros(size(P, 1) * 3, 1);
    pressure = zeros(size(M, 1), 1);
else
    load(dir + sprintf("geo%d.mat", start), ...
     "M", "P", "velocity", "pressure", "p", "o", "r"); 
    geo = Geometry(M, P);
end

%% main loop (not supposed to be modified)
for t = (start + 1):p.T
    %%% incremental potential minimization
    [~, K, ~, div, KTK, DTD] = geo.evolving_operators();           
    mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
    mass0 = blkdiag(mass0, mass0, mass0);
    P0 = P(:); P(:) = P0 + p.dt * velocity;
    u0 = bulk_velocity(P0, p.amp, p.type, p.scale);
    eps_f = Inf; eps_d = Inf; j = 0;
    while ((eps_f > o.tol_f) || (eps_d > o.tol_d)) && (j < o.max_iter)
        %%% update bending force
        fb = geo.bending_force(p.kappa);
        b = - (2 * KTK + p.alpha * mass0) * (P(:) - P0) + p.dt * ( fb(:) + p.alpha * mass0 * u0 + div' * pressure);
        expan = div * (P(:) - P0);
        %%% measure residual
        eps_f = norm_f(b, geo.v_area); eps_d = norm_d(expan, geo.f_area);
        if verbose fprintf("t = %d, j = %d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d); end
        %%% gradient descent/ascent
        P(:) = P(:) + o.h * ((preconditioner(geo, o.metric) + 1e-5 * sparse(eye(3*geo.mesh.n_v))) ...
                            \ (b - o.k * DTD * (P(:) - P0)));
        pressure = pressure - o.eta * o.h * expan ./ geo.f_area;
        %%% update geometry
        geo = Geometry(M, P);
        j = j + 1;
    end
    %%% save data
    velocity = (P(:) - P0) / p.dt;
    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "pressure", "fb", "u0", "p", "o", "r");
    fprintf("Save geo%d.mat at j =%d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d);
    %%% update geometry, remesh if needed
    geo = Geometry(M, P);
    if ~geo.is_delaunay(0)
        geo_pre = geo; 
        fprintf("Remeshing. t = %d \n", t);
        [P, M] = IO.remesh(P, M, r, dir, t);
        geo = Geometry(M, P);
        [velocity, pressure] = map_data(geo, geo_pre, velocity, pressure);
    end
end

%% helper functions
function H = preconditioner(geo, metric)
    % preconditioner for incremental potential minimization
    switch  metric 
        case "bih"
            mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
            mass0_inv = spdiags(1./geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
            H = mass0 * geo.lap * mass0_inv * geo.lap;
            H = blkdiag(H, H, H);
        case "lap"
            H = geo.lap;
            H = blkdiag(H, H, H);
        case "stks"
            H = KTK;
        case "l2"
            H = speye(3*geo.mesh.n_v);
    end
end

function u = bulk_velocity(coord, amp, type, scale)
    coord = reshape(coord, [], 3);
    nv = size(coord, 1);
    u = zeros(3 * nv, 1);
    switch type
        case "shear"
            u = [zeros(nv, 1); coord(:, 3) * amp; zeros(nv, 1)];
        case "abc"
            A = amp; B = amp; C = amp;
            X = coord(:, 1); Y = coord(:, 2); Z = coord(:, 3);
            u = [A * sin(scale * Z) + C * cos(scale * Y); ...
                 B * sin(scale * X) + A * cos(scale * Z); ...
                 C * sin(scale * Y) + B * cos(scale * X)];
        case "strainrot"
            P = carte2cylin(coord);
            A = amp;
            u_rot = 5 * [zeros(nv, 1), A * P(:, 1) .* P(:, 3), zeros(nv, 1)];
            u_strain = [- 0.5 * P(:, 1), zeros(nv, 1), A * P(:, 3)];
            u = u_rot + u_strain;
            [u, P] = cylin2carte(u, P);
            u = u(:);
    end
end

function e = norm_f(b, area)
    % normalized norm of vector field
    e = sqrt(sum(b.^2 ./ [area; area; area])) / sum(area);
end

function e = norm_d(expan, area)
    % normalized norm of scalar field
    e = sqrt(sum(expan.^2 ./ area)) / sum(area);
end

function [velocity, pressure] = map_data(geo, geo_pre, velocity_pre, pressure_pre)
    % interpolate data from previous geometry to current geometry
    kdtree = KDTreeSearcher(geo_pre.f_center);
    %%% interpolate vertex data - velocity
    [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.V, kdtree, 6);
    if fail
        error("projection failed.");
    end
    velocity = interpolate(geo_pre.F, face, uv, reshape(velocity_pre, [], 3));
    velocity = velocity(:);
    %%% interpolate face data - pressure
    [face, uv, count, fail] = project(geo_pre.V, geo_pre.F, geo.f_center, kdtree, 6);
    if fail
        error("projection failed.");
    end
    [pressure_pre_v, neighbor] = geo_pre.mesh.face_to_vertex(pressure_pre);
    pressure_pre_v = pressure_pre_v ./ neighbor;
    pressure = interpolate(geo_pre.F, face, uv, pressure_pre_v);
end

function [velocity, pressure] = aug_lag_guess(face, vertex, parameter)
    % instantaneous Stokes solution using augmented lagrangian
    geo = Geometry(face, vertex);
    id = sparse(eye(3*size(vertex, 1)));
    force = geo.bending_force(parameter.kappa);
    pressure = zeros(size(face, 1), 1);
    [~, ~, ~, div, KTK, DTD] = geo.evolving_operators();
    for i=1:30
        b = (force(:) + div' * pressure);
        velocity = (KTK +  1e-5 * id + 1 * DTD) \ b;
        % [velocity, fail] = cgs(KTK + rho * DTD +  1e-2 * id , b, 1e-4, 10000, [], [], velocity);
        pressure = pressure - 1 * (div * velocity) ./ geo.f_area;
    end
    [~, velocity] = rm_rigid(vertex, velocity, geo.v_area);
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
