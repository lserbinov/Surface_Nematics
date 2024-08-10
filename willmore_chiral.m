verbose = true;
%%% directory
dir = "./data/chiral/sphere/"; 
[status, msg, msgID] = mkdir(dir); 

%% system init
%%% 0 for new simulation, otherwise continue from geo"start".mat
start = 0; 
if start == 0
    %%% geometry
    % [P,M,~,~,~,~] = readOBJ("./assets/sphere.obj");
    % P(:, 1) = P(:, 1) * 2;
    load("./assets/gut_smooth.mat");
    geo = Geometry(M, P);
    %%% parameters
    p.dt = 1e-3; % time 
    p.kappa = 1e-1; % bending
    p.T = 4000; % total time (frames)
    p.alpha = 50; % chiral friction
    p.w = 10; % angular velocity
    %%% optimizer
    o.h = 0.05; o.eta = 1e2; o.k = 0; o.metric = "lap"; 
    o.tol_f = 1e-3; o.tol_d = 1e-3; o.max_iter = 10000;
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
    [pgrad, K, W, div, KTK, DTD] = geo.evolving_operators();
    P0 = P(:); P(:) = P0 + p.dt * velocity;
    eps_f = Inf; eps_d = Inf; j = 0;
    omega = spin_tensor(geo.f_normal, p.w);
    mass2_inv = sptensor([(1:geo.mesh.n_f)', (1:geo.mesh.n_f)'], ...
                            1./geo.f_area, [geo.mesh.n_f, geo.mesh.n_f]);
    while ((eps_f > o.tol_f) || (eps_d > o.tol_d)) && (j < o.max_iter)
        %%% update chiral force 
        WP = ttt(W, sptensor(P - reshape(P0, [], 3)), [2, 4], [1, 2]);
        WP = ttt(mass2_inv, sptensor(WP), 2, 1);
        disp("WP")
        WP(1, :, :)
        disp("omega")
        omega(1, :, :)
        fa = p.alpha * double(ttt(pgrad, sptensor(WP - sptensor(omega) * p.dt), [1, 3, 5], [1, 2, 3]));
        %%% update bending force
        fb = geo.bending_force(p.kappa);
        b = - 2 * KTK * (P(:) - P0) - p.dt * (- fb(:) - div' * pressure) - fa(:);
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
    [P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "pressure", "fb", "fa", "p", "o", "r");
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
            H = sparse(eye(3*geo.mesh.n_v));
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

function omega = spin_tensor(f_normal, angular_velocity)
    vorticty = 2 * angular_velocity * f_normal;
    omega = zeros(size(f_normal, 1) ,3, 3);
    omega(:, 1, 2) = -vorticty(:, 3);
    omega(:, 1, 3) = vorticty(:, 2);
    omega(:, 2, 1) = vorticty(:, 3);
    omega(:, 2, 3) = -vorticty(:, 1);
    omega(:, 3, 1) = -vorticty(:, 2);
    omega(:, 3, 2) = vorticty(:, 1);
    omega = omega / 2;
end