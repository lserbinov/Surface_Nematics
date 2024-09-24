verbose = false;
%%% directory
dir = "./data/gut1/"; 
[status, msg, msgID] = mkdir(dir); 

%% system init
%%% 0 for new simulation, otherwise continue from geo"start".mat
start = 0; 
if start == 0
    %%% geometry
    % load("./assets/gut_smooth.mat")
    geo = Geometry(M, P);
    % p.rank = 1; % nematic init based on eigen rank
    %%% parameters
    p.dt = 1e-2; % time 
    % p.kappa = 1e-2; % bending
    % p.alpha = -10; % activity: extensile (+), contractile (-)
    p.T = 2000; % total time (frames)
    %%% optimizer
    o.h = 1e-1; o.eta = 1e2; o.k = 0; o.metric = "lap"; 
    o.tol_f = 5e-4; o.tol_d = 5e-4; o.max_iter = 10000;
    %%% remesh
    r.edge_length = mean(geo.he_length); % target edge length
    r.n_iter = 50; % remeshing iterations
    r.smooth = 0.01; % smoothing factor
    %%% initialize
    velocity = zeros(size(P, 1) * 3, 1);
    pressure = zeros(size(M, 1), 1);
    nematic = initialize_nematic(M, P, 2);
else
    load(dir + sprintf("geo%d.mat", start), ...
     "M", "P", "velocity", "pressure", "nematic_", "p", "o", "r"); 
    geo = Geometry(M, P);
    nematic = v2c(V2v(geo, nematic_)).^2;
end

%% main loop (not supposed to be modified)
for t = (start + 1):p.T
    %%% incremental potential minimization
    [~, K, ~, div, KTK, DTD] = geo.evolving_operators();
    P0 = P(:); P(:) = P0 + p.dt * velocity;
    fa = sign(p.alpha) * double(ttt(K, sptensor(veronese(geo, nematic)), [1, 3, 5], [1, 2, 3]));
    eps_f = Inf; eps_d = Inf; j = 0;
    while ((eps_f > o.tol_f) || (eps_d > o.tol_d)) && (j < o.max_iter)
        %%% update bending force
        fb = geo.bending_force(p.kappa);
        b = - 2 * KTK * (P(:) - P0) - p.dt * (- fa(:) - fb(:) - div' * pressure);
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
    %%% diffusion
    nematic = (geo.mass2 + p.dt / (1e-8 + abs(p.alpha)) * geo.bochner_laplacian(2)) \ (geo.mass2 * nematic); 
    %%% advection
    nematic = advect(M, reshape(P0, [], 3), reshape(P, [], 3), nematic, 2);
    %%% normalize
    nematic = nematic ./ abs(nematic);
    %%% save data
    [P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    nematic_ = v2V(geo, c2v(nematic.^(1/2)));
    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "pressure", "nematic_", "fa", "fb", "p", "o", "r");
    fprintf("Save geo%d.mat at j =%d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d);
    %%% update geometry, remesh if needed
    geo = Geometry(M, P);
    if ~geo.is_delaunay(0)
        geo_pre = geo; 
        fprintf("Remeshing. t = %d \n", t);
        [P, M] = IO.remesh(P, M, r, dir, t);
        geo = Geometry(M, P);
        [velocity, pressure, nematic] = map_data(geo, geo_pre, velocity, pressure, nematic);
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

function q_bra_hat = initialize_nematic(face, vertex, rank)
    % initialize nematic field by rayleigh quotient
    k = 2; %% k-atic
    geo = Geometry(face, vertex);
    L = geo.bochner_laplacian(k);
    [V, ~] = eigs(L, geo.mass2, rank, "smallestabs");
    q_bra_hat = V(:, rank) ./ vecnorm(V(:, rank), 2, 2);
end

function [velocity, pressure, q_bra_hat] = map_data(geo, geo_pre, velocity_pre, pressure_pre, q_bra_hat_pre)
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
    %%% interpolate face data - nematic
    QQ_pre = veronese(geo_pre, q_bra_hat_pre);
    [QQ_pre_v, neighbor] = geo_pre.mesh.face_to_vertex(reshape(QQ_pre, geo_pre.mesh.n_f, []));
    QQ_pre_v = QQ_pre_v ./ neighbor;
    QQ = interpolate(geo_pre.F, face, uv, QQ_pre_v);
    q_bra_hat = veronese_inv(geo, QQ);
end

function QQ = veronese(geo, q_bra_hat)
    % veronese map from complex to tensor
    Q = v2V(geo, c2v(q_bra_hat.^(1/2)));
    QQ = Q(:, :, ones(1,3)) .* permute(Q(:,:,ones(1,3)), [1,3,2]);
end

function q_bra_hat = veronese_inv(geo, QQ)
    % inverse veronese map from tensor to complex
    eigen = @(i) eigs(reshape(QQ(i, :), 3, 3), 1 ,'largestabs');
    [Q, ~] = arrayfun(eigen, 1:geo.mesh.n_f, 'UniformOutput', false);
    Q = cell2mat(Q)';
    q = V2v(geo, Q);
    q_bra_hat = v2c(q ./ vecnorm(q, 2, 2)).^2;
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

function Fq_bra_hat = advect(mesh, vertex1, vertex2, q_bra_hat, k)
    % Lie advection of nematic field
    q_hat = q_bra_hat.^(1/k);
    F = pushforward(mesh, vertex1, vertex2);
    q = c2v(q_hat);
    Fq = squeeze(pagemtimes(permute(F, [2, 3, 1]), ...
                           permute(q, [2, 3, 1])))';
    Fq = Fq ./ vecnorm(Fq, 2, 2); % normalize
    Fq_bra_hat = v2c(Fq).^k;
end

function q = V2v(geo, Q)
    % R3 realization to local chart
    q = [dot(geo.f_basis_u, Q, 2), dot(geo.f_basis_v, Q, 2)];
end

function Q = v2V(geo, q)
    % R3 realization of tangent vector
    Q = geo.f_basis_u .* q(:, 1) + geo.f_basis_v .* q(:, 2);
end

function q = c2v(q_hat)
    % complex to vector
    q = [real(q_hat), imag(q_hat)];
end
function q_hat = v2c(q)
    % vector to complex
    q_hat = q(:,1) + 1i * q(:,2);
end

