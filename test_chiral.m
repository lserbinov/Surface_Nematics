close all; 
clc;
load("./assets/gut_smooth.mat")
geo = Geometry(M, P);
p.alpha = 1;
p.w = 1;
[pgrad, K, W, div, KTK, DTD] = geo.evolving_operators();
mass2_inv = sptensor([(1:geo.mesh.n_f)', (1:geo.mesh.n_f)'], 1./geo.f_area, [geo.mesh.n_f, geo.mesh.n_f]);

%   pgrad: 5d sparse tensor, (F, V, T*M, R3*, T*M) extrinsic gradient 
%   K: 5d sparse tensor, (F, V, [T*M, R3*, T*M]) Killing
%   operator (symmetric part of pgrad)
%   W: 5d sparse tensor, (F, V, \T*M, R3*, T*M\) antisymmetric part of pgrad
%   div: F by 3*V sparse matrix, (F, V, R3*) divergence operator
%   KTK: 3*V by 3*V sparse matrix, (V, R3, V, R3) viscosity Laplacian
%   DTD: 3*V by 3*V sparse matrix, gradient of divergence operator
P0 = P(:); P(:) = P0 + rand(size(P0));

WU = ttt(W, sptensor(P - reshape(P0, [], 3)), [2, 4], [1, 2]);
WU = ttt(mass2_inv, sptensor(WU), 2, 1);
omega = spin_tensor(geo.f_normal, p.w);
fa1 = p.alpha * double(ttt(K, sptensor(omega - WU), [1, 3, 5], [1, 2, 3]));
fa2 = p.alpha * double(ttt(pgrad, sptensor(omega - WU), [1, 3, 5], [1, 2, 3]));
fa3 = p.alpha * double(ttt(W, sptensor(omega - WU), [1, 3, 5], [1, 2, 3]));

KU = ttt(K, sptensor(P - reshape(P0, [], 3)), [2, 4], [1, 2]);
KU = ttt(mass2_inv, sptensor(KU), 2, 1);
fv1 = double(ttt(K, sptensor(KU), [1, 3, 5], [1, 2, 3]));
fv_ = reshape(KTK * (P(:) - P0), [], 3);
fv2 = double(ttt(pgrad, sptensor(KU), [1, 3, 5], [1, 2, 3]));
fv3 = double(ttt(W, sptensor(KU), [1, 3, 5], [1, 2, 3]));


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