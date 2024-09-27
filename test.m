% Area = geo.area;
% volume = geo.volume;
close all; clc; clear all;
[V, F] = subdivided_sphere(3);
IO.show(F,V,ones(size(V, 1), 1));
hold on;
geo = Geometry(V,F);
% load("data/willmore/geo1.mat")

% IO.show(M, P, rand(size(P, 1), 1));
% qvr(P, rand(size(P)),'off');
q_ = (initialize_nematic(V,F, 10));
p_ = (initialize_nematic(V,F, 5));
% inner product 
inner = real(conj(q_) .* p_);
% L2 inner product 
L2_inner = sum(real(conj(q_) .* (geo.mass2 * p_)));
a = sqrt(q_);
q = [real(a),imag(a)];
Q = geo.f_basis_u .* q(:, 1) + geo.f_basis_v .* q(:, 2);
qvr(geo.f_center,Q)
% figure 
load("data/willmore/geo100.mat")

% IO.show(M, P, ones(size(P, 1), 1));

%%
function q_bra_hat = initialize_nematic(face, vertex, rank)
    % initialize nematic field by rayleigh quotient
    k = 2; %% k-atic
    geo = Geometry(face, vertex);
    L = geo.bochner_laplacian(k);
    [V, ~] = eigs(L, geo.mass2, rank, "smallestabs");
    q_bra_hat = V(:, rank);
    % q_bra_hat = V(:, rank) ./ vecnorm(V(:, rank), 2, 2);
end

