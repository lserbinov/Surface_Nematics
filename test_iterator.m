close all;
clc;
clear all;

file = "/Users/cuncheng/Dev/active_nematics_shell/matlab/data/nematic/bob/diffusion8/geo5.mat";
load(file);
geo = Geometry(M, P);

cumsum([1, 2, 3]);

% Define the vector
vec = [1, 6, 9, 10, 1, 6];

% Remove repeating indices while preserving the order
unique_vec = unique(vec, 'stable');

% the assignment will follow the order!!
a(1, [1,1]) = [2, 3]; 


% 1) find the max number N of neighbors per vertex in the mesh by
%        N = accumarray(he_src,1)
% 2) find the first he_src ocurrance of he for each src point
%        v_he1 = accumarray(he_src, (1:n_he)', [n_v, 1], @min);
% 3) find the list of oriented neighborhood list by doing
%                                                        he_next(he_flip())
%       for N times, call "iterator"
% 4) do cumsum, say on corner angle A, along the second dimension (with
% length N) AA
% 5) flip along the second dimension (with N),get this n_v * N polar angle
% matrix AA_FLIP
% 6) do he_polar_angle(flip(iterator(:))) = AA_FLIP(:) so that the turning
% aroud version will be overwritten by the first-timer value. 



% 4) given attribute, say corner angle A, have function 
%       polar_angle = @(i) cumsum(A(unique(iterator(1, :),'stable'));
% 5) run for each vertex: arrayfun(polar_angle, 1:n_v)

%% polar angle 
    % int he0 = i@halfedge;
    % int he = he0;
    % int he_flip;
    % float phi = 0;
    % do{
    %   setvertexattrib(geoself(),"point_polar_angle",he,-1,phi); // argument list is (geohandle, string name, int prim_num, int vertex_num, <type>value)，set vertex_num = 1 and use prim_num as a linear index
    %   int he_flip = vertex(0, "flip", he);
    %   if (he_flip == -1) break;
    %   he = vertex(0,"next",he_flip); 
    %   phi += f@polar_turning_factor*vertex(0,"angle",he);
    % }while(he!=he0);

%% connection 
% int point = vertex(0, "point", @vtxnum);
% int next_point = vertex(0, "next_point", @vtxnum);
% int n = `chs("../../n_symmetric")`;
% if (point(0, "is_boundary", point) || point(0, "is_boundary", next_point)){
%     // f@n_symmetric_connection = 0;
%     setvertexattrib(0, sprintf("symmetric_%d_connection",n), @vtxnum, -1, 0.0, "set");
% } else {
%     float phi_flip = vertex(0,"point_polar_angle",i@flip);
%     // f@n_symmetric_connection = phi_flip - f@point_polar_angle + PI;
%     // f@n_symmetric_connection*= n; // for four-fold symmetric cross field 
%     float connection = phi_flip - f@point_polar_angle + PI;
%     connection *= n; // for four-fold symmetric cross field 
%     setvertexattrib(0, sprintf("symmetric_%d_connection",n), @vtxnum, -1, connection, "set");
% }

%% Bochner Laplacian
% import numpy as np
% import scipy.sparse as sp
% import scipy.sparse.linalg as la
% 
% node = hou.pwd()
% geo = node.geometry()
% 
% cache_node = hou.node(hou.ch("../cache_path")) # node to store the resulting matrices
% 
% root_path = hou.ch("../root_path")
% n_symm = hou.ch("../n_symmetric")
% 
% dec_node = hou.node(root_path + "../cuncheng_build_dec_matrices/")
% star1 = dec_node.cachedUserData("star1")
% 
% # Read geo info
% npoints = hou.Geometry.intrinsicValue(geo,"pointcount")
% nhedges = hou.Geometry.intrinsicValue(geo,"vertexcount")
% nprims = hou.Geometry.intrinsicValue(geo, "primitivecount")
% 
% # Read attributes
% point = np.array(geo.vertexIntAttribValues("point"))
% next_point = np.array(geo.vertexIntAttribValues("next_point"))
% vertex_connection = np.array(geo.vertexFloatAttribValues("symmetric_" + str(n_symm) + "_connection")) 
% 
% # Useful arrays
% point_ind = np.arange(npoints) # [0,1,...,npoints-1]
% hedge_ind = np.arange(nhedges) # [0,1,...,nhedges-1]
% prim_ind = np.arange(nprims)
% 
% # Build d0
% d0row = np.append(hedge_ind,hedge_ind)
% d0col = np.append(point,next_point)
% d0val_pt_conn = np.append(-np.ones(nhedges),np.exp(-1j*vertex_connection))
% d0_pt_conn = sp.csr_matrix((d0val_pt_conn,(d0row,d0col)),(nhedges,npoints))
% 
% # Build laplacian0
% L0_pt_conn = d0_pt_conn.getH().dot(star1.dot(d0_pt_conn)) # getH is the Hermitian transpose
% 
% # cache matrices
% cache_node.setCachedUserData("L0_connection",L0_pt_conn)
% cache_node.setCachedUserData("d0_connection",d0_pt_conn)



%% defect detection 
% int he0 = i@halfedge;
% int he = he0;
% float R = 0;
% float turning = 0;
% // int counter = 0;
% // float angle_sum = 0;
% int n_symm = `chs("../n_symm")`;
% if (i@is_boundary){
%     f@poincare = 0;
% } else {
%     do {
%         // counter++;
%         float phase1 = prim(0, "phase", vertex(0, "primitive", he));
%         float connection = vertex(0,sprintf("symmetric_%d_perpendicular_connection", n_symm),he);
%         he = vertex(0, "next", vertex(0, "flip", he));
%         float phase2 = prim(0, "phase", vertex(0, "primitive", he));
%         float eta = (phase2 - phase1 - connection);
%         eta = atan2(sin(eta), cos(eta));
%         R += connection;
%         turning += eta; 
%     } while (he != he0);
%     R = atan2(sin(R),cos(R));
%     f@poincare = (R + turning)/(2*PI);
%     // f@counter = counter;
% }

function q_bra_hat = initialize_nematic(face, vertex, rank)
    % initialize nematic field by rayleigh quotient
    k = 2; %% k-atic
    geo = Geometry(face, vertex);
    L = geo.bochner_laplacian(k);
    [V, ~] = eigs(L, geo.mass2, rank, "smallestabs");
    q_bra_hat = V(:, rank) ./ vecnorm(V(:, rank), 2, 2);
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

