function F = pushforward(mesh, vertex1, vertex2)
    % face-based pushforward (deformation gradient)
    % 
    % Inputs:
    %   mesh: #F by 3 face matrix
    %   vertex1: #V by 3 vertex matrix 1
    %   vertex2: #V by 3 vertex matrix 2
    % Outputs:
    %   F: #F by 2 by 2, pushforward from vertex1 to vertex2 in local frame
    
    n_f = size(mesh, 1);
    function v = get_basis(vertex)
        v = zeros(2, 2, n_f);
        e_a = vertex(mesh(:, 2), :) - vertex(mesh(:, 1), :);
        e_b = vertex(mesh(:, 3), :) - vertex(mesh(:, 1), :);
        l_a = vecnorm(e_a, 2, 2);
        l_b = vecnorm(e_b, 2, 2);
        t = acos(dot(e_a ./ l_a, e_b ./ l_b, 2));
        v(1, 1, :) = l_a;
        v(:, 2, :) = [cos(t) .* l_a, sin(t) .* l_b]';
    end
    v1 = get_basis(vertex1);
    v2 = get_basis(vertex2);
    F = permute(pagemrdivide(v2, v1), [3, 1, 2]);
end