classdef IO
    methods(Static)
        function res = show(mesh, vertex, data)
            t = tsurf(mesh, vertex, ...
                'cData', data, ...
                'FaceColor','interp', ...
                'FaceLighting','phong', ...
               'EdgeColor','none',...
                    'FaceAlpha', 1, ...
                fphong,'SpecularStrength',0.1);
             
                % 'EdgeColor',[0.2 0.2 0.2]);
            % colormap(flag(8))
            % t = tsurf(mesh,vertex,'EdgeColor','none',fphong,'SpecularStrength',0.1);
            l = [
                light('Position',[-10 -10 13],'Style','local');
                light('Position',[10 -10  13],'Style','local')];
            camproj('persp');
            axis equal;
            axis off;
            % res = add_shadow(t,l);
            hold on;
        end
        function h = scatter(vertex)
            h = scatter3(vertex(:, 1), vertex(:, 2), vertex(:, 3), 10, 'red', 'filled');
        end
        function vertex = smoothing(mesh, vertex, strength, duration)
            geo = Geometry(mesh, vertex);
            kappa = (sum(geo.he_length)/geo.mesh.n_he)^2 * strength;
            % id = sparse(eye(3*size(vertex, 1)));
            % P_dot = zeros(3*size(vertex, 1), 1);
            for i=1:round(duration / kappa)
                geo = Geometry(mesh, vertex);
                Hi = geo.v_mean_curvature(geo.mesh.he_src) ./ geo.v_area(geo.mesh.he_src);
                Hj = geo.v_mean_curvature(geo.mesh.he_dst) ./ geo.v_area(geo.mesh.he_dst);
                he_force = - kappa * (geo.he_schlafli_vec1 .* Hi + geo.he_schlafli_vec2 .* Hj);
                v_force = zeros(geo.mesh.n_v, 3);
                for j = 1:3
                    v_force(:, j) = accumarray(geo.mesh.he_src, he_force(:, j), [geo.mesh.n_v, 1]);
                end
                % H = geo.lap' * geo.lap ./ geo.v_area;
                % H = blkdiag(H, H, H);
                % P_dot = (H + 1 * id) \ v_force(:); 
                P_dot = v_force(:);
                vertex(:) = vertex(:) + P_dot;
                %vertex = vertex + geo.v_normal .* dot(P_dot, geo.v_normal, 2);
            end
        end
        
        function [vertex_new, mesh_new] = remesh(vertex, mesh, settings, dir, t)
            %%% default
            houdini = "./remesh.hipnc";
            %%% formatted strings
            file_in = dir + sprintf("remesh/remesh%d.mat", t);
            file_out = dir + sprintf("remesh/remesh%d.obj", t);
            matfile = sprintf('opparm -C obj/geo1/python1/ matfile ''%s''; ', file_in);
            sopoutput = sprintf('opparm -C obj/geo1/rop_geometry1/ sopoutput ''%s''; ', file_out);
            targetsize = sprintf('opparm -C obj/geo1/remesh/ targetsize %f; ', settings.edge_length);
            iterations = sprintf('opparm -C obj/geo1/remesh/ iterations %d; ', settings.n_iter);
            smoothing = sprintf('opparm -C obj/geo1/remesh/ smoothing %f; ', settings.smooth);
            call = 'render obj/geo1/rop_geometry1; quit';
            hip_command = [matfile, sopoutput, targetsize, iterations, smoothing, call];
            %%% commands
            [status, msg, msgID] = mkdir(dir + "remesh");
            save(file_in, "vertex", "mesh");
            system(sprintf('hbatch -j 12 -c "%s" %s', hip_command, houdini));
            [vertex_new,mesh_new,~,~,~,~] = readOBJ(file_out);
        end

        function [vertex, mesh] = readcsv(ptfile, primfile)
            point = readtable(ptfile);
            vertex = point{:,{ 'P_x', 'P_y', 'P_z' }};
            prim = readtable(primfile);
            mesh = prim{:,{ 'points0', 'points1', 'points2' }} + 1; % 1-indexed
            % vel_f = prim{:,{'velocity_3d_x', 'velocity_3d_y', 'velocity_3d_z'}};
            % area_f = prim{:,{'area'}};
        end
    end
end