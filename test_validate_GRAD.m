close all;
clc;
clear all;

%% set up mesh
global a c e;
a = 1;
c = 0.5;
e = sqrt(1 - c^2/a^2);
[M, P, theta, parametric_beta] = spheroid(1, 3);
geo = Geometry(M, P);
[~, K, ~, div, KTK, DTD] = geo.evolving_operators();

%% validation 1: rigid body kernel
v_area = geo.v_area;
rank = 10;
mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
mass0 = blkdiag(mass0, mass0, mass0);
[V, D] = eigs(KTK, mass0, rank, "smallestabs");

%% validation 2: mean curvature
% analytical 
geocentric_beta = para2geo(parametric_beta);
normal = normal_vector(geocentric_beta, theta);
H = mean_curvature(parametric_beta);
% numerical
Hn = reshape(div' * ones(geo.mesh.n_f, 1), [], 3) / 2;
% benchmark
eps1 = vecnorm(Hn + geo.v_mean_curvature_vec, 2, 2) ./ vecnorm(Hn, 2, 2);
eps2 = vecnorm(Hn + normal .* H .* geo.v_area, 2, 2) ./ vecnorm(Hn, 2, 2);
eps1_s = norm(Hn + geo.v_mean_curvature_vec) / norm(Hn);
eps2_s = norm(Hn + normal .* H .* geo.v_area) / norm(Hn);

%% validation 3: area conservation
% Load data
dir = "./data/willmore/ellipsoid/"; 
frame = 1:1:2000;
load(dir + "post_process.mat", "area");
dt = 5e-2;

%% plotting

% define font size 
small = 6;
medium = 6;
large = 8;
axiscolor = [0, 0, 0];

% fig = figure('Renderer', 'painters', 'Position', [10 10 1000 400]);
fig = figure;
subplot(1, 3, 1)
b = bar(abs(diag(D)));
b.FaceColor = '#ACB8B0';
set(gca,'YScale','log')
Xlabel = xlabel("$i$", "Interpreter", "latex");
Ylabel = ylabel("$|\lambda_{i}|$", "Interpreter", "latex");
set(gca, 'FontName', 'Helvetica', 'FontSize', medium);
set([Xlabel, Ylabel], 'FontName', 'Helvetica', 'FontSize', large);
set(gca, 'Box', 'on',  ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
        'XColor', axiscolor, 'YColor', axiscolor);

subplot(1, 3, 2)
analytical_beta = linspace(-pi/2, pi/2, 1000);
analytical_H = mean_curvature(analytical_beta);
plot(analytical_beta, analytical_H, "LineWidth", 0.7, 'Color', 'k');
% scatter(parametric_beta, H, 3, 'filled', 'MarkerFaceColor', 'k', 'LineWidth', 0.01);
Xlabel = xlabel("$\beta$", "Interpreter","latex", "FontSize", 24);
hold on;
scatter(parametric_beta, vecnorm(Hn, 2, 2) ./ geo.v_area, 4, "LineWidth",0.3, 'MarkerEdgeColor', '#ACB8B0');
Ylabel1 = ylabel("$H$", "Interpreter","latex");
set(gca, 'FontName', 'Helvetica', 'FontSize', medium);
% yyaxis right
% rightcolor = "#A2142F";
% ax = gca; ax.YAxis(2).Color = rightcolor;
% eps_line = scatter(parametric_beta, eps2, "filled", 'd','MarkerFaceColor', rightcolor);
% eps_line.AlphaData = 0.5;
% eps_line.MarkerFaceAlpha = 'flat';
% eps_line.MarkerFaceAlpha = 0.3;
% eps_line.MarkerEdgeAlpha = 0.3;
xticks([-pi/2, 0, pi/2]);
xticklabels({'-\pi/2','0','\pi/2'});
% Ylabel2 = ylabel("$ \epsilon_{H}$",  "Interpreter","latex");
% set(gca, 'FontName', 'Helvetica', 'FontSize', small, "YScale", "log");
set([Xlabel, Ylabel1], 'FontName', 'Helvetica', 'FontSize', large);
% leg = legend("Analytical", "Numerical", 'FontSize', small);
% leg = legend("$H(\beta)$", "$\mathrm{GRAD} \mathbf 1$","Interpreter","latex",'FontSize', small);
% leg.LineWidth = 0.01; 
% leg.BoxFace.ColorType='truecoloralpha';
% leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
% legend('boxoff');
% legend('Location','eastoutside')
set(gca, 'Box', 'on',  ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
        'XColor', axiscolor, 'YColor', axiscolor);

subplot(1, 3, 3)
line_plot = line(frame * dt, (area - area(1)) / area(1));
Xlabel = xlabel('$t$','Interpreter','latex');
Ylabel = ylabel('$(A - A_0) / A_0$','Interpreter','latex');
set(line_plot, 'LineWidth', 0.3, 'LineStyle', '-', 'Color', '#ACB8B0');
set(gca, 'FontName', 'Helvetica', 'FontSize', medium);
set([Xlabel, Ylabel], 'FontName', 'Helvetica', 'FontSize', large);
% ylim([21.4, 21.5]);
set(gca, 'Box', 'on',  ...
        'XColor', axiscolor, 'YColor', axiscolor);


% % Adjust the figure properties to ensure it doesn't get cut off
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPosition', [0 0 6 2]); % Adjust the size as needed
set(fig, 'PaperSize', [6 2]); % Adjust the size as needed
% saveas(fig, "validate_plot.png");
print(gcf, 'validate_plot.pdf', '-dpdf', '-r300','-fillpage')

%%
figure;
load(dir + "geo1.mat", "M", "P"); 
geo = Geometry(M, P);
H_elip = geo.v_mean_curvature ./geo.v_area;
IO.show(M, P, H_elip);
clim([min(H_elip), max(H_elip)]);
targetColor = [172, 184, 176]/255;
% Number of colors in the colormap
nColors = 256;
% Create the colormap transitioning from white to the target color
cmap = flip([linspace(1, targetColor(1), nColors)', ...
        linspace(1, targetColor(2), nColors)', ...
        linspace(1, targetColor(3), nColors)']);
colormap(cmap);
% colormap sky
saveas(gcf, "validate_frame1.png");
figure;
load(dir + "geo2000.mat", "M", "P"); 
geo = Geometry(M, P);
IO.show(M, P, geo.v_mean_curvature ./geo.v_area);
clim([min(H_elip), max(H_elip)]);
% colormap sky
colormap(cmap);
saveas(gcf, "validate_frame2000.png");




function [topo, coord, theta, parametric_beta] = spheroid(radius, nSub)
    global a c;
    % Obtain a spheroid mesh
    [coord, topo] = subdivided_sphere(nSub);
    coord = radius * coord;
    x = coord(:, 1);
    y = coord(:, 2);
    z = coord(:, 3);
    theta = atan2(y, x);
    r = sqrt(x.^2 + y.^2);
    % not a correct definition!!, but simply assign a parametric beta to the mesh
    parametric_beta = atan(z./r);
    coord(:, 1) = a * cos(parametric_beta) .* cos(theta);
    coord(:, 2) = a * cos(parametric_beta) .* sin(theta);
    coord(:, 3) = c * sin(parametric_beta);
end

function radius = spherical_radius(geocentric_beta)
    global a c;
    t1 = c * cos(geocentric_beta);
    t2 = a * sin(geocentric_beta);
    radius = a * c ./ sqrt(t1.^2 + t2.^2);
end

function n = normal_vector(geocentric_beta, theta)
    global a c;
    l = spherical_radius(geocentric_beta);
    r_xy = l .* cos(geocentric_beta) ./ a.^2;
    
    n = [r_xy .* cos(theta), r_xy .* sin(theta), l .* sin(geocentric_beta) ./ c.^2];
    n = n ./ vecnorm(n, 2, 2);
end

function geocentric_beta = para2geo(parametric_beta)
    global e;
    geocentric_beta = atan(tan(parametric_beta) .* sqrt(1-e^2));
end


function H = mean_curvature(parametric_beta)
    global a c
    H = c * ((2 * a^2 + (c^2 - a^2) * cos(parametric_beta).^2)) ./ ...
        (2 * a * (a^2 + (c^2 - a^2) * cos(parametric_beta).^2).^(1.5));
end

function [area, dt] = load_area(dir, frame) 
    % dir = "./data/willmore/ellipsoid/"; 
    % frame = 1:20:2000;
    
    energy = zeros(length(frame), 1);
    area = zeros(length(frame), 1);
    for i = 1:length(frame)
        load(dir + sprintf("geo%d.mat", frame(i)), ...
            "M", "P", "p"); 
        geo = Geometry(M, P);
        energy(i) = geo.willmore_energy(1);
        area(i) = geo.area;
        progressbar(i, length(frame));
    end
    dt = p.dt;
end










