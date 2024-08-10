close all;
clear all;
clc;

% Load data
% dir = "./data/willmore/genus6/"; 
% dir = "./data/willmore/ellipsoid/"; 
dir = "./data/bulk/sphere/strainrot_alpha5e-1_kappa1e-2/";
gap = 1;
frame = 1:gap:400;
% 
% dissipation = zeros(length(frame), 1);
% energy = zeros(length(frame), 1);
% area = zeros(length(frame), 1);
% work = zeros(length(frame), 1);
% for i = 1:length(frame)
%     data = load(dir + sprintf("geometry/geo%d.mat", frame(i)));
%     geo = Geometry(data.M, data.P);
%     mass0 = spdiags(geo.v_area, 0, geo.mesh.n_v, geo.mesh.n_v);
%     mass0 = blkdiag(mass0, mass0, mass0);
%     [~, K, ~, div, KTK, DTD] = geo.evolving_operators();
% 
%     dissipation(i) = data.velocity' * KTK * data.velocity * data.p.dt;
%     energy(i) = geo.willmore_energy(data.p.kappa);
%     force = data.p.alpha * (data.velocity - data.u0);
%     % force = data.fb;
%     work(i) = data.velocity' * mass0 * force * data.p.dt;
%     area(i) = geo.area;
%     progressbar(i, length(frame));
% end
% save(dir + "post_process.mat", "energy", "area", "dissipation", "work");

%% plotting 1
load(dir + "post_process.mat", "energy", "area", "dissipation","work");
work = smoothdata(work, 'movmean', 15);
bending_work = [0; smoothdata(diff(energy), 'movmean', 15)];
dissipation = smoothdata(dissipation, 'movmean', 15);
% define font size 
small = 6;
medium = 6;
large = 8;
axiscolor = [0, 0, 0];

% fig = figure('Renderer', 'painters', 'Position', [10 10 1000 400]);
fig = figure;
subplot(1, 3, 3)
hold on;     
dt = 0.01;        
time = frame * dt;
l4 = plot(time, 2 / dt * dissipation(frame));
l5 = plot(time, (work(frame))/dt);
l6 = plot(time, bending_work(frame)/dt);
% plot(time(1:end-1), smoothdata(diff(energy(frame)) + work(frame(1:end-1)) + dissipation(frame(1:end-1))));
% plot(time, energy(frame));
% plot(time, 2 * cumsum(dissipation(frame)) * gap);
% plot(time, cumsum(work(frame)) * gap);
% plot(time, energy(frame) + 2 * cumsum(dissipation(frame)) * gap + cumsum(work(frame)) * gap);
% legend("$V$", "dissipation", "work", "sum");
% legend("$\mathbf T$", "$ \mathbf B_{\alpha}$ ", "$ \mathbf B_{\kappa} $", 'Interpreter','latex');
gap = 1000;frame = 50:gap:400;time = frame * dt;
l1 = scatter(time, 2 / dt * dissipation(frame), 15, "filled", "d", "LineWidth",0.6,'MarkerEdgeColor',  'k', 'MarkerFaceColor',  '#ACB8B0');
l2 = scatter(time, work(frame)/dt, 15, 's', "LineWidth",0.6,'MarkerEdgeColor', 'k','MarkerFaceColor',  '#ACB8B0');
l3 = scatter(time, bending_work(frame)/dt, 25, "pentagram", "filled", "LineWidth",0.6,'MarkerEdgeColor', 'k','MarkerFaceColor',  '#ACB8B0');
set([l4, l5, l6], 'LineWidth', 1);
set([l4, l5, l6], 'Color', '#ACB8B0');
Xlabel = xlabel('$t$','Interpreter','latex');
Ylabel = ylabel('$\dot E$','Interpreter','latex');
% line_plot = line(frame / 100, area);
set(gca, 'FontName', 'Helvetica', 'FontSize', medium);
set([Xlabel, Ylabel], 'FontName', 'Helvetica', 'FontSize', large);
set(gca, 'Box', 'on',  ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
        'XColor', axiscolor, 'YColor', axiscolor);
% leg = legend("", "", "", "$\mathbf T$", "$ \mathbf B_{\alpha}$ ", "$ \mathbf B_{\kappa} $",...
%      'Interpreter','latex','FontSize', small, 'Orientation','horizontal');
% legend('boxoff');
% legend('Location','northoutside')

% % Adjust the figure properties to ensure it doesn't get cut off
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPosition', [0 0 6 2]); % Adjust the size as needed
set(fig, 'PaperSize', [6 2]); % Adjust the size as needed
% saveas(fig, "validate_plot.png");
print(gcf, 'bulk_plot.pdf', '-dpdf', '-r300','-fillpage')


%% plotting 2 
dir = "./data/willmore/genus6/"; 
gap = 1;
frame = 1:gap:4000;
load(dir + "post_process.mat", "energy", "area", "dissipation");
dissipation = 2 * cumsum(dissipation(frame));

% define font size 
small = 6;
medium = 6;
large = 8;
axiscolor = [0, 0, 0];

% fig = figure('Renderer', 'painters', 'Position', [10 10 1000 400]);
fig = figure;
subplot(1, 3, 3)
hold on;     
dt = 0.01;        
time = frame * dt;
l4 = plot(time, energy(frame));
l5 = plot(time, dissipation(frame));
% plot(time(1:end-1), smoothdata(diff(energy(frame)) + work(frame(1:end-1)) + dissipation(frame(1:end-1))));
% plot(time, energy(frame));
% plot(time, 2 * cumsum(dissipation(frame)) * gap);
% plot(time, cumsum(work(frame)) * gap);
% plot(time, energy(frame) + 2 * cumsum(dissipation(frame)) * gap + cumsum(work(frame)) * gap);
% legend("$V$", "dissipation", "work", "sum");
% legend("$\mathbf T$", "$ \mathbf B_{\alpha}$ ", "$ \mathbf B_{\kappa} $", 'Interpreter','latex');
gap = 10000;frame = 3000:gap:4000;time = frame * dt;
l1 = scatter(time, energy(frame), 15, "filled", "d", "LineWidth",0.6,'MarkerEdgeColor',  'k', 'MarkerFaceColor',  '#ACB8B0');
l2 = scatter(time, dissipation(frame), 15, 's', "LineWidth",0.6,'MarkerEdgeColor', 'k','MarkerFaceColor',  '#ACB8B0');
set([l4, l5], 'LineWidth', 1);
set([l4, l5], 'Color', '#ACB8B0');
Xlabel = xlabel('$t$','Interpreter','latex');
Ylabel = ylabel('$E$','Interpreter','latex');
% line_plot = line(frame / 100, area);
set(gca, 'FontName', 'Helvetica', 'FontSize', medium);
set([Xlabel, Ylabel], 'FontName', 'Helvetica', 'FontSize', large);
set(gca, 'Box', 'on',  ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
        'XColor', axiscolor, 'YColor', axiscolor);
% leg = legend("", "", "$V$", "$E_{\mu}$ ",...
%  'Interpreter','latex','FontSize', small, 'Orientation','horizontal');
% legend('boxoff');
% legend('Location','northoutside')

% % Adjust the figure properties to ensure it doesn't get cut off
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPosition', [0 0 6 2]); % Adjust the size as needed
set(fig, 'PaperSize', [6 2]); % Adjust the size as needed
% saveas(fig, "validate_plot.png");
print(gcf, 'genus6_plot.pdf', '-dpdf', '-r300','-fillpage')

