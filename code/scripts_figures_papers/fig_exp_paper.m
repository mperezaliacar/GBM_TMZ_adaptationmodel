clear all
close all
clc

%% figures CMBBE
% colors
TMZ = [0.4940, 0.1840, 0.5560];	
TMZ_BIG = [0.9290, 0.6940, 0.1250];
DMSO_BIG = [0, 0.4470, 0.7410];
color_marker = [250, 55, 243]/255;

cd .. 
cd ..
cd ..
load expDATA/expDATA_paper_area
cd ajuste/Nivel1/scripts_figures_papers

posX = 10;
posY = 10;
width = 11; width = 8; % for the paper (2 columns)
heigth = 7;
lw = 1.5;
s = 10;
fig = figure(1);
set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);

x = DATA_paper.days;
y = DATA_paper.mean;
err = DATA_paper.std;
IND = [1 2 3]; 

colors = {DMSO_BIG,TMZ,TMZ_BIG};
%colors = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250]};

name = {'Control','TMZ-sensitive','TMZ-resistant'};

for i=1:3
    errorbar(x,y(:,IND(i)),err(:,IND(i)),'color',colors{i},'linewidth',lw,'DisplayName',name{i}); hold on;
end

xx = [0 1 2 3 4 28 29 30 31 32];
yy = y(1,1)*ones(size(xx));
plot(xx(1),yy(1),'o','MarkerSize',4,'MarkerEdgeColor',color_marker,'MarkerFaceColor',color_marker,'DisplayName','TMZ dose');
plot(xx(2:(end)),yy(2:(end)),'o','MarkerSize',4,'MarkerEdgeColor',color_marker,'MarkerFaceColor',color_marker,'HandleVisibility','off');

xlim([0 1.01*x(end)]);

ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';

% ingles
xlabel('Days','interpreter','latex','fontsize',s)
ylabel('Growth rate [$\%$]','interpreter','latex','fontsize',s)


leg = legend('interpreter','latex','fontsize',s,'location','northwest','Box','on');
leg.ItemTokenSize = [11 16];

print(fig,'Fig_results_exp','-dpng','-r600')

%% figures scheme

% posX = 10;
% posY = 10;
% width = 6;
% heigth = 5;
% lw = 2;
% s = 12;
% fig = figure(2);
% set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
% 
% x = DATA_clinicALL.days;
% y = DATA_clinicALL.mean;
% err = DATA_clinicALL.std;
% IND = [1 2 3]; 
% 
% colors = {DMSO_BIG,TMZ,TMZ_BIG};
% 
% %name = {'Control (27)','TMZ-sensibles (23)','TMZ-resistentes (8)'};
% 
% for i=1:3
%     fig = figure(i);
%     set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
%     errorbar(x,y(:,IND(i)),err(:,IND(i)),'color',colors{i},'linewidth',lw); hold on;
%     set(gca, 'Ytick', [], 'box', 'off');
%     set(gca, 'Xtick', [], 'box', 'off');
%     set(gca,'LineWidth',0.5*lw);
% end
% 
% xlim([0 1.01*x(end)]);
% 
% %print(fig,'Exp_español','-dpng','-r600')