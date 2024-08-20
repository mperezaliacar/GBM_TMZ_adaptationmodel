clear 
close all
clc

cd .. 
cd ..
cd ..
load expDATA/expDATA

x1 = DATA_clinical3.days;
y1 = DATA_clinical3.mean(:,5);
err1 = DATA_clinical3.std(:,5);

load expDATA/expDATA_paper_N

x2 = DATA_paperN.days;
y2 = DATA_paperN.mean(:,3);
err2 = DATA_paperN.std(:,3);

%% figures CMBBE
% colors
TMZ = [0, 0.4470, 0.7410];	
TMZ_BIG = [0.4940, 0.1840, 0.5560];
DMSO_BIG = [0.4660, 0.6740, 0.1880];
color_marker = [250, 55, 243]/255;



cd ajuste/Nivel1/scripts_figures_papers
posX = 10;
posY = 10;
width = 11; width = 8; % for paper (2 columns)
heigth = 6;
lw = 1.5;
s = 9;

fig = figure(1);
set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);

colors = {[237, 141, 31]/255, [0.9290, 0.6940, 0.1250]};

name = {'2 cycles','1 cycle'};

errorbar(x1,y1,err1,'color',colors{1},'linewidth',lw,'DisplayName','TMZ-resistant 1 cycle'); hold on;

xx = [0 1 2 3 4 ];
yy = y1(1,1)/2*ones(size(xx));
plot(xx(1),yy(1),'o','MarkerSize',4,'MarkerEdgeColor',color_marker,'MarkerFaceColor',color_marker,'DisplayName','TMZ dose');
plot(xx(2:(end)),yy(2:(end)),'o','MarkerSize',4,'MarkerEdgeColor',color_marker,'MarkerFaceColor',color_marker,'HandleVisibility','off');

xlim([0 1.01*x1(end)]);

ylabel('Growth rate [$\%$]','interpreter','latex','fontsize',s)
xlabel('Days','interpreter','latex','fontsize',s)

ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';

leg = legend('interpreter','latex','fontsize',s,'location','northwest');
leg.ItemTokenSize = [11 16];

print(fig,'exp_1cycle','-dpng','-r600')
