clear 
close all
clc

cd ..
load data_media_betalambda_1
cd ..
cd ..
load expDATA/expDATA_ind
name = 'paperTea';
cd ajuste/Nivel1/scripts_figures_papers

% posX = 3;
% posY = 3;
% width = 8;
% heigth = 6;
% lw = 1.5;
% s = 9;

% fig1 = figure();
% set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
% histogram(vectorLambda,10);
% ax = gca;
% ax.FontSize = floor(0.9*s);
% ax.TickLabelInterpreter = 'latex';
% xlabel('$\lambda$','interpreter','latex','fontsize',s);
% ylim([0 90])
% print(fig1,'hist_lambda','-dpng','-r600');
% 
% fig2 = figure();
% set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
% histogram(vectorBeta,10);
% ax = gca;
% ax.FontSize = floor(0.9*s);
% ax.TickLabelInterpreter = 'latex';
% xlabel('$\beta$','interpreter','latex','fontsize',s);
% ylim([0 90])
% print(fig2,'hist_beta','-dpng','-r600');

%% kmeans
rng(107)
idx_true = zeros(size(datos_tratados,1),1);
idx_true(1:69) = 1;
idx_true(70:end) = 2;

% transformation option 1 (maxmin)
% transVectorBeta = (vectorBeta-min(vectorBeta))/(max(vectorBeta)-min(vectorBeta));
% transVectorLambda = (vectorLambda-min(vectorLambda))/(max(vectorLambda)-min(vectorLambda));

% transformation option 2 (meanstd) - SELECTED
transVectorBeta = (vectorBeta-mean(vectorBeta))/std(vectorBeta);
transVectorLambda = (vectorLambda-mean(vectorLambda))/std(vectorLambda);

X = [transVectorLambda, transVectorBeta];
opts = statset('Display','final');
[idx,C] = kmeans(X,2,'replicates',10,'options',opts);
CM = confusionmat(idx_true,idx);

% figure kmeans
posX = 3;
posY = 3;
width = 8;
heigth = 8;

colorSENS = [0.9290, 0.6940, 0.1250];
colorRES = [0.4940, 0.1840, 0.5560];
fig3 = set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
plot(X(idx==1,1),X(idx==1,2),'.','color',colorSENS,'MarkerSize',10)
hold on
plot(X(idx==2,1),X(idx==2,2),'diamond','color',colorRES,'MarkerSize',3,'MarkerFaceColor',colorRES)
plot(C(:,1),C(:,2),'kx','MarkerSize',10,'LineWidth',2) 
legend('Cluster 1 (TMZ-S)','Cluster 2 (TMZ-R)','Centroids','Location','NW','interpreter','latex','fontsize',s)
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
xlabel('$\hat{\lambda}$','interpreter','latex','fontsize',s);
ylabel('$\hat{\beta}$','interpreter','latex','fontsize',s);
ylim([-1 8.5]);
xlim([-0.8 2.3]);

print(fig3,'kmeans_betalambda','-dpng','-r600');

% figure profiles
posX = 3;
posY = 3;
width = 12;
heigth = 14;
lw = 0.5;
datos_resistentes_bien=datos_resistentes;
datos_resistentes_bien(21,:) = [];
figure()
fig4 = set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
plot(time, datos_sensibles,'-o','linewidth',lw,'handlevisibility','off','color',colorSENS,'markersize',5,'markeredgecolor',colorSENS); hold on;
plot(time, datos_resistentes_bien,'-diamond','linewidth',lw,'handlevisibility','off','color',colorRES,'markersize',5,'markeredgecolor',colorRES,'MarkerFaceColor',colorRES);
plot(time, datos_tratados(90,:),'-diamond','linewidth',lw,'handlevisibility','off','color',colorSENS,'markersize',5,'markeredgecolor',colorRES,'MarkerFaceColor',colorRES);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
xlabel('$t$ [days]','interpreter','latex','fontsize',s);
ylabel('Cell number','interpreter','latex','fontsize',s);
xlim([0 time(end)]);
ylim([0 5.5e4]);
