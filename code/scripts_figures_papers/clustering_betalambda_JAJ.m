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
s = 9;
% 
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

% transformation (meanstd) 
transVectorBeta = (vectorBeta-mean(vectorBeta))/std(vectorBeta);
transVectorLambda = (vectorLambda-mean(vectorLambda))/std(vectorLambda);

X = [transVectorLambda, transVectorBeta];
opts = statset('Display','final');
[idx,C] = kmeans(X,2,'replicates',10,'options',opts);
CM = confusionmat(idx_true,idx);

% indices by categories
realS = find(idx_true == 1);
realR = find(idx_true == 2);
predS = find(idx == 1);
predR = find(idx == 2);

pSrS = find(idx_true == 1 & idx == 1); % predicted sensitive, real sensitive
pRrR = find(idx_true == 2 & idx == 2); % predicted resistant, real resistant
pSrR = find(idx_true == 2 & idx == 1); % predicted sensitive, real resistant
pRrS = find(idx_true == 1 & idx == 2); % predicted resistant, real sensitive


%%
% figure kmeans
posX = 3;
posY = 3;
width = 6.5;
heigth = 8;

colorRES = [0.9290, 0.6940, 0.1250];
colorSENS = [0.4940, 0.1840, 0.5560];
fig3 = set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
plot(X(pSrS,1),X(pSrS,2),'o','color',colorSENS,'MarkerSize',3,'MarkerFaceColor',colorSENS)
hold on
plot(X(pRrR,1),X(pRrR,2),'diamond','color',colorRES,'MarkerSize',3,'MarkerFaceColor',colorRES)
plot(X(pSrR,1),X(pSrR,2),'o','color',colorRES,'MarkerSize',3,'MarkerFaceColor',colorRES)
plot(C(1,1),C(1,2),'x','color',[60,4,74]/255,'MarkerSize',8,'LineWidth',1.5) 
plot(C(2,1),C(2,2),'x','color',[156, 115, 3]/255,'MarkerSize',8,'LineWidth',1.5) 
% legend('Cluster 1 (TMZ-S)','Cluster 2 (TMZ-R)','Centroids','Location','NW','interpreter','latex','fontsize',s)
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
xlabel('$\hat{\lambda}$','interpreter','latex','fontsize',s);
ylabel('$\hat{\beta}$','interpreter','latex','fontsize',s);
ylim([-1 8.5]);
xlim([-0.8 2.3]);

print(fig3,'kmeans_betalambda','-dpng','-r600');

%% figure profiles
posX = 3;
posY = 3;
width = 11;
heigth = 13;
lw = 0.5;

% figure()
% fig4 = set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
% semilogy(time, datos_sensibles_bien,':o','linewidth',lw,'handlevisibility','off','color',colorSENS,'markersize',5,'markeredgecolor',colorSENS); hold on;
% %semilogy(time, datos_sensibles_mal,':d','linewidth',lw,'handlevisibility','off','color',colorSENS,'markersize',5,'markeredgecolor',colorSENS);
% semilogy(time, datos_resistentes_bien,':d','linewidth',lw,'handlevisibility','off','color',colorRES,'markersize',5,'markeredgecolor',colorRES);
% semilogy(time, datos_resistentes_mal,':o','linewidth',lw,'handlevisibility','off','color',colorRES,'markersize',5,'markeredgecolor',colorRES);
% ax = gca; grid on;
% ax.FontSize = floor(0.9*s);
% ax.TickLabelInterpreter = 'latex';
% xlabel('$t$ [days]','interpreter','latex','fontsize',s);
% ylabel('Cell number','interpreter','latex','fontsize',s);
% xlim([0 time(end)]);
% ylim([0 5.5e4]);
% 
% %% figure profiles
% figure()
% fig4 = set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
% semilogy(time, datos_sensibles_bien,':o','linewidth',lw,'handlevisibility','off','color','b','markersize',5,'markeredgecolor',colorSENS,'MarkerFacecolor',colorSENS); hold on;
% %semilogy(time, datos_sensibles_mal,':d','linewidth',lw,'handlevisibility','off','color',colorSENS,'markersize',5,'markeredgecolor',colorSENS,'MarkerFacecolor',colorSENS);
% semilogy(time, datos_resistentes_bien,':d','linewidth',lw,'handlevisibility','off','color','b','markersize',5,'markeredgecolor',colorRES,'MarkerFacecolor',colorRES);
% semilogy(time, datos_resistentes_mal,':o','linewidth',lw,'handlevisibility','off','color','r','markersize',5,'markeredgecolor',colorRES,'MarkerFacecolor',colorRES);
% ax = gca; grid on;
% ax.FontSize = floor(0.9*s);
% ax.TickLabelInterpreter = 'latex';
% xlabel('$t$ [days]','interpreter','latex','fontsize',s);
% ylabel('Cell number','interpreter','latex','fontsize',s);
% xlim([0 time(end)]);
% ylim([0 5.5e4]);

%% figure profiles
figure()
colorBIEN = [0.4667    0.6745    0.1882];
colorMAL = [1 0 0];
fig4 = set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
semilogy(time, datos_tratados(pSrS,:),':o','linewidth',lw,'handlevisibility','off','color',colorBIEN,'markersize',5,'markeredgecolor',colorSENS,'MarkerFacecolor',colorSENS); hold on;
%semilogy(time, datos_tratados(pRrS),':d','linewidth',lw,'handlevisibility','off','color',colorMAL,'markersize',5,'markeredgecolor',colorSENS,'MarkerFacecolor',colorSENS);
semilogy(time, datos_tratados(pRrR,:),':d','linewidth',lw,'handlevisibility','off','color',colorBIEN,'markersize',5,'markeredgecolor',colorRES,'MarkerFacecolor',colorRES);
semilogy(time, datos_tratados(pSrR,:),'-o','linewidth',lw,'handlevisibility','off','color',colorMAL,'markersize',5,'markeredgecolor',colorRES,'MarkerFacecolor',colorRES);
ax = gca; grid on;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
xlabel('$t$ [days]','interpreter','latex','fontsize',s);
ylabel('Cell number','interpreter','latex','fontsize',s);
xlim([0 time(end)]);
ylim([0 5.5e4]);

print(fig4,'kmeans_profiles','-dpng','-r600');

%% kmeans for time series
num = 1000;
errores = zeros(1,num);
for i=1:num
    rng(i)
    [idx_t, C_t] = kmeans(datos_tratados,2);

    % figure()
    % semilogy(time, datos_tratados(idx_t == 1,:),'-','linewidth',lw,'handlevisibility','off','color',colorSENS); hold on;
    % semilogy(time, datos_tratados(idx_t == 2,:),'-','linewidth',lw,'handlevisibility','off','color',colorRES);
    CM_t = confusionmat(idx_true,idx_t);
    diag1 = CM_t(1,1)+CM_t(2,2);
    diag2 = CM_t(1,2)+CM_t(2,1);
    errores(i) = min(diag1,diag2);
end
hist(errores,20)