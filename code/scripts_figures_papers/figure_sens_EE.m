clear 
close all
clc

%% load experimental data
cd ..
cd ..
cd .. 
load expDATA/expDATA_paper_N
A = DATA_paperN;
ind = [2 3];
cd ajuste/Nivel1/scripts_figures_papers

mu = A.mean;    
sigma = A.std;
days = A.days;

%%
cd .. 
load results_optimal
YSOL_ref = YSOL;
cells_ref = YSOL_ref(:,1,:)+ YSOL_ref(:,2,:);
uref = cells_ref;

%% load optimal (ref) value of parameters
filename2 = 'results_modelselection/resultsBIC_lsqnonlin_f2_UBbeta2_paperTea5';

load(filename2,'parameters');

num_PoI = 8; % parameters of interest (L,r,kappa,dth,DeltaD,tauR,Lambda,H)
ind_PoI = [1,2,3,4,5,6,7,9];

filename = 'local_sens_analysis_delta10.mat';
load(filename);
cd .. 
cd Nivel1/scripts_figures_papers

%% figure
nombres = {'$L$','$r$','$\kappa$','$S_{th}$','$\Delta S$','$\tau_\mathrm{s}$','$\lambda$','$\gamma$'};
titles = {'Sensibles','Resistentes','Promedio'};
colors = {[0.4940, 0.1840, 0.5560],[0.9290, 0.6940, 0.1250],[0.4660, 0.6740, 0.1880]};
paramcolors = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],[0.75, 0.75, 0],[0, 0.5, 0],[0.75, 0, 0.75]};	
% barplot
posX = 3;
posY = 3;
width = 6;
height = 6; height = 3.5;
lw = 1.5;
s = 9;
% fig = figure(1);
% set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
for i=1:3
    fig = figure();
    set(gcf,'units','centimeters','position',[posX,posY,width,height]);
    bar(EE_norm(:,i),'FaceColor',colors{i})
    xticks(1:8);
%     title(titles{i},'interpreter','latex','fontsize',s)
    xticklabels(nombres)
    xlim([0.4,8.6]);
    ylabel('SI','fontsize',s,'interpreter','latex');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', s);
    ylim([-1.2,1.2]);
    figname = strcat('barplot',titles{i});
    print(fig,figname,'-dpng','-r600');
end


% curva EE(T)
posX = 3;
posY = 3;
width = 12;
heigth = 25;
lw = 1;
s = 9;
fig = figure();
set(gcf,'units','centimeters','position',[posX,posY,width,heigth],'color','white');
for i=1:2
    for j=1:8
        if ((i~=1)||(j~=7))
            subplot(8,2,2*j-(2-i))
            plot(tsol,cells_ref(:,1,i),'k','displayname','REF','linewidth',lw); hold on;
            sol_EE = cells_ref(:,1,i)'+EE(j,:,i)*delta*parameters(ind_PoI(j));
            plot(tsol,sol_EE,'color',paramcolors{j},'displayname','SI');
            yl = cells_ref(:,1,i)';
            yu = sol_EE;
            fill([tsol' fliplr(tsol')], [yu fliplr(yl)], paramcolors{j}, 'linestyle', 'none','FaceAlpha',0.2,'handlevisibility','off')
            leg = legend('interpreter','latex','fontsize',s,'location','northwest','NumColumns',2);
            leg.ItemTokenSize = [11,14];
%             title(nombres{j},'interpreter','latex','fontsize',s);
            ax = gca;
            ax.FontSize = floor(0.9*s);
            ax.TickLabelInterpreter = 'latex';
            ylim([0 2.3e4]);
            xlim([0 tsol(end)]);
            xlabel('$t$ [h]','interpreter','latex','fontsize',s)
            ylabel('Cell number','interpreter','latex','fontsize',s)
        end
    end   
end
print(fig,'EE_evolution','-dpng','-r600');

% a_sens = annotation('textbox',[0.23 0.92 0.3 0.05],'string','Sensitive','interpreter','latex','fontsize',1.5*s,'linestyle','none');
% a_res = annotation('textbox',[0.67 0.92 0.3 0.05],'string','Resistant','interpreter','latex','fontsize',1.5*s,'linestyle','none');