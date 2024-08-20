clear all
close all
clc

cd ..
filename1 = 'results_modelselection/results_lsqnonlin_f1_paperTea5';
filename2 = 'results_modelselection/resultsBIC_lsqnonlin_f2_UBbeta2_paperTea5';

load(filename2,'parameters');
L = parameters(1); 
n = parameters(2);  
kappa = parameters(3);
dth = parameters(4);
DeltaD = parameters(5); 
tauR = parameters(6); 
Lambda = parameters(7); 
Tmin_ini = parameters(8);
H = parameters(9);
dthA = parameters(10);
DeltaDA = parameters(11); 
Beta = parameters(12);
num = nnz(parameters);

refparam = [kappa,dth,DeltaD,H];
kappa = refparam(1);
dth = refparam(2);
DeltaD = refparam(3); 
H = refparam(4);

name = 'loop_chi2_4param_50points.mat';
load(name)
cd scripts_figures_papers

N = 50;
pert = linspace(-0.8,0.8,N); % vary parameters in a range of +-50% their optimal value

%% figure
REFparam = [kappa,dth,DeltaD,H];

vector_kappa = kappa*(1+pert);
vector_dth = dth*(1+pert);
vector_DeltaD = DeltaD*(1+pert);
vector_H = H*(1+pert);

matriz_vectores = [vector_kappa',vector_dth',vector_DeltaD',vector_H'];

param_names = {'$\kappa \, [\mathrm{h}^{-1}]$','$S_\mathrm{th} \, [-]$','$\Delta S \, [-]$','$\gamma \, [\mu \mathrm{M} \cdot \mathrm{h}^{-1}]$'};

colors = [0.3010, 0.7450, 0.9330;0.75, 0.75, 0;0.8500, 0.3250, 0.0980];

%% 2d isosurface
posX = 3;
posY = 3;
width = 15;
heigth = 10;
lw = 1.5;
s = 10;

fig = figure(3);
set(gcf,'units','centimeters','position',[posX,posY,width,heigth],'color','white');
colormap parula

combinations = nchoosek(1:4,2);
for i=1:length(combinations)
    subplot(2,3,i) %todos menos H
    indices = combinations(i,:);
    matriz_mesh_aux = matriz_vectores(:,indices);
    names_aux = param_names(indices);
    [X,Y] = meshgrid(matriz_mesh_aux(:,1),matriz_mesh_aux(:,2));
     
    idx1 = 1:N;
    idx2 = 1:N;
    idx3 = 1:N;
    idx4 = 1:N;
    
    if ~ismember(1,indices)
        idx1 = ceil(N/2);
    end
    
    if ~ismember(2,indices)
        idx2 = ceil(N/2);
    end
    
    if ~ismember(3,indices)
        idx3 = ceil(N/2);
    end
    
    if ~ismember(4,indices)
        idx4 = ceil(N/2);
    end
    
    V = CHI2matrix(idx1,idx2,idx3,idx4);
    V = reshape(V,[N,N]);
    Vmax = max(max(max(V)));
    Vmin = min(min(min(V)));
    
    f1 = chi2inv(0.5,8);
    f2 = chi2inv(0.7,8);
    f3 = chi2inv(0.9,8);
    levels = [f1,f2,f3];
    

    contour(X,Y,V,levels,'linewidth',1); hold on;
    
    xlim([min(matriz_mesh_aux(:,1)),max(matriz_mesh_aux(:,1))]);
    ylim([min(matriz_mesh_aux(:,2)),max(matriz_mesh_aux(:,2))]);
    xlabel(names_aux{1},'interpreter','latex','fontsize',s)
    ylabel(names_aux{2},'interpreter','latex','fontsize',s)

    ref = REFparam(indices); 
    plot(ref(1),ref(2),'x','markersize',5,'color','red');
    
    ax = gca;
    ax.FontSize = floor(0.9*s);
    ax.TickLabelInterpreter = 'latex';
%     leg = legend('$p=0.5$','$p=0.7$','$p=0.9$','interpreter','latex','fontsize',s);

end
print(fig,'contour_chi2_2d','-dpng','-r600');