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
posX = 3;
posY = 3;
width = 15;
heigth = 12;
lw = 1.5;
s = 10;
fig = figure(1);
set(gcf,'units','centimeters','position',[posX,posY,width,heigth],'color','white');

numColumnas = size(CHI2matrix, 2);
REFparam = [kappa,dth,DeltaD,H];

vector_kappa = kappa*(1+pert);
vector_dth = dth*(1+pert);
vector_DeltaD = DeltaD*(1+pert);
vector_H = H*(1+pert);

matriz_vectores = [vector_kappa',vector_dth',vector_DeltaD',vector_H'];

param_names = {'$\kappa \, [\mathrm{h}^{-1}]$','$S_\mathrm{th} \, [-]$','$\Delta S \, [-]$','$\gamma \, [\mu \mathrm{M} \cdot \mathrm{h}^{-1}]$'};


view_angles = [150,20;60,20;60,20;60,20];
legend_position = [0.13,0.91,0.35,0.044;0.57,0.91,0.35,0.044;0.13,0.44,0.35,0.044;0.57,0.44,0.35,0.044];
%% 3D isosurfaces
for i=1:4
    subplot(2,2,i) %todos menos H
    indices = 1:4;
    indices(i) = [];
    matriz_mesh_aux = matriz_vectores(:,indices);
    names_aux = param_names(indices);
    [X,Y,Z] = meshgrid(matriz_mesh_aux(:,1),matriz_mesh_aux(:,2),matriz_mesh_aux(:,3));
     
    idx1 = 1:N;
    idx2 = 1:N;
    idx3 = 1:N;
    idx4 = 1:N;
    if i==1
        idx1 = ceil(N/2);
    elseif i==2
        idx2 = ceil(N/2);
    elseif i==3
        idx3 = ceil(N/2);
    elseif i==4
        idx4 = ceil(N/2);
    end
    
    V = CHI2matrix(idx1,idx2,idx3,idx4);
    V = reshape(V,[N,N,N]);
    Vmax = max(max(max(V)));
    Vmin = min(min(min(V)));
    
    f1 = chi2inv(0.5,8);
    f2 = chi2inv(0.7,8);
    f3 = chi2inv(0.9,8);
    
    surf1 = isosurface(X,Y,Z,V,f1); hold on;
    p1 = patch(surf1);
    isonormals(X,Y,Z,V,p1);
    set(p1,'FaceColor',[91, 60, 201]/255,'EdgeColor','none','FaceAlpha',0.3); % set the color, mesh and transparency level of the surface
    %daspect([1,1,1])
    view(3); axis tight
    camlight; lighting gouraud
    surf2=isosurface(X,Y,Z,V,f2);
    p2 = patch(surf2);
    isonormals(X,Y,Z,V,p2);
    set(p2,'FaceColor',[0.3010, 0.7450, 0.9330],'EdgeColor','none','FaceAlpha',0.2);
    surf3 = isosurface(X,Y,Z,V,f3);
    p3 = patch(surf3);
    isonormals(X,Y,Z,V,p3);
    set(p3,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.15);
    
    xlim([min(matriz_mesh_aux(:,1)),max(matriz_mesh_aux(:,1))]);
    ylim([min(matriz_mesh_aux(:,2)),max(matriz_mesh_aux(:,2))]);
    zlim([min(matriz_mesh_aux(:,3)),max(matriz_mesh_aux(:,3))]);
    xlabel(names_aux{1},'interpreter','latex','fontsize',s)
    ylabel(names_aux{2},'interpreter','latex','fontsize',s)
    zlabel(names_aux{3},'interpreter','latex','fontsize',s)

    ref = REFparam(indices); 
    plot3(ref(1),ref(2),ref(3),'x','markersize',6,'color','red');
    
    ax = gca;
    ax.FontSize = floor(0.9*s);
    ax.TickLabelInterpreter = 'latex';
    view(view_angles(i,:))
    leg = legend('$p=0.5$','$p=0.7$','$p=0.9$','interpreter','latex','fontsize',s,'numcolumns',3,'position',legend_position(i,:));
    leg.ItemTokenSize = [11,14];
end
print(fig,'contour_chi2_3d','-dpng','-r600');
