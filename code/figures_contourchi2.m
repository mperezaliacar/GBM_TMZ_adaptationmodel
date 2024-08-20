function figures_contourchi2(name,refparam,N)

pert = linspace(-0.8,0.8,N); % vary parameters in a range of +-50% their optimal value

kappa = refparam(1);
dth = refparam(2);
DeltaD = refparam(3); 
H = refparam(4);

load(name)

%% figure
posX = 3;
posY = 3;
width = 20;
heigth = 18;
lw = 1.5;
s = 15;
fig = figure(1);
set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);

numColumnas = size(CHI2matrix, 2);
REFparam = [kappa,dth,DeltaD,H];

vector_kappa = kappa*(1+pert);
vector_dth = dth*(1+pert);
vector_DeltaD = DeltaD*(1+pert);
vector_H = H*(1+pert);

matriz_vectores = [vector_kappa',vector_dth',vector_DeltaD',vector_H'];

param_names = {'$\kappa$','$S_\mathrm{th}$','$\Delta S$','$\gamma$'};

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
    set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
    view(3); axis tight
    camlight; lighting gouraud
    surf2=isosurface(X,Y,Z,V,f2);
    p2 = patch(surf2);
    isonormals(X,Y,Z,V,p2);
    set(p2,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2);
    surf3 = isosurface(X,Y,Z,V,f3);
    p3 = patch(surf3);
    isonormals(X,Y,Z,V,p3);
    set(p3,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.3);
    
    xlim([min(matriz_mesh_aux(:,1)),max(matriz_mesh_aux(:,1))]);
    ylim([min(matriz_mesh_aux(:,2)),max(matriz_mesh_aux(:,2))]);
    zlim([min(matriz_mesh_aux(:,3)),max(matriz_mesh_aux(:,3))]);
    xlabel(names_aux{1},'interpreter','latex','fontsize',s)
    ylabel(names_aux{2},'interpreter','latex','fontsize',s)
    zlabel(names_aux{3},'interpreter','latex','fontsize',s)

    ref = REFparam(indices); 
    plot3(ref(1),ref(2),ref(3),'*','markersize',10,'color','red');
end
figname = strcat('3d_isosurface_',num2str(N),'points.fig');
savefig(fig,figname);

%% 2d isosurface
posX = 3;
posY = 3;
width = 20;
heigth = 18;
lw = 1.5;
s = 15;
fig = figure(3);
set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);

combinations = nchoosek(1:4,2);
for i=1:length(combinations)
    subplot(3,2,i) %todos menos H
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
    
    contour(X,Y,V,levels); hold on;
    
    xlim([min(matriz_mesh_aux(:,1)),max(matriz_mesh_aux(:,1))]);
    ylim([min(matriz_mesh_aux(:,2)),max(matriz_mesh_aux(:,2))]);
    xlabel(names_aux{1},'interpreter','latex','fontsize',s)
    ylabel(names_aux{2},'interpreter','latex','fontsize',s)

    ref = REFparam(indices); 
    plot(ref(1),ref(2),'*','markersize',10,'color','red');
end

figname = strcat('2d_isosurface_',num2str(N),'points.fig');
savefig(fig,figname);
