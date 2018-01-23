%%%%%%%% Lecture des images de la sequence %%%%%%%%%%

SEQUENCE = './seq1/';
START = 1;
% charge le nom des images de la sequence
filenames = dir([SEQUENCE '*.png']);
filenames = sort({filenames.name});
T = length(filenames);


%%%%%%%% Initialisation des parametres %%%%%%%%%%%%%%

global  c1 c2 c3 Cmap
N = 100
Nb = 10
lambda = 20

c1 = 800
c2 = 500
c3 = 10

%%%%%%%%%%%%%% Option affichage des particules dans l'algo sir %%%%%%%%%%%

%affichage = 0  % Afficher seulement la moyenne ponderee - rectangle rouge (l'estimateur)
%affichage = 1  % afficher la moyenne ponderee et les centres des autres particules
affichage = 2 % afficher toutes les particules et l'estimateur 


%%%%%%%%%%% Affichage de la premier image %%%%%%%%%%%%
% charge la premiere image dans 'im'
tt = START;
im = imread([SEQUENCE filenames{tt}]);
% affiche 'im'
figure;
set(gcf,'DoubleBuffer','on');
imagesc(im);

%%%%%%Selection et affichage de la zone a suivre%%%%%%%%%

disp('Cliquer 4 points dans l''image pour definir la zone a suivre.');
% on recupere la zone a tracker
zone = zeros(2,4);
compteur=1;
while(compteur ~= 5)
    [x,y,button] = ginput(1);
    zone(1,compteur) = x;
    zone(2,compteur) = y;
    text(x,y,'X','Color','r');
    compteur = compteur+1;
end
newzone = zeros(2,4);
newzone(1,:) = sort(zone(1,:));
newzone(2,:) = sort(zone(2,:));

% definition de la zone a tracker
% x haut gauche, y haut gauche, largeur, hauteur
zoneAT = zeros(1,4);
zoneAT(1) = newzone(1,1);
zoneAT(2) = newzone(2,1);
zoneAT(3) = newzone(1,4)-newzone(1,1);
zoneAT(4) = newzone(2,4)-newzone(2,1);
% affichage du rectangle � suivre 
rectangle('Position',zoneAT,'EdgeColor','r','LineWidth',3);


%%%%%% Calcul de l?histogramme de couleur associe %%%%%%%% 
% Histogramme initial
littleim = imcrop(im,zoneAT(1:4));
[tmp, Cmap] = rgb2ind(littleim, Nb, 'nodither');
littleim = rgb2ind(littleim,Cmap,'nodither');
histoRef = imhist(littleim,Cmap);
histoRef = histoRef/sum(histoRef);


%%%%%%%%%% Initialisation des particules  %%%%%%%%

X_0 = [zoneAT(1)+zoneAT(3)/2,zoneAT(2)+zoneAT(4)/2, 100];
N_particules_init = repmat(X_0,N,1) + normrnd(0,repmat([1 1 0],N,1));


%%%%%%%%%% Affichage N particules

zoneN = zeros(N,4);
zoneN(:,1) = N_particules_init(:,1)-N_particules_init(:,3)*zoneAT(3)/200;
zoneN(:,2) = N_particules_init(:,2)-N_particules_init(:,3)*zoneAT(4)/200;
zoneN(:,3) = N_particules_init(:,3)*zoneAT(3)/100;
zoneN(:,4) = N_particules_init(:,3)*zoneAT(4)/100;

for i = 1:N  
    rectangle('Position',zoneN(i,:),'EdgeColor','r','LineWidth',3);
end     

%%%%%%%% vraisemblance intiale %%%%%%%%%%

w = zeros(1,N);
for k = 1:N
        w(k) = vraisemblance(lambda,im, Nb, N_particules_init(k,:),zoneAT, histoRef) ; %
end

ponderation = w./sum(w) ;


%%%%%%%%%%%% Algorithme SIR  %%%%%%%%%%%%%%%%%%%%
N_finales = ponderation*N_particules_init;

t = size(filenames, 2); % nb d'images
w = ones(1,N);
current_vraissemblance = ones(1,N);

for j=2:(t+1)
    im = imread([SEQUENCE filenames{j-1}]); % 
    for k = 1:N
        current_vraissemblance(k) = vraisemblance(lambda,im, Nb, N_particules_init(k,:),zoneAT, histoRef) ; % Calcul des poids
    end
    
    ponderation = current_vraissemblance./sum(current_vraissemblance) ; 
    [N_choisi,~,~]= fct_multi(transpose(N_particules_init),current_vraissemblance,N);
    N_finales= ponderation*(N_particules_init); 
    N_particules_init = update_particules(transpose(N_choisi)) ;
    
    
    
if isempty(im) 
        break;
    end
    % affiche l'image et le rectangle de l'estimateur (moyenne des particules)
    imagesc(im);
    zone = zeros(1,4);
    zone(:,1) = N_finales(:,1)-N_finales(:,3)*zoneAT(3)/200;
    zone(:,2) = N_finales(:,2)-N_finales(:,3)*zoneAT(4)/200;
    zone(:,3) = N_finales(:,3)*zoneAT(3)/100;
    zone(:,4) = N_finales(:,3)*zoneAT(4)/100;

    % Affichage des centres des particules 
    switch(affichage)
        case 2
            hold on;
            for i = 1:N  
            zoneN = zeros(1,4);
            zoneN(:,1) = N_particules_init(i,1)- N_particules_init(i,3)*zoneAT(3)/200;
            zoneN(:,2) = N_particules_init(i,2)- N_particules_init(i,3)*zoneAT(4)/200;
            zoneN(:,3) = N_particules_init(i,3)*zoneAT(3)/100;
            zoneN(:,4) = N_particules_init(i,3)*zoneAT(4)/100;
            rectangle('Position',zoneN,'EdgeColor','b','LineWidth',0.1);
            end
        case 1
             % Affichage des centres des particules 
            hold on;
            plot(N_particules_init(:,1),N_particules_init(:,2),'+');
            legend('Centres des particules')
        case 0 
    end        
    %% Pour l'affichage de l'estimateur et des particules et des images
    text(10,10,num2str(j-1),'Color', 'r');
    disp(['---------- ' filenames{j-1} ' ---------------' ]);
    rectangle('Position',zone,'EdgeColor','r','LineWidth',3);

    drawnow;    
    %pause; 
    pause(0.2)
    hold off;
    
end 
