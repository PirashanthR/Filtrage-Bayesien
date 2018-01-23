X1MIN = -10000;
X1MAX= 10000;
X2MIN = -10000;
X2MAX= 10000;
r_0 =  [-6000 -2000];
v_0 = [120 0];
delta = 1

figure(1)
map = load('./nct/mnt.data');
[N1 N2] = size(map);

imagesc([X1MIN X1MAX], [X2MIN X2MAX],transpose(map));

hold on;
colormap('default');
axis square;
axis off;

load('./nct/traj.mat','rtrue','vtrue');
nmax = size(rtrue,2);


% Question 1, tracer de la position
plot(rtrue(1,:),rtrue(2,:),'r-','DisplayName','r reelle');

load('./nct/ins.mat','a_INS')

r_INS(:,1) = r_0; v_INS(:,1) = v_0; % initialisation
for n=2:nmax %evolution de la position grace a la centrale inertielle qui mesure l'acceleration
r_INS(:,n) = r_INS(:,n-1)+delta*v_INS(:,n-1); 
v_INS(:,n) = v_INS(:,n-1)+delta*a_INS(:,n-1);
end

%Question 2 tracer de la position estimée
%plot(r_INS(1,:),r_INS(2,:),'m-','DisplayName','r mesure');
legend('show')
title('Trajectoire de l avion')
xlabel('x')
ylabel('y')

hold off 
%Definition des fonctions qui permettent d'associer a une
%position ses indices dans la matrice map associe
get_pos_i_quadri= @(pos_reel)floor((pos_reel-X1MIN)/(X1MAX-X1MIN)*N1) + 1;
get_pos_j_quadri = @(pos_reel) floor((pos_reel-X2MIN)/(X2MAX-X2MIN)*N2) + 1 ;


%La liste des position r_t = x,y en terme d'indices dans le maillage 
pos_i = get_pos_i_quadri(rtrue(1,:));
pos_j = get_pos_j_quadri(rtrue(2,:));

%Les hauteurs reelles
h_true = arrayfun(@(i,j)map(i,j),pos_i,pos_j);

figure(2)
%On trace la vraie position sur un graphie
plot(h_true,'blue','DisplayName','h true')
hold on

load('./nct/alt.mat','h_ALT');

%On trace la position mesuree sur le meme graphique
plot(h_ALT, 'r+','DisplayName','h mesure');
legend('show')
title('Comparaison hauteurs reelles et hauteurs mesurees')
ylabel('Hauteur')
xlabel('Temps')
hold off

