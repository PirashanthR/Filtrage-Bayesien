global delta rtrue  r_INS  v_INS  sigma_mes_halt X1MIN X1MAX X2MIN X2MAX N1 N2 h_ALT map update_matrix sigma_INS nb_particule

X1MIN = -10000;
X1MAX= 10000;
X2MIN = -10000;
X2MAX= 10000;
r_0 =  [-6000 -2000];
v_0 = [120 0];
delta = 1;
sigma_r0 = 100;
sigma_v0 = 10 ;
sigma_INS = 7 ;
sigma_ALT = 10;
sigma_BAR = 20 ; 

map = load('./nct/mnt.data');
[N1,N2] = size(map); %trace la carte de la map

get_pos_i_quadri= @(pos_reel)floor((pos_reel-X1MIN)/(X1MAX-X1MIN)*N1) + 1;
get_pos_j_quadri = @(pos_reel) floor((pos_reel-X2MIN)/(X2MAX-X2MIN)*N2) + 1 ;


load('./nct/traj.mat','rtrue','vtrue'); %prend rtrue et vtrue
nmax = size(rtrue,2); %nombre de pas de temps

load('./nct/ins.mat','a_INS')

r_INS(:,1) = r_0; v_INS(:,1) = v_0; % initialisation

for n=2:nmax
    r_INS(:,n) = r_INS(:,n-1)+delta*v_INS(:,n-1);
    v_INS(:,n) = v_INS(:,n-1)+delta*a_INS(:,n-1);
end

load('./nct/alt.mat','h_ALT');


%Tirage aléatoire pour commencer 

sigma_mes_halt = sqrt(sigma_ALT^2 + sigma_BAR^2);
update_matrix = [eye(2),delta*eye(2);zeros(2),eye(2)] ; 
nb_particule= 5000

% Initialialisation des donnees
delta_r_v = transpose(normrnd(0,repmat([sigma_r0 sigma_r0 sigma_v0 sigma_v0],nb_particule,1)));
w = vraisemblance(delta_r_v,1);


ponderation = w./sum(w) ;
delta_r_v_final= transpose(ponderation*transpose(delta_r_v));

traj = zeros(size(r_INS));

traj(:,1) = r_INS(:,1) + delta_r_v_final(1:2)  ;


for t=2:nmax
    [delta_r_v_choisi,~,~]= fct_multi(delta_r_v,w,nb_particule);
    delta_r_v = update_delta_r_v(delta_r_v_choisi);
    w = vraisemblance(delta_r_v,t);
    ponderation = w./sum(w) ; 
    delta_r_v_final= transpose(ponderation*transpose(delta_r_v));
    
    x_estim = r_INS(1,t) + delta_r_v_final(1);
    y_estim = r_INS(2,t) + delta_r_v_final(2);  
    
    traj(1,t) = x_estim;
    traj(2,t) = y_estim; 
    
    imagesc([X1MIN X1MAX], [X2MIN X2MAX],transpose(map));

    hold on;
    colormap('default');
    axis square;
    axis off;
    plot(rtrue(1,1:t),rtrue(2,1:t),'black')
    plot(r_INS(1,1:t),r_INS(2,1:t),'white')
    plot(traj(1,1:t),traj(2,1:t),'red')
    legend('r true','r INS', 'Position estimee SIR')
    %plot(r_INS(1,t)+delta_r_v(1,:),r_INS(2,t)+delta_r_v(2,:))
    drawnow;
    hold off
end 
