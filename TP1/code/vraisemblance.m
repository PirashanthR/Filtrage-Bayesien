function value=vraisemblance(pos_particule,time_step)

global r_INS  sigma_mes_halt X1MIN X1MAX X2MIN X2MAX N1 N2 h_ALT

get_pos_i_quadri= @(pos_reel)floor((pos_reel-X1MIN)/(X1MAX-X1MIN)*N1) + 1;
get_pos_j_quadri = @(pos_reel) floor((pos_reel-X2MIN)/(X2MAX-X2MIN)*N2) + 1 ;
%donc vraisemblance nulle

pos_i = get_pos_i_quadri(pos_particule(1,:)+r_INS(1,time_step));
pos_j = get_pos_j_quadri(pos_particule(2,:)+r_INS(2,time_step));


value = exp(-1/2*(h_ALT(time_step)  - arrayfun(@(i,j)read_map(i,j),pos_i,pos_j)).^2./sigma_mes_halt^2);
end 
