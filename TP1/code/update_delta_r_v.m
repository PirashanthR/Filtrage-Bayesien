function updated_matrix = update_delta_r_v(delta_r_v_previous) 

global update_matrix sigma_INS nb_particule delta

updated_matrix = update_matrix*delta_r_v_previous - delta*transpose(normrnd(0,repmat([0 0 sigma_INS,sigma_INS],nb_particule,1)));
end
