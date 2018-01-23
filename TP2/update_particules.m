function Xnew = update_particules(Xprev) 

global c1 c2 c3

Xnew = Xprev + normrnd(0,repmat([sqrt(c1) sqrt(c2) sqrt(c3)],size(Xprev,1),1));

end