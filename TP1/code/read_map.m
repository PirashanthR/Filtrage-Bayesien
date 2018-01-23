function hauteur = read_map(i,j)

global map 

if (i>0)&&(j>0)&&(i<469)&&(j<359)
    hauteur = map(i,j);
else 
    hauteur = Inf(1);
end 
