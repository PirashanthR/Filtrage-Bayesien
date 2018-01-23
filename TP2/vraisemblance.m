function [V] = vraisemblance (lambda,ImageNew, Nb, N_particule_init,zoneRef,histoRef)


zoneNew = zeros(1,4);
zoneNew(:,1) = N_particule_init(:,1)- N_particule_init(:,3)*zoneRef(3)/200;
zoneNew(:,2) = N_particule_init(:,2)- N_particule_init(:,3)*zoneRef(4)/200;
zoneNew(:,3) = N_particule_init(:,3)*zoneRef(3)/100;
zoneNew(:,4) = N_particule_init(:,3)*zoneRef(4)/100;

%Test si le rectangle sort
x = size(ImageNew,2);
y = size(ImageNew,1);

if((zoneNew(1)+zoneNew(3) < x+1) && (zoneNew(2)+zoneNew(4) < y+1) && (zoneNew(2)> 0) && (zoneNew(1) > 0) )
    
    histoNew = lecture_histo(ImageNew,Nb,zoneNew);
    Distance = sqrt(1 - sum(sqrt(histoNew.*histoRef)));
    V = exp(-lambda*Distance^2);
else
    V = 0;
end

end