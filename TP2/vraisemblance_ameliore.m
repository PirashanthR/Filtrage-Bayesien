function [V] = vraisemblance_ameliore(lambda,ImageNew, Nb, N_particule_init,im_prev,N_particules_prev,zoneRef,histoRef)


zoneNew = zeros(1,4);
zoneNew(:,1) = N_particule_init(:,1)-N_particule_init(:,3)*zoneRef(3)/200;
zoneNew(:,2) = N_particule_init(:,2)-N_particule_init(:,3)*zoneRef(4)/200;
zoneNew(:,3) = N_particule_init(:,3)*zoneRef(3)/100;
zoneNew(:,4) = N_particule_init(:,3)*zoneRef(4)/100;

if N_particules_prev~=0 
    zonePrev = zeros(1,4);
    zonePrev(:,1) = N_particules_prev(:,1)-N_particules_prev(:,3)*zoneRef(3)/200;
    zonePrev(:,2) = N_particules_prev(:,2)-N_particules_prev(:,3)*zoneRef(4)/200;
    zonePrev(:,3) = N_particules_prev(:,3)*zoneRef(3)/100;
    zonePrev(:,4) = N_particules_prev(:,3)*zoneRef(4)/100;
else
    zonePrev = zoneNew;
end


%Test si le rectangle sort
x = size(ImageNew,2);
y = size(ImageNew,1);

if((zoneNew(1)+zoneNew(3)>x | zoneNew(2)+zoneNew(4) > y | zoneNew(2) < 0 | zoneNew(1) < 0 )|(zonePrev(1)+zonePrev(3)>x | zonePrev(2)+zonePrev(4) > y | zonePrev(2) < 0 | zonePrev(1)<0 ))
    V = 0;
else
    histoNew = lecture_histo(ImageNew,Nb,zoneNew);
    histoPrev = lecture_histo(im_prev,Nb,zonePrev);
    DistanceN = sqrt(1 - sum(sqrt(histoNew.*histoRef)));
    DistanceP = sqrt(max(0,1 - sum(sqrt(histoPrev.*histoNew))));
    Distance = DistanceN + DistanceP;
    V = exp(-lambda*Distance^2);
end

end