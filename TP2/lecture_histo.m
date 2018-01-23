function [histo] = lecture_histo(image, Nb, zone)

global Cmap

littleim = imcrop(image,zone(1:4));
littleim = rgb2ind(littleim,Cmap,'nodither');
histo = imhist(littleim,Cmap);
histo = histo/sum(histo);


end
