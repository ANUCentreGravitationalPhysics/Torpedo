function z = gaussbeam(dists,offset,waist)


z = waist*sqrt(1+((dists-offset)./(pi*waist^2/(1064*10^(-9)))).^2);