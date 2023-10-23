function [ xx,yy ] = randDH2(ll, nn1, nn2, SD )
% This function generate samples xx and yy from F1 and F2
% according to DH2 configuration
temp = rand(nn1,1);
yy=(temp.^2).*(temp<=ll)./ll + temp.*(temp>ll); %from F2
temp = rand(nn2,1);
xx = temp.*(temp <= ll) + (ll + ((temp - ll).^2)./(1-ll)).*(temp>ll);
if (SD==0) 
    tmpxx=xx;
    xx=yy;
    yy=tmpxx;
end

end

