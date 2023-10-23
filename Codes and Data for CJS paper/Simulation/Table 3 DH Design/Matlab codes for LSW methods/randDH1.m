function [ xx,yy ] = randDH1( ll, nn1, nn2, SD )
% This function generate samples xx and yy from F1 and F2
% defined by Whang 2019, page85,design 1 
temp = rand(nn1,1);
yy=(temp.^2).*(temp<=ll)./ll + temp.*(temp>ll);
xx = rand(nn2,1);
if (SD==0) 
    tmpxx=xx;
    xx=yy;
    yy=tmpxx;
end
end

