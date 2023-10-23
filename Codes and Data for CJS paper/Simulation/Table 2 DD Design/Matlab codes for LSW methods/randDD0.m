function [ xx,yy ] = randDD0( seg, nn1, nn2, SD )
% generate two samples: both taking values in the unit interval
% Sample sizes are nn1 and nn2.
% G is uniform and distribution of YY.
% It decides the percentage of observations
% falls into 8 segments
      
if (SD==1)
    mm = mnrnd(nn1, seg,1);
    xx = zeros(nn1,1);
    tmpa=0; %#ok<*NASGU>
    tmpb=0;
    for i =1:8
        tmpa=tmpb+1;
        tmpb=tmpb+mm(i);
        xx(tmpa:tmpb) = (i-1+rand(mm(i),1))/8;
    end
    yy = rand(nn2,1);
end
if (SD==0)
    mm = mnrnd(nn2, seg,1);
    yy = zeros(nn2,1);
    tmpa=0; %#ok<*NASGU>
    tmpb=0;
    for i =1:8
        tmpa=tmpb+1;
        tmpb=tmpb+mm(i);
        yy(tmpa:tmpb) = (i-1+rand(mm(i),1))/8;
    end
    xx = rand(nn1,1);
end
end

