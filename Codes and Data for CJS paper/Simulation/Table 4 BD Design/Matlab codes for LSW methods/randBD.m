function [ xx,yy ] = randBD( nn1, nn2,BD )
% generate two samples according BD1-4 configurations
% Sample sizes are nn1 and nn2.
% 
  z1 = normrnd(0,1,nn1,1);  %Z_{1i}
  z2 = normrnd(0,1,nn2,1);  %Z_{2i}
%%%%%BD1  
  if (BD==1)
    mu1=0.85;
    sd1=0.6;
    mu2=0.6;
    sd2=0.8;
    xx = exp(sd1*z1+mu1);
    yy = exp(sd2*z2+mu2);
  end
%%%%%BD2  
  if (BD==2)
    mu1=0.85;
    sd1=0.6;
    mu2=1.2;
    sd2=0.2;
    xx = exp(sd1*z1+mu1);
    yy = exp(sd2*z2+mu2);
  end
%%%%%BD3  
  if (BD==3)
    mu1=0.85;
    sd1=0.6;
    mu2=0.8;
    sd2=0.5;
    mu3=0.9;
    sd3=0.9;
    z3 = normrnd(0,1,nn2,1); %% Z_{3i}
    u_i= rand(nn2,1); %u_i a uniform [0,1]
    xx = exp(sd1*z1+mu1);
    %%%%
    yy=zeros(nn2,1);
    for i = 1:nn2
        if (u_i(i)>=0.1)
            yy(i)=exp(sd2*z2(i)+mu2);
        else
            yy(i)=exp(sd3*z3(i)+mu3);
        end
    end
  end
%%%%%BD4
  if (BD==4)
    mu1=0.85;
    sd1=0.6;
    mu2=0.85;
    sd2=0.4;
    mu3=0.4;
    sd3=0.9; 
    z3 = normrnd(0,1,nn2,1); %% Z_{3i}
    u_i= rand(nn2,1); %u_i a uniform [0,1]
    xx = exp(sd1*z1+mu1);
    %%%%
     yy=zeros(nn2,1);
    for i = 1:nn2
        if (u_i(i)>=0.1)
            yy(i)=exp(sd2*z2(i)+mu2);
        else
            yy(i)=exp(sd3*z3(i)+mu3);
        end
    end
  end
end

