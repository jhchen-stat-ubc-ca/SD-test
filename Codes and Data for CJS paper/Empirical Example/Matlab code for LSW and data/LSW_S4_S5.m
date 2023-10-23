tic;
%-------------------------- Input ---------------------------%

s1 = csvread('S4.csv');  %sample 1 data for 1983
s2 = csvread('S5.csv');  %sample 2 data for 1985



B=11;
SDorder=1;
ngrid=50;
c=3;


%n1=200;
%n2=200;

% Loop the whole thing and record pvalues n number of times
% Store pvalues in a n x 5 matrix
n = 500;
allstorage = zeros(n,2);

BBstat0=zeros(n*B,2);
T_N0=zeros(1,2);
pvalue_all0=ones(1,2);
cvalues = zeros(3,2);

for j = 1:n % Loop through number of times to run calculations
    
    
    % Storage
  %  pvalues = zeros(1,2);
   
    for i = 1:2
        
        switch i
            case 1
                sample1 = s1;
                sample2 = s2;               
            case 2
                sample1 = s2;
                sample2 = s1;            
        end
        %------------------------ Subroutine ------------------------%
        % Measure the size of data
        N1=size(sample1,1);
        N2=size(sample2,1);
        N=N1+N2;
        % Construct a Support
        pooled=sort([sample1; sample2]);
        grid=linspace(min(pooled),max(pooled),ngrid)';
        % Setting up the tuning parameter
        seq_c=c*log(log(N))/sqrt(N);
        %seq_c=c*log(log(200))/sqrt(200);
        % Compute the Test Statistic
        operator=@(X,z)(X<= z).*(z-X).^(SDorder-1)/factorial(SDorder-1);
        rawcdf1=mean(bsxfun(operator, sample1, reshape(grid,[1,1,ngrid])),1);
        rawcdf2=mean(bsxfun(operator, sample2, reshape(grid,[1,1,ngrid])),1);
        cdf1=squeeze(rawcdf1);
        cdf2=squeeze(rawcdf2);
        base=(max(pooled)-min(pooled))/ngrid;
        cmstat=N*trapz( ((cdf1-cdf2)>0).*(cdf1-cdf2).^ 2)*base;
        if (j==1)
            T_N0(1,i)=cmstat;
        end
        % Bootstrap
        index1=randi(N1,N1,B);
        index2=randi(N2,N2,B);
        bsample1=sample1(index1); % bootstrap sample
        bsample2=sample2(index2); % bootstrap sample
        bcdf1=mean(bsxfun(operator, bsample1, reshape(grid, [1, 1,ngrid])), 1) - repmat(rawcdf1, [1, B, 1]);
        bcdf2=mean(bsxfun(operator, bsample2, reshape(grid, [1, 1,ngrid])), 1) - repmat(rawcdf2, [1, B, 1]);
        
        contact=(abs(cdf1-cdf2)<seq_c); % contact set, i.e. how many points between cdf1 and cdf2 have a very smol difference, less than seq_c
        if sum(contact)==0
            bcmstat=N*trapz( ((bcdf1 - bcdf2)>0).*(bcdf1 - bcdf2).^ 2, 3)*base;
        else
            ct_rep=repmat(contact(contact==1),[1,B]);
            ct_num=size(ct_rep,1);
            bcdf1=transpose(squeeze(bcdf1));
            bcdf2=transpose(squeeze(bcdf2));
            bcdf1=reshape(bcdf1(ct_rep==1),[ct_num, B]);
            bcdf2=reshape(bcdf2(ct_rep==1),[ct_num, B]);
            bcmstat=N*trapz( ((bcdf1-bcdf2)>0).*(bcdf1-bcdf2).^ 2)*base;
        end
        pvalue= sum(bcmstat>cmstat)/B;
   %     pvalues(i) = pvalue;
        clearvars sample1 sample2
        tmpn1=(j-1)*B+1;
        tmpn2=j*B;
        BBstat0([tmpn1:tmpn2],i)=bcmstat;
    end
    
   % allstorage(j,:) = pvalues;
    
end
%%%%%%%%%%%%%%%%%%%%%%



for kk=1:2
    pvalue_all0(1,kk)= sum(BBstat0(1:5001,kk)>T_N0(1,kk))/5001; 
    cvalues(:,kk) = prctile(BBstat0(1:5001,kk),[90,95,99]);
end



%-------------------------- Output --------------------------%
% Get the number of values that are less than some alpha value

%lessthan001 = [sum(allstorage<0.01,1)]./n;
%lessthan005 = [sum(allstorage<0.05,1)]./n;
%lessthan01 = [sum(allstorage<0.1,1)]./n;

%lessthanall = [lessthan001;lessthan005;lessthan01];

%t = table([lessthan001;lessthan005;lessthan01]);

csvwrite('pvaluesS4_S5.txt',round(pvalue_all0,4));
csvwrite('cvaluesS4_S5.txt',round(cvalues,4));
csvwrite('S4_S5_N0bnew1.txt',round(T_N0,4));
csvwrite('S4_S5BBstat0b.txt',BBstat0);


toc;

