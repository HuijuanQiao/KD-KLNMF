clc
clear 
S1=fastaread('\Nucleolus.txt');
S2=fastaread('\Cytoplasm.txt');
S3=fastaread('\Ribosome.txt');
S4=fastaread('\Exosome.txt');
tt=[S1;S2;S3;S4];
G=10;
n=length(tt);
feature2=zeros(n,G*225);
for i=1:n
    SSS{1,i}=tt(i).Sequence;
    Seq=SSS{1,i};    
    %生成DNA样本的特征向量 
    value_na = cell(1,16);
    value_na{1} = 'AA'; value_na{2} = 'TT'; value_na{3} = 'AT'; value_na{4} = 'AG';
    value_na{5} = 'CT'; value_na{6} = 'AC'; value_na{7} = 'GT'; value_na{8} = 'TA';
    value_na{9} = 'TG'; value_na{10} = 'CA'; value_na{11} = 'TC'; value_na{12} = 'GA';
    value_na{13} = 'GG'; value_na{14} = 'CC'; value_na{15} = 'GC'; value_na{16} = 'CG';
    
    Pro = [0.041 0.041 0.054 0.042 0.042 0.065 0.065 0.031 0.035 0.035 0.049 0.049 0.041 0.041 0.054 0.039;
     0.078 0.078 0.098 0.058 0.058 0.071 0.071 0.065 0.056 0.056 0.065 0.065 0.057 0.057 0.065 0.059;
     0.069 0.069 0.071 0.053 0.053 0.064 0.064 0.048 0.052 0.052 0.057 0.057 0.055 0.055 0.059 0.051;
     6.689 6.689 9.611 3.472 3.472 6.803 6.803 1.853 2.003 2.003 4.268 4.268 2.99 2.99 4.206 2.713;
     6.239 6.239 4.658 2.801 2.801 2.911 2.911 4.107 2.882 2.882 3.58 3.58 2.67 2.67 2.655 3.019;
    21.34 21.34 24.792 17.477 17.477 21.977 21.977 14.239 14.512 14.512 18.41 18.41 14.252 14.252 17.311 14.655;
    1.051 1.051 0.612 3.6 3.6 2.005 2.005 3.499 5.6 5.6 2.444 2.444 4.682 4.682 1.697 6.015;
    -1.261 -1.261 0 -1.655 -1.655 0.334 0.334 0 0.137 0.137 1.437 1.437 -0.77 -0.77 0 0;
    35.02 35.02 30.72 32.29 32.29 31.53 31.53 36.94 35.43 35.43 35.67 35.67 33.54 33.54 34.07 33.67;
    -0.176 -0.176 -0.679 -0.223 -0.223 -0.593 -0.593 0.044 0.481 0.481 -0.046 -0.046 -0.166 -0.166 -0.19 0.443;
    0.013 0.013 0 -0.023 -0.023 -0.018 -0.018 0 0.009 0.009 -0.011 -0.011 0.026 0.026 0 0;
    3.253 3.253 3.208 3.322 3.322 3.243 3.243 3.389 3.366 3.366 3.299 3.299 3.361 3.361 3.267 3.291;
    -1 -1 -0.88  -1.28 -1.28 -1.44 -1.44 -0.58 -1.45 -1.45 -1.3 -1.3 -1.84 -1.84 -2.24 -2.17;
    -7.6 -7.6 -7.2 -7.8 -7.8 -8.4 -8.4 -7.2 -8.5 -8.5 -8.2 -8.2 -8 -8 -9.8 -10.6;
    -21.3 -21.3 -20.4 -21 -21 -22.4 -22.4 -21.3 -22.7 -22.7 -22.2 -22.2 -19.9 -19.9 -24.4 -27.2];
    for ii = 1:size(Pro,1)
        Pro(ii,:) = (Pro(ii,:)-min(Pro(ii,:)))/(max(Pro(ii,:))-min(Pro(ii,:)));
    end
    
    value_matrix = zeros(size(Pro,1),length(Seq)-1);
    for ii = 1:size(Pro,1)
        for jj = 1:16
            value_matrix(ii,strfind(Seq,value_na{jj})) = Pro(ii,jj);
        end
    end
    value_matrix = value_matrix';
    M = zeros(15*15,G);
    FF = zeros(15);
    for g = 1:G
        for s = 1:15
            for t = 1:15
            Csg = value_matrix(1:length(Seq)-1-g,s);
            Csg1= value_matrix(1+g:length(Seq)-1,s) ;
            Ctg= value_matrix(1:length(Seq)-1-g,t) ;
            Ctg1 = value_matrix((1+g):(length(Seq)-1),t);
            C1=value_matrix(:,s);
            C2=value_matrix(:,t);       
            t1=sum((Csg-Csg1).*(Ctg-Ctg1));
            l1=2*(length(Seq)-1-g);
            t2=sum(C1-mean(C1).*(C2-mean(C2)));
            l2=(length(Seq)-1-1);            
            FF(s,t) =(t1/l1)/(t2/l2);

            end
        end
        M(:,g) = FF(:);
    end
    feature2(i,:) = M(:)';
end
%csvwrite('\DSA2250.csv',feature2);