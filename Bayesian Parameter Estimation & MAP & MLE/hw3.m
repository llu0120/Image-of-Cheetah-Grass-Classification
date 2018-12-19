clear all; 
load('TrainingSamplesDCT_subsets_8.mat');
alpha = load('Alpha.mat');
prior_1 = load('Prior_1.mat');
prior_2 = load('Prior_2.mat');

%%
Data_BG = D4_BG;
Data_FG = D4_FG;
[m,n] = size(Data_BG);
[p,s] = size(Data_FG);
[a,b] = size(alpha.alpha);
mu_BG = mean(Data_BG);
mu_FG = mean(Data_FG);
p_BG = size(Data_BG,1)/(size(Data_BG,1)+size(Data_FG,1));
p_FG = 1-p_BG;

%Compute covariance BG & FG 
cov_BG = zeros(n,n);
cov_FG = zeros(n,n);
for i = 1:m
    cov_BG = cov_BG + (Data_BG(i,:) - mu_BG).' * (Data_BG(i,:) - mu_BG);
end
    cov_BG = cov_BG/m;
for i = 1:p
    cov_FG = cov_FG + (Data_FG(i,:) - mu_FG).' * (Data_FG(i,:) - mu_FG);
end
    cov_FG = cov_FG/p;

cov_0 = zeros(n,n);
for k = 1:b %9
    %Compute covariance0 
    for h = 1:n %64
        cov_0(h,h) = alpha.alpha(k) * prior_2.W0(h);
    end
    
    %Compute mu1,covariance1 of BG & FG 
    mu_1_BG = ((cov_0/(cov_0+(1/m)*cov_BG))*mu_BG.'+(1/m)*(cov_BG/(cov_0+(1/m)*cov_BG))*prior_2.mu0_BG.');          
    cov_1_BG = (cov_0/(cov_0+(1/m)*cov_BG))*(1/m)*cov_BG;
    mu_1_FG = ((cov_0/(cov_0+(1/p)*cov_FG))*mu_FG.'+(1/p)*(cov_FG/(cov_0+(1/p)*cov_FG))*prior_2.mu0_FG.');
    cov_1_FG = (cov_0/(cov_0+(1/p)*cov_FG))*(1/p)*cov_FG;
    
%% 
%input bmp file & ZigZag pattern
A = im2double(imread('cheetah.bmp'));
Z = load('Zig-Zag Pattern.txt');
B = padarray(A,[7,7],'symmetric','post');% original size 255*270 ---> need to use padarray() to fill 255+7, 270+7 
[q,l] = size(B);
matrix_zigzag=[];
BAYES = zeros(q,l);
MAP = zeros(q,l);
ML = zeros(q,l);
for i=1:q-7
    for j=1:l-7
        matrix_dct2 = dct2(B(i:i+7,j:j+7));
        matrix_zigzag(Z+ones(8,8)) = matrix_dct2; %store matrix_dct into zigzag which is from 0~63 so plus a ones(8,8)
        
%% plug in three solutions to mvnpdf ---> BDR
        cov_sum_BG = cov_1_BG + cov_BG;
        cov_sum_FG = cov_1_FG + cov_FG;
        Pxy_xgrass_BAYES= mvnpdf((matrix_zigzag.'), mu_1_BG, cov_sum_BG);
        Pxy_xcheeta_BAYES= mvnpdf((matrix_zigzag.'), mu_1_FG, cov_sum_FG);
        if  (Pxy_xcheeta_BAYES*p_FG) >  (Pxy_xgrass_BAYES*p_BG)
           BAYES(i,j) = 1;
        else% Pxy_xcheeta*p_foreground < Pxy_xgrass*p_background 
           BAYES(i,j) = 0;
        end
        
        
        Pxy_xgrass_MAP= mvnpdf((matrix_zigzag.'), mu_1_BG, cov_BG);
        Pxy_xcheeta_MAP= mvnpdf((matrix_zigzag.'), mu_1_FG, cov_FG);
        if  (Pxy_xcheeta_MAP*p_FG) >  (Pxy_xgrass_MAP*p_BG)
           MAP(i,j) = 1;
        else% Pxy_xcheeta*p_foreground < Pxy_xgrass*p_background 
           MAP(i,j) = 0;
        end
        
        
        Pxy_xgrass_ML= mvnpdf(matrix_zigzag, mu_BG, cov_BG);
        Pxy_xcheeta_ML= mvnpdf(matrix_zigzag, mu_FG, cov_FG);
        if  (Pxy_xcheeta_ML*p_FG) >  (Pxy_xgrass_ML*p_BG)
            ML(i,j) = 1;
        else% Pxy_xcheeta*p_foreground < Pxy_xgrass*p_background 
            ML(i,j) = 0;
        end   
    end
end
%% Error Rate of BAYES
groundtruth_mask = im2double(imread('cheetah_mask.bmp'));
diff_grass = 0;
diff_cheeta = 0; 
count_grass = 0;
count_cheeta = 0; 
for i=1:(q-7)
    for j=1:(l-7)
        if BAYES(i,j) == 1  
            if groundtruth_mask(i,j) == 0 
                diff_cheeta = diff_cheeta + 1; 
            end
        elseif BAYES(i,j) == 0
            if groundtruth_mask(i,j) == 1 
                diff_grass = diff_grass + 1; 
            end
        end   
     end
end 

for i=1:(q-7)
    for j=1:(l-7)
        if groundtruth_mask(i,j) == 1 
            count_cheeta = count_cheeta + 1; 
        else 
            count_grass = count_grass + 1;
        end     
    end
end 
error_cheeta = (diff_cheeta / count_cheeta)*p_FG ; 
error_grass = (diff_grass / count_grass)*p_BG;
error_BAYES(k) = error_cheeta + error_grass;

%% Error Rate of MAP
%groundtruth_mask = im2double(imread('cheetah_mask.bmp'));
diff_grass = 0;
diff_cheeta = 0; 
count_grass = 0;
count_cheeta = 0; 
for i=1:(q-7)
    for j=1:(l-7)
        if MAP(i,j) == 1  
            if groundtruth_mask(i,j) == 0 
                diff_cheeta = diff_cheeta + 1; 
            end
        elseif MAP(i,j) == 0
            if groundtruth_mask(i,j) == 1 
                diff_grass = diff_grass + 1; 
            end
        end   
     end
end 

for i=1:(q-7)
    for j=1:(l-7)
        if groundtruth_mask(i,j) == 1 
            count_cheeta = count_cheeta + 1; 
        else 
            count_grass = count_grass + 1;
        end     
    end
end 
error_cheeta = (diff_cheeta / count_cheeta)*p_FG ; 
error_grass = (diff_grass / count_grass)*p_BG;
error_MAP(k) = error_cheeta + error_grass;

%% Error Rate of ML
%groundtruth_mask = im2double(imread('cheetah_mask.bmp'));
diff_grass = 0;
diff_cheeta = 0; 
count_grass = 0;
count_cheeta = 0; 
for i=1:(q-7)
    for j=1:(l-7)
        if ML(i,j) == 1  
            if groundtruth_mask(i,j) == 0 
                diff_cheeta = diff_cheeta + 1; 
            end
        elseif ML(i,j) == 0
            if groundtruth_mask(i,j) == 1 
                diff_grass = diff_grass + 1; 
            end
        end   
     end
end 

for i=1:(q-7)
    for j=1:(l-7)
        if groundtruth_mask(i,j) == 1 
            count_cheeta = count_cheeta + 1; 
        else 
            count_grass = count_grass + 1;
        end     
    end
end 
error_cheeta = (diff_cheeta / count_cheeta)*p_FG ; 
error_grass = (diff_grass / count_grass)*p_BG;
error_ML(k) = error_cheeta + error_grass;

end


x = alpha.alpha;
yB = error_BAYES;
plot(x,yB,'r');
set(gca, 'XScale', 'log');
hold on;

yMAP = error_MAP;
plot(x,yMAP,'b');
set(gca, 'XScale', 'log');
hold on;

yML = error_ML;
plot(x,yML,'g');
set(gca, 'XScale', 'log');
hold off;
xlabel('a');
ylabel('probability of error');
legend('Bayes','MAP','ML');