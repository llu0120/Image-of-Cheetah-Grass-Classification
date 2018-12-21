clear all; 
%% Compute the priors Py(cheeta) & Py(grass) 
load('TrainingSamplesDCT_8_new.mat');
p_background = size(TrainsampleDCT_BG,1) / (size(TrainsampleDCT_BG,1) + size(TrainsampleDCT_FG,1)); % PY(grass)
p_foreground = 1 - p_background; %PY(cheeta) 

%% plot the 64 features in Guassian 
mean_background = mean(TrainsampleDCT_BG); 
sd_background = sqrt(var(TrainsampleDCT_BG));
mean_foreground = mean(TrainsampleDCT_FG);
sd_foreground = sqrt(var(TrainsampleDCT_FG));

figure(1);
for i = 1:64
    x_min = min(mean_background(1,i) - 3*sd_background(1,i), mean_foreground(1,i) - 3*sd_foreground(1,i));
    x_max = max(mean_background(1,i) + 3*sd_background(1,i), mean_foreground(1,i) + 3*sd_foreground(1,i));
    x = linspace(x_min,x_max);
    norm_background = normpdf(x, mean_background(1,i), sd_background(1,i));
    norm_foreground = normpdf(x, mean_foreground(1,i), sd_foreground(1,i));
    
    subplot(8,8,i);
    plot(x, norm_background,'B', x, norm_foreground,'R');
    title(i);  
    hold on;
end 
legend('background','foreground'); 



%% Best 8 features (1,8,18,24,34,40,41,44)
figure(2);
best = 1;
best_features = [1,8,18,24,34,40,41,44];
for i = best_features
    x_min = min(mean_background(1,i) - 3*sd_background(1,i), mean_foreground(1,i) - 3*sd_foreground(1,i));
    x_max = max(mean_background(1,i) + 3*sd_background(1,i), mean_foreground(1,i) + 3*sd_foreground(1,i));
    x = linspace(x_min,x_max);
    norm_background = normpdf(x, mean_background(1,i), sd_background(1,i));
    norm_foreground = normpdf(x, mean_foreground(1,i), sd_foreground(1,i));
    
    subplot(2,4,best);
    best = best + 1;
    plot(x, norm_background,'B', x, norm_foreground,'R');
    title(i);  
    hold on;
end 
legend('background','foreground'); 

%% Worst 8 features (2,3,4,5,6,59,63,64)
figure(3);
worst = 1;
worst_features = [2,3,4,5,6,59,63,64];
for i = worst_features
    x_min = min(mean_background(1,i) - 3*sd_background(1,i), mean_foreground(1,i) - 3*sd_foreground(1,i));
    x_max = max(mean_background(1,i) + 3*sd_background(1,i), mean_foreground(1,i) + 3*sd_foreground(1,i));
    x = linspace(x_min,x_max);
    norm_background = normpdf(x, mean_background(1,i), sd_background(1,i));
    norm_foreground = normpdf(x, mean_foreground(1,i), sd_foreground(1,i));
    
    subplot(2,4,worst);
    worst = worst + 1;
    plot(x, norm_background,'B', x, norm_foreground,'R');
    title(i);  
    hold on;
end 
legend('background','foreground'); 

%If Py|x(cheeta|x) > Py|x(grass|) ----> the pixel is cheeta
%same as Px|y(x|cheeta)Py(cheeta) > Px|y(x|grass)Py(grass) 
%input bmp file & ZigZag pattern
A = im2double(imread('cheetah.bmp'));
Z = load('Zig-Zag Pattern.txt');

%% original size 255*270 ---> need to use padarray() to fill 255+7, 270+7 
B = padarray(A,[7,7],'symmetric','post');
[q,l] = size(B);
matrix_zigzag=[];
mask = zeros(q,l);

%% make the picture into 8*8 pixels and use mask determine which pixel is cheeta/grass
var_background = var(TrainsampleDCT_BG);
var_foreground = var(TrainsampleDCT_FG);

for i=1:q-7
    for j=1:l-7
        matrix_dct2 = dct2(B(i:i+7,j:j+7));
        matrix_zigzag(Z+ones(8,8)) = matrix_dct2; %store matrix_dct into zigzag which is from 0~63 so plus a ones(8,8)
%         [V,I] = max(abs(matrix_zigzag(2:64)));
%%         compute class-conditioned Px|y(x|cheeta)
%         index_cheeta = zeros(1,64);
%         TrainsampleDCT_FG(:,1) = []; 
%         for ci = 1: size(TrainsampleDCT_FG,1)
%             [VC,C] = max(abs(TrainsampleDCT_FG(ci,:))); 
%             index_cheeta(C) = index_cheeta(C) + 1;
%         end
%         index_cheeta = index_cheeta / size(TrainsampleDCT_FG,1);

%%         compute class-conditioned Px|y(x|grass)
%         index_grass = zeros(1,64);
%         TrainsampleDCT_BG(:,1) = []; 
%         for gi = 1: size(TrainsampleDCT_BG,1)
%             [VG,G] = max(abs(TrainsampleDCT_BG(gi,:))); 
%             index_grass(G) = index_grass(G) + 1;
%         end
%         index_grass = index_grass / size(TrainsampleDCT_BG,1);

%% index_grass ----> Pxy_xgrass; index_cheeta -----> Pxy_xcheeta
   %Compute class-conditioned Px|y(x|cheeta) & compute class-conditioned Px|y(x|grass)
        Pxy_xgrass= mvnpdf(matrix_zigzag, mean_background, var_background);
        Pxy_xcheeta = mvnpdf(matrix_zigzag, mean_foreground, var_foreground);
        if (Pxy_xcheeta*p_foreground) >  (Pxy_xgrass*p_background)
            mask(i,j) = 1;
        else% Pxy_xcheeta*p_foreground < Pxy_xgrass*p_background 
            mask(i,j) = 0;
        end
    end
end
figure(4);
imagesc(mask);
colormap(gray(255));

%% Error Rate of all 64 features
diff_grass = 0;
diff_cheeta = 0; 
count_grass = 0;
count_cheeta = 0; 

groundtruth_mask = im2double(imread('cheetah_mask.bmp'));
for i=1:(q-7)
    for j=1:(l-7)
        if mask(i,j) == 1  
            if groundtruth_mask(i,j) == 0 
                diff_cheeta = diff_cheeta + 1; 
            else 
                diff_cheeta = diff_cheeta;
            end
        elseif mask(i,j) == 0
            if groundtruth_mask(i,j) == 1 
                diff_grass = diff_grass + 1; 
            else 
                diff_grass = diff_grass; 
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
error_cheeta64 = (diff_cheeta / count_cheeta)*p_foreground ; 
error_grass64 = (diff_grass / count_grass)*p_background;
error_rate64 = error_cheeta64 + error_grass64;


%% Best Features git practice
for i=1:q-7
    for j=1:l-7
        matrix_dct2 = dct2(B(i:i+7,j:j+7));
        matrix_zigzag(Z+ones(8,8)) = matrix_dct2; %store matrix_dct into zigzag which is from 0~63 so plus a ones(8,8)
%% index_grass ----> Pxy_xgrass; index_cheeta -----> Pxy_xcheeta
    %Compute class-conditioned Px|y(x|cheeta) & compute class-conditioned Px|y(x|grass)
        Pxy_xgrass = mvnpdf(matrix_zigzag(best_features), mean_background(best_features), var_background(best_features)); 
        Pxy_xcheeta = mvnpdf(matrix_zigzag(best_features), mean_foreground(best_features), var_foreground(best_features));
        if (Pxy_xcheeta*p_foreground) >  (Pxy_xgrass*p_background)
            mask(i,j) = 1;
        else% Pxy_xcheeta*p_foreground < Pxy_xgrass*p_background 
            mask(i,j) = 0;
        end
    end
end
figure(5);
imagesc(mask);
colormap(gray(255));

%% Error Rate of 8 best features
diff_grass = 0;
diff_cheeta = 0; 
count_grass = 0;
count_cheeta = 0; 

groundtruth_mask = im2double(imread('cheetah_mask.bmp'));
for i=1:(q-7)
    for j=1:(l-7)
        if mask(i,j) == 1  
            if groundtruth_mask(i,j) == 0 
                diff_cheeta = diff_cheeta + 1; 
            else 
                diff_cheeta = diff_cheeta;
            end
        elseif mask(i,j) == 0
            if groundtruth_mask(i,j) == 1 
                diff_grass = diff_grass + 1; 
            else 
                diff_grass - diff_grass; 
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
error_cheetab = (diff_cheeta / count_cheeta)*p_foreground ; 
error_grassb = (diff_grass / count_grass)*p_background;
error_rateb = error_cheetab + error_grassb;

%% Worst Features
for i=1:q-7
    for j=1:l-7
        matrix_dct2 = dct2(B(i:i+7,j:j+7));
        matrix_zigzag(Z+ones(8,8)) = matrix_dct2; %store matrix_dct into zigzag which is from 0~63 so plus a ones(8,8)
%% index_grass ----> Pxy_xgrass; index_cheeta -----> Pxy_xcheeta
    %Compute class-conditioned Px|y(x|cheeta) & compute class-conditioned Px|y(x|grass)
        Pxy_xgrass = mvnpdf(matrix_zigzag(worst_features), mean_background(worst_features), var_background(worst_features)); 
        Pxy_xcheeta = mvnpdf(matrix_zigzag(worst_features), mean_foreground(worst_features), var_foreground(worst_features));
        if (Pxy_xcheeta*p_foreground) >  (Pxy_xgrass*p_background)
            mask(i,j) = 1;
        else% Pxy_xcheeta*p_foreground < Pxy_xgrass*p_background 
            mask(i,j) = 0;
        end
    end
end
figure(6);
imagesc(mask);
colormap(gray(255));

%% Error Rate of 8 worst features
diff_grass = 0;
diff_cheeta = 0; 
count_grass = 0;
count_cheeta = 0; 

groundtruth_mask = im2double(imread('cheetah_mask.bmp'));
for i=1:(q-7)
    for j=1:(l-7)
        if mask(i,j) == 1  
            if groundtruth_mask(i,j) == 0 
                diff_cheeta = diff_cheeta + 1; 
            else 
                diff_cheeta = diff_cheeta;
            end
        elseif mask(i,j) == 0
            if groundtruth_mask(i,j) == 1 
                diff_grass = diff_grass + 1; 
            else 
                diff_grass - diff_grass; 
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
error_cheetaw = (diff_cheeta / count_cheeta)*p_foreground ; 
error_grassw = (diff_grass / count_grass)*p_background;
error_ratew = error_cheetaw + error_grassw;
