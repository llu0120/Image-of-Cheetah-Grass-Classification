clear all; 
%% TrainingSamplesDCT_8.mat
%compute the priors Py(cheeta) & Py(grass) 
load('TrainingSamplesDCT_8.mat');
p_background = size(TrainsampleDCT_BG,1) / (size(TrainsampleDCT_BG,1) + size(TrainsampleDCT_FG,1)); % PY(grass)
p_foreground = 1 - p_background; %PY(cheeta) 

%compute class-conditioned Px|y(x|cheeta)
index_cheeta = zeros(1,64);
TrainsampleDCT_FG(:,1) = []; 
for ci = 1: size(TrainsampleDCT_FG,1)
    [VC,C] = max(abs(TrainsampleDCT_FG(ci,:))); 
    index_cheeta(C) = index_cheeta(C) + 1;
end
index_cheeta = index_cheeta / size(TrainsampleDCT_FG,1);

%compute class-conditioned Px|y(x|grass)
index_grass = zeros(1,64);
TrainsampleDCT_BG(:,1) = []; 
for gi = 1: size(TrainsampleDCT_BG,1)
    [VG,G] = max(abs(TrainsampleDCT_BG(gi,:))); 
    index_grass(G) = index_grass(G) + 1;
end
index_grass = index_grass / size(TrainsampleDCT_BG,1);

%plot on histogram
figure(1);
subplot(1,2,1);
bar(index_cheeta);
title('Px|y(x|cheeta)');

subplot(1,2,2);
bar(index_grass);
title('Px|y(x|grass)');

%If Py|x(cheeta|x) > Py|x(grass|x) ----> the pixel is cheeta
%same as Px|y(x|cheeta)Py(cheeta) > Px|y(x|grass)Py(grass) 
%input bmp file & ZigZag pattern
A = imread('cheetah.bmp');
Z = load('Zig-Zag Pattern.txt');
%% original size 255*270 ---> need to use padarray() to fill 255+7, 270+7 
B = padarray(A,[7,7],'symmetric','post');
[q,l] = size(B);
matrix_zigzag=[];
mask = zeros(q,l);
%% make the picture into 8*8 pixels and use mask determine which pixel is cheeta/grass
for i=1:q-7
    for j=1:l-7
        matrix_dct2 = dct2(B(i:i+7,j:j+7));
        matrix_zigzag(Z+ones(8,8)) = matrix_dct2; %store matrix_dct into zigzag which is from 0~63 so plus a ones(8,8)
        [V,I] = max(abs(matrix_zigzag(2:64)));
        if index_cheeta(I)*p_foreground > index_grass(I)*p_background 
            mask(i,j) = 1;
        elseif index_cheeta(I)*p_foreground < index_grass(I)*p_background 
            mask(i,j) = 0;
        else           %if index_cheeta(I) == index_grass(I) == 0 ---> use outer 2*2  
            white = 0; 
            black = 0;
            for z=i-2:i+2
                for w=j-2:j+2 
                    if mask(z,w)== 1
                        white = white + 1;
                    elseif mask(z,w) == 0 
                        black = black + 1;
                    end
                end
            end
            if black > white
               mask(i,j) = 1;
            else
               mask(i,j) = 0;
            end
        end
    end
end
figure(2);
imagesc(mask);
colormap(gray(255));

%% compute error rate 
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
error_cheeta = (diff_cheeta / count_cheeta)*p_foreground ; 
error_grass = (diff_grass / count_grass)*p_background;
error_rate = error_cheeta + error_grass;






