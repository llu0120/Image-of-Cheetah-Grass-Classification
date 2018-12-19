
% load('sigmaFG_cell.mat');
% load('sigmaBG_cell.mat');
% load('priorFG_cell.mat');
% load('priorBG_cell.mat');
% load('muFG_cell.mat');
% load('muBG_cell.mat');

% %% Training Part
% load('./TrainingSamplesDCT_8_new.mat');
% [sample_FG, dim] = size(TrainsampleDCT_FG);
% [sample_BG, dim] = size(TrainsampleDCT_BG);
%
% priorFG_cell = cell(6, 1);
% priorBG_cell = cell(6, 1);
% muFG_cell = cell(6, 1);
% muBG_cell = cell(6, 1);
% sigmaFG_cell = cell(6, 1);
% sigmaBG_cell = cell(6, 1);
% class = [1,2,4,8,16,32];
% Dim = 64;
% Iter = 30;
% countc=1;
% for c = class
% %% Compute 5 sets of FG & BG parameters through EM algorithm (30 iterations)
% Mixture = 1;
% for M = 1
%     %E-step initialize
%     priorFG = ones(1, c)./c;
%     priorBG = ones(1, c)./c;
%     muFG = rand(c, Dim);
%     muBG = rand(c, Dim);
%     diagFG = cell(c, 1);
%     diagBG = cell(c, 1);
%     for j = 1: c
%         temp = rand(1, Dim);
%         temp(temp<0.001) = 0.001;
%         diagFG{j} = diag(temp);
%         temp = rand(1, Dim);
%         temp(temp<0.001) = 0.001;
%         diagBG{j} = diag(temp);
%     end
%
%     %E-step
%     h_FG = zeros(sample_FG, c);
%     h_BG = zeros(sample_BG, c);
%     muFGprev = muFG;
%     muBGprev = muBG;
%     priorFGprev = priorFG;
%     priorBGprev = priorBG;
%     diagFGprev = diagFG;
%     diagBGprev = diagBG;
%     for i = 1: Iter
%         i
%         for j = 1: sample_FG
%             for z = 1: c
%             %compute hij
%                 pdf = mvnpdf(TrainsampleDCT_FG(j, 1:Dim), muFGprev(z, :), diagFGprev{z});
%                 h_FG(j, z) = pdf*priorFGprev(1, z);
%             end
%             h_FG(j, :) = h_FG(j, :)/sum(h_FG(j, :));
%         end
%
%         for z = 1: c
%            muFG(z, :) = h_FG(:, z)'*TrainsampleDCT_FG(:, 1:Dim)./sum(h_FG(:, z));
%            priorFG(1, z) = sum(h_FG(:, z))/sample_FG;
%         end
%
%         for z = 1: c
%             sig_tmp = zeros(1, dim);
%             sum_1 = 0;
%             sum_2 = 0;
%             for j = 1: sample_FG
%                 sum_1 = sum_1 + h_FG(j, z)* (TrainsampleDCT_FG(j, 1:Dim)-muFG(z, :)).^2;
%                 sum_2 = sum_2 + h_FG(j, z);
%
%             end
%             diag_tmp = sum_1/sum_2;
%             diag_tmp(diag_tmp<0.002) = 0.002;
%             diagFG{z} = diag(diag_tmp);
%         end
%         if max(max(abs(muFG - muFGprev)./abs(muFGprev))) < 0.01
%            break;
%         end
%
%         muFGprev = muFG;
%         priorFGprev = priorFG;
%         diagFGprev = diagFG;
%     end
%
%     for i = 1: Iter
%             i
%             for j = 1: sample_BG
%                 for z = 1: c
%                 %compute hij
%                     pdf = mvnpdf(TrainsampleDCT_BG(j, 1:Dim), muBGprev(z, :), diagBGprev{z});
%                     h_BG(j, z) = pdf*priorBGprev(1, z);
%                 end
%                 h_BG(j, :) = h_BG(j, :)/sum(h_BG(j, :));
%             end
%             %BG
%
%             for z = 1: c
%                muBG(z, :) = h_BG(:, z)'*TrainsampleDCT_BG(:, 1:Dim)./sum(h_BG(:, z));
%                priorBG(1, z) = sum(h_BG(:, z))/sample_BG;
%             end
%
%             for z = 1: c
%                 sum_1 = 0;
%                 sum_2 = 0;
%                 for j = 1: sample_BG
%                     sum_1 = sum_1 + h_BG(j, z)* (TrainsampleDCT_BG(j, 1:Dim)-muBG(z, :)).^2;
%                     sum_2 = sum_2 + h_BG(j, z);
%                 end
%                 diag_tmp = sum_1/sum_2;
%                 diag_tmp(diag_tmp<0.002) = 0.002;
%                 diagBG{z} = diag(diag_tmp);
%             end
%             if max(max(abs(muBG - muBGprev)./abs(muBGprev))) < 0.01
%                break;
%             end
%
%             muBGprev = muBG;
%             priorBGprev = priorBG;
%             diagBGprev = diagBG;
%
%     end
%     priorFG_cell{countc} = priorFG;
%     priorBG_cell{countc} = priorBG;
%     muFG_cell{countc} = muFG;
%     muBG_cell{countc} = muBG;
%     sigmaFG_cell{countc} = diagFG;
%     sigmaBG_cell{countc} = diagBG;
%     countc =countc + 1;
% end
% end
% %% prior probability
%  P_BG = sample_BG/(sample_BG + sample_FG);
%  P_FG = 1-P_BG;
%% 
%input bmp file & ZigZag pattern
%store matrix_dct into zigzag which is from 0~63 so plus one to each number
class=8;
 A = im2double(imread('cheetah.bmp'));
 Z = load('Zig-Zag-Pattern.txt');
 Z = Z+1;
 B = padarray(A,[7,7],'symmetric','post');% original size 255*270 ---> need to use padarray() to fill 255+7, 270+7 
 [q,l] = size(B);
 matrix_zigzag=[];
 mask = zeros(q,l);
 dimension = [1,2,4,8,16,24,32,40,48,56,64];
 feature_map = zeros(255*270, 64);
% 
%      
% 
% % sliding window
count = 1;
 for i=1:q-7
            for j=1:l-7
                
            matrix_dct2 = dct2(B(i:i+7,j:j+7));
            matrix_zigzag(Z)= matrix_dct2;
            matrix_zigzag=reshape(matrix_zigzag,[1,64]);
            feature_map(count,:) = matrix_zigzag;
            count=count+1;
            end
end
% % count # of FG & BG from the original mask 
count_grass = 0;
count_cheeta = 0;
groundtruth_mask = im2double(imread('cheetah_mask.bmp'));
for o=1:(q-7)
            for p=1:(l-7)
                if groundtruth_mask(o,p) == 1 
                    count_cheeta = count_cheeta + 1; 
                else 
                    count_grass = count_grass + 1;
                end     
            end
end
% 
% % 5 FG & BG mixture model
% % plug in solutions to mvnpdf ---> BDR
error_vector=[];
for n1 = 1:5
    n1
    for n2 = 1:5
        n2
        poe_error=[];
        for dim = dimension
            dim
            count = 1;
            PI_BG = pi_BG{n1};
            MU_BG = mu_BG{n1};
            SIGMA_BG = sigma_BG{n1};
            PI_FG = pi_FG{n2};
            MU_FG = mu_FG{n2};
            SIGMA_FG = sigma_FG{n2};
        for i=1:255
             for j=1:270
            Pxy_xgrass = 0;
            Pxy_xcheeta = 0;
            
            for c = 1:class
             
                
                G_grass = mvnpdf((feature_map(count,1:dim)), MU_BG(c,1:dim), diag(SIGMA_BG(c,1:dim)));
            Pxy_xgrass = Pxy_xgrass + PI_BG(1,c)*G_grass;
        
                G_cheeta = mvnpdf((feature_map(count,1:dim)), MU_FG(c,1:dim), diag(SIGMA_FG(c,1:dim)));
            Pxy_xcheeta= Pxy_xcheeta + PI_FG(1,c)*G_cheeta;
            end
       
        if  (Pxy_xcheeta*P_FG) >  (Pxy_xgrass*P_BG)
            mask(i,j) = 1;
        else% Pxy_xcheeta*p_foreground < Pxy_xgrass*p_background 
            mask(i,j) = 0;
        end
        count=count+1;
        end
    end
        %Error Rate
        error = [];
        
        diff_grass = 0;
        diff_cheeta = 0; 
         
        for o=1:(q-7)
            for p=1:(l-7)
                if mask(o,p) == 0  
                    if groundtruth_mask(o,p) == 1
                        diff_cheeta = diff_cheeta + 1; 
                    end
                elseif mask(o,p) == 1
                    if groundtruth_mask(o,p) == 0 
                        diff_grass = diff_grass + 1; 
                    end
                end   
            end
        end 

        
        error_cheeta = (diff_cheeta / count_cheeta)*P_FG ; 
        error_grass = (diff_grass / count_grass)*P_BG;
        error = error_cheeta + error_grass;
        
        poe_error = [poe_error error];
        end
        error_vector = [error_vector; poe_error];

        
    end
end

%% plot
[row,col]=size(error_vector);
figure;
for r = 1:row
    plot(dimension,error_vector(r,:),'-o');
    hold on;
end
title('POE and Dimemsion with different parameters')
xlabel('Dimensions');
ylabel('probability of error');
legend;
hold off; 
