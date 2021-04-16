function [TD, TG] = downsample(M,scalar,rate,mode)
M_size = size(M);

%max_num = max(max(max(M)));
%=======vertical and horizontal sampling==================================%
if mode == 'RGB'
    if (mod(M_size(1), rate) ~= 0) || (mod(M_size(2), rate) ~= 0)
        error('The size of input tensor cannot be exactly divided by rate')
    end
    TG = scalar*ones(M_size);
    TD = M(1:rate:M_size(1), 1:rate:M_size(2), :);
    TG(1:rate:M_size(1), 1:rate:M_size(2), :) = TD;
else
    if (mod(M_size(1), rate) ~= 0) || (mod(M_size(2), rate) ~= 0) || (mod(M_size(3), rate) ~= 0)
        error('The size of input tensor cannot be exactly divided by rate')
    end
    TG = scalar*ones(M_size);
    TD = M(1:rate:M_size(1), 1:rate:M_size(2), 1:rate:M_size(3));
    TG(1:rate:M_size(1), 1:rate:M_size(2), 1:rate:M_size(3)) = TD;
end
%=========interpolation===================================================%
% T = M(1:2:M_size(1), 1:2:M_size(2), :);
% T = imresize(T, 2, 'bicubic');
%=======diagonal sampling=================================================%
% T = scalar*ones(M_size);
% for i=1:M_size(1);
%     for j=1:M_size(2);
%         if mod(i+j,2)==0
%             T(i,j,:) = M(i,j,:);
%         end
%     end
% end
%=============circle sampling=============================================%
% TG = scalar*ones(M_size);
% max_size = max(M_size);
% for i=1:2:max_size
%     distance = round((i-1)*sqrt(2));
%     for j=1:(distance + 1)
%         index_i = round(sqrt(distance^2 - (j-1)^2)) + 1;
%         if (index_i <= M_size(1)) && (j <= M_size(2))
%             TG(index_i, j, :) = M(index_i, j, :);
%         end
%     end
% end
end