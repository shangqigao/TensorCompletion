function filterX = mygaussfilter(X, mode)
%This function is used to calculate the gaussian filter of X with sigma=1
%input:
%   X: The 3D tensor need to filter
%   mode: color image or not, if color image, mode='RGB', else, mode='NOT'
%output:
%   filterX: The tensor filtered by gaussian filter
dim = size(X);
filterX = zeros(dim);
if mode == 'RGB'
    H = fspecial('gaussian', 3, 1);
    for i=1:dim(3)
        filterX(:,:,i) = imfilter(X(:,:,i), H, 'conv', 'circular');
    end
else
    filterX = gauss3filter(X, 1);
end
    
end