function[rank] = rankestimation(M)
%M=double(imread('real_and_fake_peppers_RGB.bmp'));
dim=size(M);
normM=norm(M(:));
ndim=length(dim);
X=cell(ndim,1);
for i=1:ndim;
    X{i}=Unfold(M,dim,i);
    d=sqrt(svd(X{i}*X{i}'));
    rank(i)=length(d(d>1e-2*max(d)));
%     rank(i)=length(d(d>5e-2*normM));
end
end