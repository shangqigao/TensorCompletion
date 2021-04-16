function U = upsample(TD, rate, mode)
D_size = size(TD);
if mode == 'RGB'
    U = zeros(rate*D_size(1), rate*D_size(2), D_size(3));
    for i=1:rate
        for j=1:rate
            U(i:rate:end, j:rate:end, :) = TD;
        end
    end
else
    U = zeros(D_size*rate);
    for i=1:rate
        for j=1:rate
            for k=1:rate
                U(i:rate:end, j:rate:end, k:rate:end) = TD;
            end
        end
    end
end
end