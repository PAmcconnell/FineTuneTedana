function J = jaccard(A,B)
    J = sum(A.*B)/(sum(A+B) - sum(A.*B));
end
