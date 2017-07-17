function X = columnnorm(X)

nvec = size(X,2);
for iVec = 1:nvec
    vecn = norm(X(:,iVec));
    X(:,iVec) = X(:,iVec)/vecn;
end

end