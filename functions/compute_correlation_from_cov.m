function corrvals = compute_correlation_from_cov(C1_all, C2_all, Filters)

nSegment = size(C1_all,3);
nComp = size(Filters,2);

for iComp = 1:nComp
    for iSegment = 1:nSegment
        val(iSegment) = Filters(:,iComp)'*C2_all(:,:,iSegment)*Filters(:,iComp);
    end
    cov1 = (Filters(:,iComp)'*mean(C1_all,3)*Filters(:,iComp));
    corrvals(iComp) = cov1/std(val);
end

end