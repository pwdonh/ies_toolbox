function corrvals = compute_correlation_from_cov(C1_all, C2_all, Filters)

nSegment = size(C2_all,3);
nComp = size(Filters,2);

fprintf('Computing correlation source %5d of %5d\n', 0, nComp);
for iComp = 1:nComp
    fprintf([repmat('\b',1,15) '%5d of %5d\n'], iComp, nComp);
    for iSegment = 1:nSegment
        val(iSegment) = Filters(:,iComp)'*C2_all(:,:,iSegment)*Filters(:,iComp);
    end
    cov1 = (Filters(:,iComp)'*mean(C1_all,3)*Filters(:,iComp));
    corrvals(iComp) = cov1/std(val);
end

end