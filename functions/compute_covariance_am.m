function Cw = compute_covariance_am(C_all, y)

yn = (y-mean(y))/std(y);
for iSegment = 1:size(C_all,3)
    Cw(:,:,iSegment) = (C_all(:,:,iSegment))*yn(iSegment);
end

end