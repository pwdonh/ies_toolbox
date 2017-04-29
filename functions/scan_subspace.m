function [sc, x_all] = scan_subspace(Us, G, nIter, iSeed, isOutOrient)

if nargin<5
    isOutOrient = 0;
end

nDip = size(G,3);
sc = zeros(nDip,1);
if ~isOutOrient
    xx = [];
    parfor iDip = 1:nDip
        sc(iDip) = subcorr(G(:,:,iDip), Us, 1);
    end
else
    x_all = zeros(3,nDip,nIter);
    parfor iDip = 1:nDip
        [sc(iDip), x_all(:,iDip,1), ~] = subcorr(G(:,:,iDip), Us, 1);
    end
end

if nIter>1
    m = size(G,1);
    for iIter = 2:nIter
        [~,iMax] = max(sc(:,iIter-1));
        if (iIter==2)&&(iSeed>0)
            iMax = iSeed;
        end
        [~, xx] = subcorr(G(:,:,iMax), Us, 1);
        G_max = G(:,:,iMax)*xx;
        % Compute projector
        P = eye(m) - G_max*pinv(G_max'*G_max)*G_max';
        % Project away subspace and leadfields
        Us = P*Us;
        if ~isOutOrient
            parfor iDip=1:nDip
                G(:,:,iDip) = P*G(:,:,iDip);
                sc(iDip,iIter) = subcorr(G(:,:,iDip), Us, 1);
            end
        else
            parfor iDip = 1:nDip
                G(:,:,iDip) = P*G(:,:,iDip);
                [sc(iDip,iIter), x_all(:,iDip,iIter), ~] = subcorr(G(:,:,iDip), Us, 1);
            end
        end
    end
else
    sc = [sc sc];
end

end



