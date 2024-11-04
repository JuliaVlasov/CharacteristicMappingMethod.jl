function varargout = HMapInterp(dom, gsz, interpType, Mjet, xb, yb, xq, yq)
% mapInterp(Xjdata, Yjdata, xb, yb, xq, yq)
% Interpolates the transformation from the displacement data.
% Interpolates the transformation whose DISPLACEMENT data is
% (Xjetdata, Yjetdata), that is, the actual transformation data
% is (x + Xjetdata, y + Yjetdata)
% Outputs the evaluation of the transformation (x+X, y+Y) at
% query points


grid = grid2D(gsz, dom);
bsz = size(xb) + [0 0];

if nargin == 6
    xq = 0;
    yq = 0;


    if nargout == 2

        fq = HInterp2(grid, interpType, Mjet, [0 0], xb(:), yb(:));
        varargout{1} = reshape(fq(:,:,:,1), bsz) + xb;
        varargout{2} = reshape(fq(:,:,:,2), bsz) + yb;

    elseif nargout == 3

        fq = HInterp2(grid, interpType, Mjet, [0 0; 1 0; 0 1], xb(:), yb(:));
        varargout{1} = reshape(fq(:,:,1,1), bsz) + xb;
        varargout{2} = reshape(fq(:,:,1,2), bsz) + yb;
        varargout{3} = cat(3, reshape(fq(:,:,[2 3], 1), [bsz 2]), reshape(fq(:,:,[2 3], 2), [bsz 2])) + cat(3, 1, 0, 0, 1);

    end

else

    if nargout == 2

        fq = HInterp2(grid, interpType, Mjet, [0 0], xb(:), yb(:), xq, yq);
        varargout{1} = fq(:,:,:,1) + xb(:) + xq;
        varargout{2} = fq(:,:,:,2) + yb(:) + yq;

    elseif nargout == 3

        fq = HInterp2(grid, interpType, Mjet, [0 0; 1 0; 0 1], xb(:), yb(:), xq, yq);
        varargout{1} = fq(:,:,1,1) + xb(:) + xq;
        varargout{2} = fq(:,:,1,2) + yb(:) + yq;
        varargout{3} = cat(3, fq(:,:,[2 3], 1), fq(:,:,[2 3], 2)) + cat(3, 1, 0, 0, 1);

    end

end


end
