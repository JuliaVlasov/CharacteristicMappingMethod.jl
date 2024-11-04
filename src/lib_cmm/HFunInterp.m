function fq = HFunInterp(dom, gsz, interpType, fjetList, derivList, xb, yb, xq, yq)

if nargin == 7
   xq = 0;
   yq = 0;
end

grid = grid2D(gsz, dom);
bsz = size(xb);

fq = squeeze(HInterp2(grid, interpType, fjetList, derivList, xb(:), yb(:), xq, yq));

fsz = size(fq);

if bsz(2) ~= 1 && nargin == 7
    fsz = [bsz fsz(2:end)];
end



fq = reshape(fq, fsz);

return