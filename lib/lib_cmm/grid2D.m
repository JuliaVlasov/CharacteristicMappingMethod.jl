classdef grid2D
    % grid2D Defines a grid structure on a 2D computational domain.
    %   Takes care of all operations and change of variables between
    %   computational domain and reference space.
    %   The computational space is assumed to be [0, LX] x [0, LY] periodic.
    
   properties
      
       gridsize     % size of the matrix (may be obsolete, remove if ok)
       NX           % Number of cells in x-direction (for rectangular)
       NY           % Number of cells in y-direction (for rectangular)
       NG           % Number of grid points
       cellwidth    % Width of each cell [dx dy]
       type         % Type of grid (boundary conditions, currently only allows flat torus)
       boxsize      % rectangular box
       domain
       
   end
    
   methods
      
       function obj = grid2D(gs, dom)
           % Construct grid object. So far, only type = 'uniform-periodic'
           % with nx by ny rectangular grid is implemented.
           
           
           nx = gs(2);
           ny = gs(1);
           type = 'uniform-periodic';
           
%            switch nargin
%                
%                case 2
%                    nx = varargin{1};
%                    ny = varargin{2};
%                    dom = [0 0 1 1];
%                    type = 'uniform-periodic';
%                    
%                case 3
%                    
%                    nx = varargin{1};
%                    ny = varargin{2};
%                    dom = varargin{3};
%                    type = 'uniform-periodic';
% 
%                case 4
%                    
%                    nx = varargin{1};
%                    ny = varargin{2};
%                    dom = varargin{3};
%                    type = varargin{4};
%            end
%            
           
           
%            if nargin == 4
%               type = 'uniform-periodic';
%            end
           
           obj.type = type;
           
           switch type
               case 'uniform-periodic'
                   obj.NX = nx;
                   obj.NY = ny;
                   obj.domain = dom;
%                    x0 = dom(1);
%                    y0 = dom(2);
                   LX = dom(3)-dom(1);
                   LY = dom(4)-dom(2);
                   
                   obj.gridsize = [ny nx];
                   obj.cellwidth = [LX/nx LY/ny];
                   obj.NG = nx*ny;
                   obj.boxsize = [LX LY];
                   
               otherwise
                   disp('grid type unavailable');
           end
           
       end
       
       
       
       function [cinds, xs, ys, dx, dy] = cellcoords(obj, xb, yb, xq, yq)
           % Maps coordinates in computation domain to the reference space
           % [0,1]^2 (unit square). 
           % Input: (xb, yb) n by 1 vectors which specify the which cell is
           % used. (xq, yq) n by nsten vectors. Each row contain the
           % position of the nsten query points RELATIVE to the same row in
           % (xb, yb).
           % Output: cind, nb by 4 the linear indices of the 4 corners of
           % the cell specified by (xb, yb) in order (--, -+, +-, ++).
           % (xs, ys) the positons of (xq, yq) in the reference space
           % [0,1]^2 corresponding to the selected cell. Note, values may
           % lie outside of [0,1]^2 if (xq, yq) are in difference cells
           % from (xb, yb).
           
           
           dom = obj.domain;
           x0 = dom(1);
           y0 = dom(2);
           cw = obj.cellwidth;
           dx = cw(1); dy = cw(2);
           gsize = obj.gridsize;
           nx = obj.NX; ny = obj.NY;
           
           % use base points to choose cell
           xdiv = (xb-x0)./dx;
           ydiv = (yb-y0)./dy;
           xfloor = floor(xdiv);
           yfloor = floor(ydiv);
           
           xc = xdiv - xfloor;
           yc = ydiv - yfloor;
           
           xgprev = xfloor;
           xgnext = xfloor+1;
           xgprev = mod(xgprev, nx) + 1;
           xgnext = mod(xgnext, nx) + 1;
           ygprev = yfloor;
           ygnext = yfloor+1;
           ygprev = mod(ygprev, ny) + 1;
           ygnext = mod(ygnext, ny) + 1;
           
           
           cmm = sub2ind(gsize, ygprev, xgprev);
           cmp = sub2ind(gsize, ygnext, xgprev);
           cpm = sub2ind(gsize, ygprev, xgnext);
           cpp = sub2ind(gsize, ygnext, xgnext);
           
           cinds = [cmm cmp cpm cpp];
           
           xs = xc + xq./dx;
           ys = yc + yq./dy;
           
           
           
       end
       
       
       function lind = shiftLinInd(obj, shiftvec, lind)
           % does the same as circshift(A, shiftvec), outputs linear
           % indexing of the shifted array.
           % useful for defining finite difference sparse matrices
           
           lind = lind(:);
           
           [i,j] = ind2sub(obj.gridsize, lind);
           i = mod(i-1 - shiftvec(1), obj.NY)+1;
           j = mod(j-1 - shiftvec(2), obj.NX)+1;
           
           lind = sub2ind(obj.gridsize, i, j);
           lind = lind(:);
           
       end
       
       function shiftOp = shiftIndsOp(obj, shiftvec)
           
          shiftOp = @(lind) shiftLinInd(obj, shiftvec, lind);
           
       end
      
       
       function [xg, yg] = gridpoints(obj)
           % Lists all grid points in matrix form.
           
           cwidth = obj.cellwidth;
           nx = obj.NX;
           ny = obj.NY;
           dx = cwidth(1);
           dy = cwidth(2);
           
           dom = obj.domain;
           x0 = dom(1);
           y0 = dom(2);
           
           [xg, yg] = meshgrid(0:nx-1, 0:ny-1);
           xg = xg*dx + x0;
           yg = yg*dy + y0;

       end
       
       function [kxn, kyn] = fftGridpoints(obj)
           % actual modes need to be multiplied by 2pi*i
           
           box = obj.boxsize;
           MX = obj.NX;
           MY = obj.NY;
           
           kxn = ( ones(1,MY)'*(mod((1:MX)-ceil(MX/2+1),MX)- floor(MX/2)) )./box(1);
           kyn = ( (mod((1:MY)'-ceil(MY/2+1),MY)-floor(MY/2))*ones(1,MX) )./box(2);

           
       end
       
   end
    
    
end