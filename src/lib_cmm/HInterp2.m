function fq = HInterp2(grid, interpType, fjetList, derivList, xb, yb, xq, yq)
% HInterp2(grid, interpType functionList, derivList, basepoints, stencilpoints)
% output 4D: fq, stencil, derivs, fun

if nargin == 6
   xq = 0;
   yq = 0;
end

ifreal = all(isreal(fjetList)) && all(isreal(xb)) && all(isreal(yb)) && all(isreal(xq)) && all(isreal(yq));

if ~ifreal
    warning('not all data are real numbers, imaginary parts discarded');
    fjetList = real(fjetList);
    xb = real(xb); yb = real(yb);
    xq = real(xq); yq = real(yq);
    keyboard();
end
 
ifFin = all(isfinite(fjetList(:))) && all(isfinite(xb(:))) && all(isfinite(yb(:))) && all(isfinite(xq(:))) && all(isfinite(yq(:)));
if ~ifFin
    error("some inputs are NaN or Inf");
end


% potentially base points going out of bounds causes problems
% also stencil points going too far might be a problem
cw = grid.cellwidth;
ds = max(cw);

if (max(abs(xq(:))) > ds) || (max(abs(yq(:))) > ds)
   warning('stencil points are very far'); 
end


switch lower(interpType)
    
    
    case 'linear'
        %% "Hermite" Linear interp
        
        gsize = grid.gridsize;
        gsize = gsize([2 1]);
                
        if (sum((derivList(:) ~= 0) & (derivList(:) ~= 1) ) > 0)
            error("invalid derivative number");
        end
        
        if ( size(fjetList,1) ~= grid.NG )
            error("grid and data sizes mismatch");
        end
        
        if ( size(fjetList,2) ~= 1 )
            error("grid data not linear interp");
        end
        
        fq = HLINTERP2_MT(grid.domain, gsize, fjetList, derivList, xb(:), yb(:), xq, yq);
%         fq = HLINTERP2(grid.domain, gsize, fjetList, derivList, xb(:), yb(:), xq, yq);



        %% problem with C code, revert to matlab
%         ng = grid.NG;
% 
% 
%         [nderiv, ~] = size(derivList);
%         nfun = size(fjetList, 3);
% 
%         xb = xb(:); yb = yb(:);
%         nb = numel(xb);
%         [~, nsten] = size(xq);
%         % if nsten > 1
%         %     error('no multiple stencil points, just interpolate separately using same base point');
%         % end
% 
%         [cinds, xs, ys, dx, dy] = grid.cellcoords(xb, yb, xq, yq);
%         clear xb yb xq yq
%         xs = xs(:); ys = ys(:);
%         % Each row of cinds corresponds to the cell containing xb, yb
%         % of that row. cinds has 4 columns indexing the gridpoints
%         % corresponding to the corners of this cell. xs, ys are same
%         % size as xq, yq and corresponds to coordinates of these points
%         % in that cell.
% 
% 
%         % Now evaluate 1D basis, for all the derivatives required
% 
%         [xderiv, ~, whichXbasis] = unique(derivList(:,1));
%         % all the required xderivs. whichXbasis will tell us which one
%         % of the evaluated basis to use for that row
%         nxderiv = numel(xderiv);
% 
%         [yderiv, ~, whichYbasis] = unique(derivList(:,2));
%         nyderiv = numel(yderiv);
% 
%         
%         BX = cell(nxderiv, 1);
%             BY = cell(nyderiv, 1);
%             
%             
%             for i = 1:nxderiv
%                 
%                 B = evalDerivLBasis1D(xs, dx, xderiv(i));
%                 BX{i} = B;
%                 
%             end
%             
%             for i = 1:nyderiv
%                 
%                 B = evalDerivLBasis1D(ys, dy, yderiv(i));
%                 BY{i} = B;
%                 
%             end
%             
%             clear xs ys
%             
%             % The above prevents repeated evaluation of the same basis
%             % functions. Still will need a copy of each for tensoring
%             
%             BXfull = cell(nderiv, 1);
%             BXfull(1:nderiv) = BX(whichXbasis(1:nderiv));
%             BYfull = cell(nderiv, 1);
%             BYfull(1:nderiv) = BY(whichYbasis(1:nderiv));
%             
%             clear B dB BX dBX BY dBY
%             
%             
%             % pick all corners for each of the 1 jet data
%             fqdata = fjetList(cinds(:), 1, :);
%             
%             fqdata = reshape(fqdata, [], 1, 4, nfun);
% 
%             fq = zeros(nb, nsten, nderiv, nfun);
%             % dim3 will be contracted
%             
%             
%             % separate sum needed for each derivative
%             for i = 1:nderiv
%                 
%                 % the basis [bb b'b bb' b'b'] for [f fx fy fxy]
%                 B = rowKron22(BXfull{i}, BYfull{i});
%                 
%                 B = reshape(B, [], nsten, 4);
%                 fq(:,:,i,:) = sum(fqdata.*B, 3);
%                                 
%             end
%             fq = reshape(fq, [], nsten, nderiv, nfun);

            

    case 'cubic'
        %% Hermite Cubic interp
        
        gsize = grid.gridsize;
        gsize = gsize([2 1]);
        
        if (sum((derivList(:) ~= 0) & (derivList(:) ~= 1) & (derivList(:) ~= 2) ) > 0)
            error("invalid derivative number");
        end

        if ( size(fjetList,1) ~= grid.NG )
            error("grid and data sizes mismatch");
        end
        
        if ( size(fjetList,2) ~= 4 )
            error("grid data not linear interp");
        end
            
        fq = HCINTERP2_MT(grid.domain, gsize, fjetList, derivList, xb(:), yb(:), xq, yq);
%         fq = HCINTERP2(grid.domain, gsize, fjetList, derivList, xb(:), yb(:), xq, yq);
  
            
    case 'quintic'
        %% Hermite Quintic interp
        
        
        ng = grid.NG;
        
        
        [nderiv, ~] = size(derivList);
        nfun = size(fjetList, 3);
        
        xb = xb(:); yb = yb(:);
        nb = numel(xb);
        [~, nsten] = size(xq);
        % if nsten > 1
        %     error('no multiple stencil points, just interpolate separately using same base point');
        % end
        
        [cinds, xs, ys, dx, dy] = grid.cellcoords(xb, yb, xq, yq);
        clear xb yb xq yq
        xs = xs(:); ys = ys(:);
        % Each row of cinds corresponds to the cell containing xb, yb
        % of that row. cinds has 4 columns indexing the gridpoints
        % corresponding to the corners of this cell. xs, ys are same
        % size as xq, yq and corresponds to coordinates of these points
        % in that cell.
        
        
        % Now evaluate 1D basis, for all the derivatives required
        
        [xderiv, ~, whichXbasis] = unique(derivList(:,1));
        % all the required xderivs. whichXbasis will tell us which one
        % of the evaluated basis to use for that row
        nxderiv = numel(xderiv);
        
        [yderiv, ~, whichYbasis] = unique(derivList(:,2));
        nyderiv = numel(yderiv);
        
            BX = cell(nxderiv, 1); dBX = cell(nxderiv, 1); ddBX = cell(nxderiv, 1);
            BY = cell(nyderiv, 1); dBY = cell(nyderiv, 1); ddBY = cell(nyderiv, 1);
            
            
            for i = 1:nxderiv
                
                [B, dB, ddB] = evalDerivQBasis1D(xs, dx, xderiv(i));
                BX{i} = B;
                dBX{i} = dB;
                ddBX{i} = ddB;
                
            end
            
            for i = 1:nyderiv
                
                [B, dB, ddB] = evalDerivQBasis1D(ys, dy, yderiv(i));
                BY{i} = B;
                dBY{i} = dB;
                ddBY{i} = ddB;
                
            end
            
            clear xs ys
            
            % The above prevents repeated evaluation of the same basis
            % functions. Still will need a copy of each for tensoring
            
            BXfull = cell(nderiv, 1);
            BXfull(1:nderiv) = BX(whichXbasis(1:nderiv));
            dBXfull = cell(nderiv, 1);
            dBXfull(1:nderiv) = dBX(whichXbasis(1:nderiv));
            ddBXfull = cell(nderiv, 1);
            ddBXfull(1:nderiv) = ddBX(whichXbasis(1:nderiv));
            BYfull = cell(nderiv, 1);
            BYfull(1:nderiv) = BY(whichYbasis(1:nderiv));
            dBYfull = cell(nderiv, 1);
            dBYfull(1:nderiv) = dBY(whichYbasis(1:nderiv));
            ddBYfull = cell(nderiv, 1);
            ddBYfull(1:nderiv) = ddBY(whichYbasis(1:nderiv));
            
            clear B dB ddB BX dBX ddBX BY dBY ddBY
            
            
            % pick all corners for each of the 9 jet data
            fqdata = fjetList(cinds(:), (1:9), :);
            fqdata = reshape(fqdata, [], 1, 36, nfun);
            
            
            fq = zeros(nb, nsten, nderiv, nfun);
            
            
            % separate sum needed for each derivative
            for i = 1:nderiv
                
                % the basis [BB dBB BdB ddBB dBdB BddB ddBdB dBddB ddBddB] for [f fx fy fxx fxy fyy fxxy fxyy fxxyy]
                B = [rowKron22(BXfull{i}, BYfull{i}) rowKron22(dBXfull{i}, BYfull{i}) rowKron22(BXfull{i}, dBYfull{i}) rowKron22(ddBXfull{i}, BYfull{i}) rowKron22(dBXfull{i}, dBYfull{i}) ...
                    rowKron22(BXfull{i}, ddBYfull{i}) rowKron22(ddBXfull{i}, dBYfull{i}) rowKron22(dBXfull{i}, ddBYfull{i}) rowKron22(ddBXfull{i}, ddBYfull{i})];
                
                B = reshape(B, [], nsten, 36);
                fq(:,:,i,:) = sum(fqdata.*B, 3);
                
            end
%             fq = squeeze(fq);

            
            
    case 'linear-op'
        %% "Hermite" Linear interp
        
        
        
        
        ng = grid.NG;
        
        
        [nderiv, ~] = size(derivList);
        nfun = size(fjetList, 3);
        
        xb = xb(:); yb = yb(:);
        nb = numel(xb);
        [~, nsten] = size(xq);
        % if nsten > 1
        %     error('no multiple stencil points, just interpolate separately using same base point');
        % end
        
        [cinds, xs, ys, dx, dy] = grid.cellcoords(xb, yb, xq, yq);
        clear xb yb xq yq
        xs = xs(:); ys = ys(:);
        % Each row of cinds corresponds to the cell containing xb, yb
        % of that row. cinds has 4 columns indexing the gridpoints
        % corresponding to the corners of this cell. xs, ys are same
        % size as xq, yq and corresponds to coordinates of these points
        % in that cell.
        
        
        % Now evaluate 1D basis, for all the derivatives required
        
        [xderiv, ~, whichXbasis] = unique(derivList(:,1));
        % all the required xderivs. whichXbasis will tell us which one
        % of the evaluated basis to use for that row
        nxderiv = numel(xderiv);
        
        [yderiv, ~, whichYbasis] = unique(derivList(:,2));
        nyderiv = numel(yderiv);
        
            % Using cells. Test the speed of these two versions
            BX = cell(nxderiv, 1);
            BY = cell(nyderiv, 1);
            
            
            for i = 1:nxderiv
                
                B = evalDerivLBasis1D(xs, dx, xderiv(i));
                BX{i} = B;
                
            end
            
            for i = 1:nyderiv
                
                B = evalDerivLBasis1D(ys, dy, yderiv(i));
                BY{i} = B;
                
            end
            
            clear xs ys
            
            % The above prevents repeated evaluation of the same basis
            % functions. Still will need a copy of each for tensoring
            
            BXfull = cell(nderiv, 1);
            BXfull(1:nderiv) = BX(whichXbasis(1:nderiv));
            BYfull = cell(nderiv, 1);
            BYfull(1:nderiv) = BY(whichYbasis(1:nderiv));
            
            clear B dB BX dBX BY dBY
            
            cinds = repmat(cinds, [nsten 1]);
            
            lincinds = cinds(:);
            
            fq = cell(nderiv, 1);

            
            % separate sum needed for each derivative
            for i = 1:nderiv
                
                B = rowKron22(BXfull{i}, BYfull{i});
                % in linear ordering of B, first list all evaluations of
                % stencils for one basis, then list all 4 basis

                fq{i} = sparse(repmat( (1:nb*nsten).', [4 1]), lincinds, B(:), nb*nsten, ng);
                
                                
            end
            
            if nderiv == 1
               fq = fq{1}; 
            end
            
            
            
            
        
                 
    case 'cubic-op'
        %% Hermite Cubic interp
        
        
        
        
        ng = grid.NG;
        
        
        [nderiv, ~] = size(derivList);
        nfun = size(fjetList, 3);
        
        xb = xb(:); yb = yb(:);
        nb = numel(xb);
        [~, nsten] = size(xq);
        % if nsten > 1
        %     error('no multiple stencil points, just interpolate separately using same base point');
        % end
        
        [cinds, xs, ys, dx, dy] = grid.cellcoords(xb, yb, xq, yq);
        clear xb yb xq yq
        xs = xs(:); ys = ys(:);
        % Each row of cinds corresponds to the cell containing xb, yb
        % of that row. cinds has 4 columns indexing the gridpoints
        % corresponding to the corners of this cell. xs, ys are same
        % size as xq, yq and corresponds to coordinates of these points
        % in that cell.
        
        
        % Now evaluate 1D basis, for all the derivatives required
        
        [xderiv, ~, whichXbasis] = unique(derivList(:,1));
        % all the required xderivs. whichXbasis will tell us which one
        % of the evaluated basis to use for that row
        nxderiv = numel(xderiv);
        
        [yderiv, ~, whichYbasis] = unique(derivList(:,2));
        nyderiv = numel(yderiv);
        
        
            BX = cell(nxderiv, 1); dBX = cell(nxderiv, 1);
            BY = cell(nyderiv, 1); dBY = cell(nyderiv, 1);
            
            
            for i = 1:nxderiv
                
                [B, dB] = evalDerivCBasis1D(xs, dx, xderiv(i));
                BX{i} = B;
                dBX{i} = dB;
                
            end
            
            for i = 1:nyderiv
                
                [B, dB] = evalDerivCBasis1D(ys, dy, yderiv(i));
                BY{i} = B;
                dBY{i} = dB;
                
            end
            
            clear xs ys
            
            % The above prevents repeated evaluation of the same basis
            % functions. Still will need a copy of each for tensoring
            
            BXfull = cell(nderiv, 1);
            BXfull(1:nderiv) = BX(whichXbasis(1:nderiv));
            dBXfull = cell(nderiv, 1);
            dBXfull(1:nderiv) = dBX(whichXbasis(1:nderiv));
            BYfull = cell(nderiv, 1);
            BYfull(1:nderiv) = BY(whichYbasis(1:nderiv));
            dBYfull = cell(nderiv, 1);
            dBYfull(1:nderiv) = dBY(whichYbasis(1:nderiv));
            
            clear B dB BX dBX BY dBY
            
            
            cinds = repmat(cinds, [nsten 1]);
            
            lincinds = cinds(:) + [0 ng 2*ng 3*ng];
            lincinds = lincinds(:);
            
            fq = cell(nderiv, 1);
            

            % separate sum needed for each derivative
            for i = 1:nderiv
                
                % the basis [BB dBB BdB dBdB] for [f fx fy fxy]
                B = [rowKron22(BXfull{i}, BYfull{i}) rowKron22(dBXfull{i}, BYfull{i}) rowKron22(BXfull{i}, dBYfull{i}) rowKron22(dBXfull{i}, dBYfull{i})];
                
                fq{i} = sparse(repmat( (1:nb*nsten).', [16 1]), lincinds, B(:), nb*nsten, 4*ng);
                
                
            end
            
            
            if nderiv == 1
               fq = fq{1}; 
            end
            
            
            
            
            
            
            
            
            
    case 'quintic-op'
        %% Hermite Quintic interp
        
        
        
        
        ng = grid.NG;
        
        
        [nderiv, ~] = size(derivList);
        nfun = size(fjetList, 3);
        
        xb = xb(:); yb = yb(:);
        nb = numel(xb);
        [~, nsten] = size(xq);
        % if nsten > 1
        %     error('no multiple stencil points, just interpolate separately using same base point');
        % end
        
        [cinds, xs, ys, dx, dy] = grid.cellcoords(xb, yb, xq, yq);
        clear xb yb xq yq
        xs = xs(:); ys = ys(:);
        % Each row of cinds corresponds to the cell containing xb, yb
        % of that row. cinds has 4 columns indexing the gridpoints
        % corresponding to the corners of this cell. xs, ys are same
        % size as xq, yq and corresponds to coordinates of these points
        % in that cell.
        
        
        % Now evaluate 1D basis, for all the derivatives required
        
        [xderiv, ~, whichXbasis] = unique(derivList(:,1));
        % all the required xderivs. whichXbasis will tell us which one
        % of the evaluated basis to use for that row
        nxderiv = numel(xderiv);
        
        [yderiv, ~, whichYbasis] = unique(derivList(:,2));
        nyderiv = numel(yderiv);
        
        
            BX = cell(nxderiv, 1); dBX = cell(nxderiv, 1); ddBX = cell(nxderiv, 1);
            BY = cell(nyderiv, 1); dBY = cell(nyderiv, 1); ddBY = cell(nyderiv, 1);
            
            
            for i = 1:nxderiv
                
                [B, dB, ddB] = evalDerivQBasis1D(xs, dx, xderiv(i));
                BX{i} = B;
                dBX{i} = dB;
                ddBX{i} = ddB;
                
            end
            
            for i = 1:nyderiv
                
                [B, dB, ddB] = evalDerivQBasis1D(ys, dy, yderiv(i));
                BY{i} = B;
                dBY{i} = dB;
                ddBY{i} = ddB;
                
            end
            
            clear xs ys
            
            % The above prevents repeated evaluation of the same basis
            % functions. Still will need a copy of each for tensoring
            
            BXfull = cell(nderiv, 1);
            BXfull(1:nderiv) = BX(whichXbasis(1:nderiv));
            dBXfull = cell(nderiv, 1);
            dBXfull(1:nderiv) = dBX(whichXbasis(1:nderiv));
            ddBXfull = cell(nderiv, 1);
            ddBXfull(1:nderiv) = ddBX(whichXbasis(1:nderiv));
            BYfull = cell(nderiv, 1);
            BYfull(1:nderiv) = BY(whichYbasis(1:nderiv));
            dBYfull = cell(nderiv, 1);
            dBYfull(1:nderiv) = dBY(whichYbasis(1:nderiv));
            ddBYfull = cell(nderiv, 1);
            ddBYfull(1:nderiv) = ddBY(whichYbasis(1:nderiv));
            
            clear B dB ddB BX dBX ddBX BY dBY ddBY
            
            
            cinds = repmat(cinds, [nsten 1]);
            
            lincinds = cinds(:) + (0:8)*ng;
            lincinds = lincinds(:);
            
            fq = cell(nderiv, 1);
            

            % separate sum needed for each derivative
            for i = 1:nderiv
                
                % the basis [BB dBB BdB ddBB dBdB BddB ddBdB dBddB ddBddB] for [f fx fy fxx fxy fyy fxxy fxyy fxxyy]
                B = [rowKron22(BXfull{i}, BYfull{i}) rowKron22(dBXfull{i}, BYfull{i}) rowKron22(BXfull{i}, dBYfull{i}) rowKron22(ddBXfull{i}, BYfull{i}) rowKron22(dBXfull{i}, dBYfull{i}) ...
                    rowKron22(BXfull{i}, ddBYfull{i}) rowKron22(ddBXfull{i}, dBYfull{i}) rowKron22(dBXfull{i}, ddBYfull{i}) rowKron22(ddBXfull{i}, ddBYfull{i})];
                
                
                fq{i} = sparse(repmat( (1:nb*nsten).', [36 1]), lincinds, B(:), nb*nsten, 8*ng);

            end
            
            if nderiv == 1
               fq = fq{1}; 
            end
        
        
end








end


%% Helper funcs

%% 2D kronecker for the columns
function C = rowKron22(A,B)
% row kronecker tensor of two matrices with 2 columns each. Corresponds to
% kronecker tensor applied to each row.


% its job is to evaluate the basis associated to each of the corners of the
% cell

[nrA, ncA] = size(A);
[nrB, ncB] = size(B);

if nrA ~= nrB
    error('A and B must have same number of rows');
end

if ncA ~= 2 || ncB ~= 2
    error('A and B must each have 2 columns');
end


C = [A(:,1).*B A(:,2).*B];

end


%% For Linear
function B = evalDerivLBasis1D(xp, dx, deriv)

xp = xp(:);

% deriv specifies which derivative of the basis we evaluate
switch deriv
    
    
    case 0
        
        B = [(1-xp) xp];

    case 1
        B = [-ones(numel(xp),1) ones(numel(xp),1)]./dx;
        
    otherwise
        
        B = zeros(numel(xp),2);
        
        disp('Linear has only derivatives 0 and 1');
        
end



end




%% For Cubics
function [B, dB] = evalDerivCBasis1D(xp, dx, deriv)

xp = xp(:);

% deriv specifies which derivative of the basis we evaluate
switch deriv
    
    
    case 0
        
        omt2 = (1-xp).^2;
        t2 = xp.^2;
        
        h0p = (1+2*xp).*omt2;
        p0p = xp.*omt2;
        h1p = t2.*(3-2*xp);
        p1p = t2.*(xp-1);
        
        B = [h0p h1p];
        dB = [p0p p1p]*dx;

    case 1
        
        st = 6*xp;
        omt = 1-xp;
        ttmo = 3*xp-1;
        
        h0p = -st.*omt;
        p0p = -omt.*ttmo;
        h1p = -h0p;
        p1p = xp.*(ttmo-1);
        
        B = [h0p h1p]./dx;
        dB = [p0p p1p];
        
    case 2
        
        
        op = 6*xp-3;
        
        h0p = 2*op;
        p0p = op-1;
        h1p = -h0p;
        p1p = op+1;
        
        B = [h0p h1p]./(dx^2);
        dB = [p0p p1p]./dx;
        
    case 3
        
        on = ones(size(xp));
        h0p = 12*on;
        p0p = 6*on;
        h1p = -12*on;
        p1p = 6*on;
        
        B = [h0p h1p]./(dx^3);
        dB = [p0p p1p]./(dx^2);
        
     
    otherwise
        
        B = [0 0];
        dB = [0 0];
        
        disp('Hermite cubics has only derivatives 0, 1, 2, 3');
        
end



end


%% for Quintics
function [B, dB, ddB] = evalDerivQBasis1D(xp, dx, deriv)

xp = xp(:);

% deriv specifies which derivative of the basis we evaluate
switch deriv
    
    
    case 0
        
        % take the powers first to save computeations
        
        t = xp;
        t2 = xp.^2;
        t3 = xp.^3;
        t4 = xp.^4;
        t5 = xp.^5;
        
        
        h0 =  1 - 10*t3 + 15*t4 - 6*t5;
        h1 =  10*t3 - 15*t4 + 6*t5;
        p0 =  t - 6*t3 + 8*t4 - 3*t5;
        p1 = -4*t3 + 7*t4 - 3*t5;
        q0 = 0.5*t2 - 1.5*t3 + 1.5*t4 - 0.5*t5;
        q1 = 0.5*t3 - t4 + 0.5*t5;
        
        B = [h0 h1];
        dB = [p0 p1]*dx;
        ddB = [q0 q1]*dx^2;

    case 1
        
        t = xp;
        t2 = xp.^2;
        t3 = xp.^3;
        t4 = xp.^4;
        
        
        h0 =  -30*t2 + 60*t3 - 30*t4;
        h1 =  -h0;
        p0 =  1 - 18*t2 + 32*t3 - 15*t4;
        p1 = -12*t2 + 28*t3 - 15*t4;
        q0 = t - 4.5*t2 + 6*t3 - 2.5*t4;
        q1 = 1.5*t2 - 4*t3 + 2.5*t4;
        
        B = [h0 h1]./dx;
        dB = [p0 p1];
        ddB = [q0 q1]*dx;
        
    case 2
        
        t = xp;
        t2 = xp.^2;
        t3 = xp.^3;
        
        
        h0 =  -60*t + 180*t2 - 120*t3;
        h1 =  -h0;
        p0 =  -36*t + 96*t2 - 60*t3;
        p1 = -24*t + 84*t2 - 60*t3;
        q0 = 1 - 9*t + 18*t2 - 10*t3;
        q1 = 3*t - 12*t2 + 10*t3;
        
        B = [h0 h1]./(dx^2);
        dB = [p0 p1]./dx;
        ddB = [q0 q1];
        
    case 3
        
        t = xp;
        t2 = xp.^2;        
        
        h0 =  -60 + 360*t - 360*t2;
        h1 =  -h0;
        p0 =  -36 + 192*t - 180*t2;
        p1 = -24 + 168*t - 180*t2;
        q0 = -9 + 36*t - 30*t2;
        q1 = 3 - 24*t + 30*t2;
        
        B = [h0 h1]./(dx^3);
        dB = [p0 p1]./(dx^2);
        ddB = [q0 q1]./dx;
        
     
    otherwise
        
        B = [0 0];
        dB = [0 0];
        ddB = [0 0];
        
        disp('I only coded up to 3 derivatives for quintics');
        
end


end

