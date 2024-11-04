function fjet = HPullback1formjet(dom, gszo, form1jet, gszi, Minjet)
% mapCompose(Xoutjet, Youtjet, Xinjet, Yinjet)
% Finds the Hermite interpolant of the composition of two
% Hermite transformations.
% The outside map (Xout, Yout) must be a Hermite define on the
% current grid, no restriction to the grid for (Xin, Yin).

nx = gszi(2); ny = gszi(1);
[xg, yg] = meshgrid((0:nx-1)*dom(3)/nx+dom(1), (0:ny-1)*dom(4)/ny+dom(2));


EPS = 1e-6;
xg = xg(:); yg = yg(:);
xq = EPS.*[0 -1 -1 1 1];
yq = EPS.*[0 -1 1 -1 1];


[Xq, Yq, Jac] = HMapInterp(dom, gszi, 'cubic', Minjet, xg, yg, xq, yq);
Xb = Xq(:,1); Yb = Yq(:,1);
Xq = Xq(:,2:5)-Xb; Yq = Yq(:,2:5)-Yb;
Jac(:,1,:) = [];

fq = HFunInterp(dom, gszo, 'cubic', form1jet, [0 0], Xb, Yb, Xq, Yq);

fq = MatTimes(Jac, fq, true);
f1 = fq(:,:,1);
f2 = fq(:,:,2);

xcoeff = [-1; -1; 1; 1]./(4*EPS);
ycoeff = [-1; 1; -1; 1]./(4*EPS);
xycoeff = [-1; 1; 1; -1]./(4*EPS.*EPS);

f1x = f1*xcoeff;
f1y = f1*ycoeff;
f1xy = f1*xycoeff;
f1 = sum(f1, 2)/4;


f2x = f2*xcoeff;
f2y = f2*ycoeff;
f2xy = f2*xycoeff;
f2 = sum(f2, 2)/4;

fjet = cat(3, [f1 f1x f1y f1xy], [f2 f2x f2y f2xy]);

end