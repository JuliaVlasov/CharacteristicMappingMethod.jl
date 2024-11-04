function Mcjet = HCompose(params, gszo, Outjet, gszi, Minjet)
% mapCompose(Xoutjet, Youtjet, Xinjet, Yinjet)
% Finds the Hermite interpolant of the composition of two
% Hermite transformations.
% The outside map (Xout, Yout) must be a Hermite define on the
% current grid, no restriction to the grid for (Xin, Yin).

nx = gszi(2); ny = gszi(1);
[xg, yg] = meshgrid((0:nx-1)*params.L(1)/nx+params.dom(1), (0:ny-1)*params.L(2)/ny+params.dom(2));

Xq = Minjet(:,1,1) + xg(:);
Yq = Minjet(:,1,2) + yg(:);
Xxq = Minjet(:,2,1)+1; Xyq = Minjet(:,3,1); Xxyq = Minjet(:,4,1);
Yxq = Minjet(:,2,2); Yyq = Minjet(:,3,2)+1; Yxyq = Minjet(:,4,2);

derivList = [0 0; 1 0; 0 1; 2 0; 1 1; 0 2];
Mod2q = HFunInterp(params.dom, gszo, 'cubic', Outjet, derivList, Xq, Yq);

Mc = Mod2q(:,1,:);
Mcx = Mod2q(:,2,:).*Xxq + Mod2q(:,3,:).*Yxq;
Mcy = Mod2q(:,2,:).*Xyq + Mod2q(:,3,:).*Yyq;
Mcxy = Mod2q(:,4,:).*Xxq.*Xyq + Mod2q(:,5,:).*(Xxq.*Yyq + Xyq.*Yxq) + Mod2q(:,6,:).*Yxq.*Yyq ...
    + Mod2q(:,2,:).*Xxyq + Mod2q(:,3,:).*Yxyq;

Mcjet = cat(2, Mc, Mcx, Mcy, Mcxy);


end