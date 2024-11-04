function pc(params,datafield)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pcolor plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcolor(params.X,params.V,datafield);
shading flat
xlabel("$x$")
ylabel("$v$")