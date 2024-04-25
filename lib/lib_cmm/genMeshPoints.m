function [xg, yg, x, y] = genMeshPoints(dom, gsz)

x = dom(1) + (dom(3)-dom(1))*(0:gsz(2)-1)./gsz(2);
y = dom(2) + (dom(4)-dom(2))*(0:gsz(1)-1)./gsz(1);

[xg, yg] = meshgrid(x,y);

return