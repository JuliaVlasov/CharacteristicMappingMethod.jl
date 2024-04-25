function [kx, ky, kn, rx, ry, kxn, kyn] = genWaveNumbers(dom, gsz)

MX = gsz(2);
MY = gsz(1);
kxn = ( ones(1,MY)'*(mod((1:MX)-ceil(MX/2+1),MX)- floor(MX/2)) );
kyn = ( (mod((1:MY)'-ceil(MY/2+1),MY)-floor(MY/2))*ones(1,MX) );

rx = 2*pi*1i./dom(3);
ry = 2*pi*1i./dom(4);

kx = kxn.*rx;
ky = kyn.*ry;

kn = sqrt(abs(kx).^2+abs(ky).^2);

return