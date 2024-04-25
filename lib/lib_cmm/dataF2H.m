function fjet = dataF2H(params, fhat, gsz, ngsz)

fhat = FourierFreqResamp(fhat, gsz, ngsz);


MX = ngsz(2);
MY = ngsz(1);
kxn = ( ones(1,MY)'*(mod((1:MX)-ceil(MX/2+1),MX)- floor(MX/2)) );
kyn = ( (mod((1:MY)'-ceil(MY/2+1),MY)-floor(MY/2))*ones(1,MX) );

rx = 2*pi*1i./params.L(1);
ry = 2*pi*1i./params.L(2);

f = ifft2(fhat, 'symmetric');
fx = ifft2(rx.*kxn.*fhat, 'symmetric');
fy = ifft2(ry.*kyn.*fhat, 'symmetric');
fxy = ifft2(rx.*kxn*ry.*kyn.*fhat, 'symmetric');

fjet = cat(3, f, fx, fy, fxy);

fjet = reshape(fjet, [], 4);


return