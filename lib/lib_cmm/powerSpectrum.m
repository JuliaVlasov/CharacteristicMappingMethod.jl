function [pow, kn] = powerSpectrum(fh, gsz)


MX = gsz(2);
MY = gsz(1);
kxn = ( ones(1,MY)'*(mod((1:MX)-ceil(MX/2+1),MX)- floor(MX/2)) );
kyn = ( (mod((1:MY)'-ceil(MY/2+1),MY)-floor(MY/2))*ones(1,MX) );

KN = round(sqrt(kxn.^2 + kyn.^2));

KNmax = max(KN, [], 'all');

kn = (1:KNmax).';

FH = sqrt(sum(abs(fh).^2, 3));

KN = KN(:);
FH = FH(:);

FH(KN == 0) = [];
KN(KN == 0) = [];

pow = accumarray(KN(:), FH(:), [KNmax 1]);
pow = pow./prod(gsz);

return