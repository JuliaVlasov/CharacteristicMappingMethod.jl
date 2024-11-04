function ph = fourierProd23(fh, gh, sz)
% does 2/3 rule by padding 50%, multiply and remove padding

fh = reshape(fh, sz);
gh = reshape(gh, sz);

nsz = ceil(1.5.*sz);

fh = FourierFreqResamp(fh, sz, nsz);
gh = FourierFreqResamp(gh, sz, nsz);

f = real(ifft2(fh)); g = real(ifft2(gh));
ph = fft2(f.*g);

ph = FourierFreqResamp(ph, nsz, sz);

return