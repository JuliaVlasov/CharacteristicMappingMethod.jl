function c = AdjTimes(A, b)

% A is nxmx4 matrix, b is nxmx2 vector

% compute Adjugate(A)*b

c = cat(3, A(:,:,4).*b(:,:,1) - A(:,:,2).*b(:,:,2),  -A(:,:,3).*b(:,:,1) + A(:,:,1).*b(:,:,2));


return