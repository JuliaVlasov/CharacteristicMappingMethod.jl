function c = MatTimes(A, b, ifT)

% A is nxmx4 matrix, b is nxmx2 vector

if nargin == 2
    ifT = false;
end

if(~ifT)
    c = cat(3, A(:,:,1).*b(:,:,1) + A(:,:,2).*b(:,:,2),  A(:,:,3).*b(:,:,1) + A(:,:,4).*b(:,:,2));
else
    c = cat(3, A(:,:,1).*b(:,:,1) + A(:,:,3).*b(:,:,2),  A(:,:,2).*b(:,:,1) + A(:,:,4).*b(:,:,2));
end
    


return