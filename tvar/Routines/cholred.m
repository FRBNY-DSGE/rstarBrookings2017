function L = cholred(A)

L = A;
[n nn]=size(L);

k = 1;

while k<=n
    if L(k,k)>(1e-3)*mean(diag(A))
    L(k,k)=abs(sqrt(L(k,k)));
    L(k+1:n,k)=L(k+1:n,k)/L(k,k);
    for j=k+1:n
        L(j:n,j)=L(j:n,j)-L(j,k)*L(j:n,k);
    end
    k = k+1;
    else
        k = 1e+6;
    end;
end

L = tril(L)';
