clear;

// Factorisation de Cholesky: obtention de L telle que A = L*transp(L)
function [linfe, ldiago] = factorise(diago, ss_diago)
    n = length(diago);
    ldiago($+1) = sqrt(diago(1));
    linfe($+1) = ss_diago(1) / ldiago(1);
    for i=2:n-1
        ldiago($+1) = sqrt(diago(i) - linfe(i-1)*linfe(i-1));
        linfe($+1) = ss_diago(i) / ldiago(i);
    end
    ldiago($+1) = sqrt(diago(n) - linfe(n-1)*linfe(n-1));
endfunction

// Descente: resolution de LZ = B
function Z = descente(linfe, ldiago, Y)
    n = length(ldiago);
    Z($+1) = Y(1)/ldiago(1);
    for i=2:n
        Z($+1) = (Y(i) - linfe(i-1) * Z(i-1))/ldiago(i);
    end
endfunction

// Remonte: resolution de transp(L)U = Z    (U=X=ce que l'on cherche)
function X = remonte(linfe, ldiago, Z)
    n = length(ldiago);
    for j=1:n
        X($+1) = 0;
    end
    X(n) = Z(n)/ldiago(n)
    for i=n-1:-1:1
        X(i) = (Z(i) - linfe(i) * X(i+1))/ldiago(i)
    end
endfunction

A = [[1,-1,0,0,0];[-1,2,-1,0,0];[0,-1,2,-1, 0];[0,0,-1,2,-1];[0,0,0,-1,2]];
disp("matrice A: ");
disp(A);
Y = [3,4,5,6,7];
disp("vecteur b: ");
disp(Y);
disp("on veut résoudre Ax=b");
diago = [1,2,2,2,2];
ss_diago = [-1,-1,-1,-1];
[linfe, ldiago] = factorise(diago, ss_diago);
Z = descente(linfe, ldiago, Y);
X = remonte(linfe, ldiago, Z);
disp("résultat calculé: ");
disp(X);
res = [65;62;55;43;25];
disp("résultat attendu: ");
disp(res);
