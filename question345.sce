function [linfe, ldiago] = factorise(diago, ss_diago)
    ldiago($+1) = sqrt(diago(1));
    for i=2:length(diago)
        ldiago($+1) = sqrt(diago(i) - ss_diago(i-1)^2);
        linfe($+1) = ss_diago(i-1)/diago(i-1);
    end
endfunction


function Z = descente(linfe, ldiago, Y)
    Z($+1) = Y(1)/ldiago(1);
    for i=2:length(ldiago)
        Z($+1) = (Y(i) - linfe(i-1) * Z(i-1)) / ldiago(i);
    end
endfunction


function X = remonte(linfe, ldiago, Z)
    n = length(ldiago);
    for j = 1:n
        X($+1) = 0;
    end
    X(n) = Z(n)/ldiago(n);
    for i=n-1:-1:1
        X(i) = (Z(i) - linfe(i) * X(i+1)) / ldiago(i);
    end
endfunction

A = [[2,1,0];[1,3,1];[0,1,2]];
disp("matrice A: ");
disp(A);
Y = [3,4,5]
disp("vecteur b: ");
disp(Y);
diago = [2,3,2];;
ss_diago = [1,1]
[linfe, ldiago] = factorise(diago, ss_diago);
Z = descente(linfe, ldiago, Y);
X = remonte(linfe, ldiago, Z);
disp("résultat: ");
disp(X);
