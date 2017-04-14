clear;

function [linfe, ldiago] = factorise(diago, ss_diago)
    n = length(diago);
    for j=1:n
        ldiago($+1) = 0;
    end
    for k=1:n-1
        linfe($+1) = 0;
    end
    ldiago(1) = sqrt(diago(1));
    for i=2:n
        ldiago(i) = sqrt(diago(i) - linfe(i-1)^2);
        linfe(i-1) = ss_diago(i-1)/ldiago(i-1);
    end
endfunction


function Z = descente(linfe, ldiago, Y)
    n = length(ldiago);
    for j=1:n
        Z($+1) = 0;
    end
    Z(1) = Y(1)/ldiago(1);
    for i=2:n
        Z(i) = (Y(i) - linfe(i-1) * Z(i-1)) / ldiago(i);
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


function y = C(x,l)
    y = exp(-x/l)
endfunction


function B = b_create(n)
    for i=1:n
        B($+1) = 0;
    end
    B(1) = C(0.5);
endfunction


function A = A_create(n,l)
    A = zeros(n,n);
    dx = 2*l / (n+1);
    X($+1) = -l;
    for i=1:n+1
        X($+1) = X(i) + dx; 
    end
    A(1,1) = C(2 + 0.5) + C(2 - 0.5);
    A(1,2) = - C(2 + 0.5);
    A(n,n-1) = - C((n+1) - 0.5);
    A(n,n) = C((n+1) + 0.5) + C((n+1) - 0.5);
    for i=2:n-1
        A(i,i-1) = - C((i+1) - 0.5);
        A(i,i) = C((i+1) + 0.5) + C((i+1) - 0.5);
        A(i,i+1) = - C((i+1) + 0.5);
    end
endfunction


function y = sol_exacte(x, l)
    y = exp(1)*(exp(1) - exp(x/l))/(exp(2)- 1);
endfunction


function diago = diagonale(A,n)
    for i=1:n
        diago($+1) = A(i,i);
    end
endfunction


function ss_diago = ss_diagonale(A,n)
    for i=1:n-1
        ss_diago($+1) = A(i+1,i);
    end
endfunction


n = 20;
l = 100;
A = A_create(n,l);
b = b_create(n);
//disp("matrice A: ");
//disp(A);
//disp("vecteur b: ");
//disp(b);
//disp("on veut résoudre Au=b");
diago = diagonale(A,n);
ss_diago = ss_diagonale(A,n);
[linfe, ldiago] = factorise(diago, ss_diago);
Z = descente(linfe, ldiago, b);
X = remonte(linfe, ldiago, Z);
//disp("résultat calculé: ");
//disp(X);
pas = 2*l/n;
X_coord = [-l:pas:l];
for i=1:n
    Y_exact($+1) = sol_exacte(X_coord(i));
end

disp("coord:");
disp(X_coord);
clf; //Efface les courbes déjà présentes
plot2d(sol_exacte(X_coord),style=2);
plot2d(X,style=1);
plot2d(sin(X_coord),style=3);
legend(["Solution exacte";"Solution numérique avec n="+string(n);"ln"]);


















