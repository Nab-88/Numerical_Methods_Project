clear;
clf;

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


// La fonction C donnée dans l'énoncé de la question 7
function y = fC(x,l)
    y = exp(-x/l);
endfunction


// Création du vecteur B du système AU=B
function B = b_create(n,l)
    for i=1:n
        B($+1) = 0;
    end
    pas = 2*l/(n+1)
    B(1) = fC(-l + pas/2,l);
endfunction


// Creation de la matrice A du système AU=B
function A = A_create(n,l)
    A = zeros(n,n);
    dx = 2*l /(n+1);
    X($+1) = -l;
    for i=1:n+1
        X($+1) = X(i) + dx; 
    end
    A(1,1) = fC(X(2) + dx/2, l) + fC(X(2) - dx/2, l);
    A(1,2) = - fC(X(2) + dx/2, l);
    A(n,n-1) = - fC(X(n+1) - dx/2, l);
    A(n,n) = fC(X(n+1) + dx/2, l) + fC(X(n+1) - dx/2, l);
    for i=2:n-1
        A(i,i-1) = - fC(X(i+1) - dx/2, l);
        A(i,i) = fC(X(i+1) + dx/2, l) + fC(X(i+1) - dx/2, l);
        A(i,i+1) = - fC(X(i+1) + dx/2, l);
    end
endfunction


// La solution exacte du systeme (13) trouvée à la main: voir compte-rendu
function y = sol_exacte(x, l)
    y = exp(1)*(exp(1) - exp(x/l))/(exp(2)- 1);
endfunction


// Extraction de la diagonale de A
function diago = diagonale(A,n)
    for i=1:n
        diago($+1) = A(i,i);
    end
endfunction


// Extraction de la sous-diagonale de A
function ss_diago = ss_diagonale(A,n)
    for i=1:n-1
        ss_diago($+1) = A(i+1,i);
    end
endfunction


// Résolution du système AX=B en fonction de n et l
function X = resol_numerique(n,l)
    A = A_create(n,l);
    b = b_create(n,l);
    diago = diagonale(A,n);
    ss_diago = ss_diagonale(A,n);
    [linfe, ldiago] = factorise(diago, ss_diago);
    Z = descente(linfe, ldiago, b);
    X = remonte(linfe, ldiago, Z);
endfunction


// Tracé des points Y en fonction de X dans la couleur choisie, en prenant
// en compte un décalage (par exemple commencer à x1 au lieu de x0
function trace_points(X,Y,decalage,n,couleur)
    for i=1:n
        plot(X(i+decalage),Y(i), couleur);
    end
endfunction


// Début du "main"

couleurs = [".r",".g",".y",".b", ".c"]; // les couleurs que l'on va utiliser
iterations = 1; // pour choisir les couleurs
l = 13; // choix du parametre l
debut = 100; // choix du plus petit n
fin = 3000; // choix du plus grand n


// Pour avoir les courbes (exacte et numerique) => 0 (false)
// Pour avoir la courbe de la norme inf de la difference => 1 (true)
Courbe_norme = 1;

for p=debut:100:fin
    n = p;
    X = resol_numerique(n,l); // on obtient X = (x1, ..., xn)
    pas = 2*l/(n+1);
    X_coord = -l:pas:l; // de x0 a x(n+1)
    for j=1:length(X)
        X_exacte($+1) = sol_exacte(X_coord(j),l);
    end
    if Courbe_norme == 0 then
        if p == fin then
            plot(X_coord(2:n+1), X_exacte, "k"); // Tracé de la solution exacte
            legende($+1) = "Solution exacte";
        end
        trace_points(X_coord, X, 1, n, couleurs(pmodulo(iterations,5) + 1)); // Tracé de la sol numerique
        legende($+1) = "Solution numérique avec n="+string(n); // ajout de la legende
    end
    iterations = iterations + 1; // pour changer de couleur
    // Prise en compte de la norme de la difference entre sol exacte et sol numerique
    if Courbe_norme == 1 then
        valeur = norm(X-X_exacte, 'inf');
        abscisse($+1) = n;
        ordonnee($+1) = valeur;
    end
    X_exacte = []; // reinitialisation de X_exacte
end
// Tracé de la courbe norme de la difference
if Courbe_norme == 1 then
    plot(abscisse, ordonnee, "r");
    legende($+1) = "Norme infinie de la différence en fontcion de n"; // ajout de la legende
end
legend(legende); // affichage de la légende










