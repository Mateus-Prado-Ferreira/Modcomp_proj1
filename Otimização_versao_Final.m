clear, clc, clf
syms x y %definindo como variáveis simbólicas
f = (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2; %Definindo a função de Beale
G = gradient(f);
H = hessian(f);

% Esquema de definição do tamanho de passo
% 1 - passo fixo
% 2 - passo descendente
% 3 - Armijo
% 4 - Goldstein
% 5 - Minimização Unidimensional
passo = 5;

% Esquema de definicao da direção
% 1 - gradiente
% 2 - Newton
direcao = 1;

% Plotagem da função
hold on
fsurf(f, [-4 4 -4 4]);

% Definir um ponto de partida
r = [2.5, 0.4];
scatter3(r(1), r(2), subs(f, [x y], [r(1) r(2)]), 100, 'g', 'filled');
xlabel('x'); ylabel('y'); zlabel('f(x,y)')

% Procedimento iterativo
i = 0;
tol = 0.0000001;
nmi = 150;  %Número máximo de iterações
df = realmax;
r_old = r;
grad_old = double(subs(G, [x y], [r(1), r(2)]));
d = -grad_old;

while abs(df) > tol
        switch direcao
            case 1  % gradiente
                d = -subs(G, [x y], [r(1) r(2)]);
            case 2  % Newton
                d = -subs(H\G, [x y], [r(1) r(2)]);
        end
        d = double(d) / norm(double(d));
    switch passo
        case 1  % passo fixo
            a = 0.0025;
        case 2  % passo descendente
            a = 0.3 * 0.99^i;
        case 3  % Armijo
            a = 0.03;
            sigma = 0.5;
            theta = 0.9;
            A = double(subs(f, [x y], [r(1) + a * d(1), r(2) + a * d(2)]));
            B = double(subs(f, [x y], [r(1), r(2)]));
            C = double(subs(G, [x y], [r(1), r(2)])' * d);
            while (A - B) > sigma * a * C
                a = theta * a;
                A = double(subs(f, [x y], [r(1) + a * d(1), r(2) + a * d(2)]));
                B = double(subs(f, [x y], [r(1), r(2)]));
                C = double(subs(G, [x y], [r(1), r(2)])' * d);
            end
        case 4  % Goldstein
            a = 0.05;
            sigma1 = 0.1;
            sigma2 = 1.1;
            theta = 0.9;
            A = double(subs(f, [x y], [r(1) + a * d(1), r(2) + a * d(2)]));
            B = double(subs(f, [x y], [r(1), r(2)]));
            C = double(subs(G, [x y], [r(1), r(2)])' * d);
            while (A - B > sigma1 * a * C) || (A - B < sigma2 * a * C)
                a = theta * a;
                A = double(subs(f, [x y], [r(1) + a * d(1), r(2) + a * d(2)]));
                B = double(subs(f, [x y], [r(1), r(2)]));
                C = double(subs(G, [x y], [r(1), r(2)])' * d);
            end
        case 5  % Minimização Unidimensional (Steepest Descent)
            amin = 0.05;
            amax = 3;
            tol2 = tol;
            df2 = realmax;
            while abs(df2) > tol2 && (amax - amin) > tol2
                amed = (amin + amax) / 2;
                A = double(subs(f, [x y], [r(1) + amin * d(1), r(2) + amin * d(2)]));
                B = double(subs(f, [x y], [r(1) + amed * d(1), r(2) + amed * d(2)]));
                if B <= A
                    amin = amed;
                else
                    amax = amed;
                end
                df2 = abs(A - B);
            end
            a = amed;
    end
    
    rold = r;
    fold = double(subs(f, [x y], [r(1), r(2)]));
    r = r + a * double(d);
    fnew = double(subs(f, [x y], [r(1), r(2)]));
    df = fnew - fold;
    
    if i == nmi || df >= 0
        r = rold;
        fnew = fold;
        break
    end
    
    scatter3(r(1), r(2), subs(f, [x y], [r(1), r(2)]), 15, 'y', 'filled');
    pause(0.05);
    i = i + 1;
    
    disp(['Iteração ' num2str(i) ': f = ' num2str(fnew) ', x = ' num2str(r(1)) ', y = ' num2str(r(2))]);
end

scatter3(r(1), r(2), fnew, 100, 'r', 'filled');
hold off

%Exibindo as métricas na command Window:

disp(['Número total de iterações: ' num2str(i)]);
disp(['Valor final da função: ' num2str(fnew)]);
disp(['Coordenadas finais: x = ' num2str(r(1)) ', y = ' num2str(r(2))]);
