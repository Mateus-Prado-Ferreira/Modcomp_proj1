clear, clc, clf
syms x y a %Definindo as variáveis de forma simbólica
f = (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2; %F de Beale
G = gradient(f);

% Plot da função
hold on
fsurf(f, [-4 4 -4 4]);

% Definindo um ponto de partida
r = [0, 0];
scatter3(r(1), r(2), double(subs(f, {x, y}, {r(1), r(2)})), 100, 'g', 'filled');
xlabel('x'); ylabel('y'); zlabel('f(x,y)')

%Variáveis da iteração
i = 0;
tol = 1e-5;
nmi = 150; %Número máximo de iterações
df = realmax;

while abs(df) > tol && i < nmi
    % Calculando o gradiente no ponto atual
    grad = double(subs(G, {x, y}, {r(1), r(2)}));
    
    % Determinando a direção de descida como a oposta ao gradiente
    d = -grad;
    
    % Definindo a função que dá o passo ótimo
    phi = subs(f, {x, y}, {r(1) + a*d(1), r(2) + a*d(2)});
    phi = matlabFunction(phi, 'vars', a); % Convertendo para função anônima no MATLAB
    
    % Otimizando o tamanho do passo usando fminbnd
    a_opt = fminbnd(phi, 0, 1);
    
    % Atualizaando a ponto de partida
    rold = r;
    fold = double(subs(f, {x, y}, {r(1), r(2)}));
    r = r + a_opt * d;
    fnew = double(subs(f, {x, y}, {r(1), r(2)}));
    df = fnew - fold;
    
    % Plotando o ponto atual
    scatter3(r(1), r(2), fnew, 15, 'y', 'filled');
    pause(0.05);
    
    i = i + 1;
    
    disp(['Iteração ' num2str(i) ': f = ' num2str(fnew) ', x = ' num2str(r(1)) ', y = ' num2str(r(2))]);
end

scatter3(r(1), r(2), fnew, 100, 'r', 'filled');
hold off

disp(['Número total de iterações: ' num2str(i)]);
disp(['Valor final da função: ' num2str(fnew)]);
disp(['Coordenadas finais: x = ' num2str(r(1)) ', y = ' num2str(r(2))]);
