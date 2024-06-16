clear, clc, clf
syms x y
f = (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2;
G = gradient(f);
H = hessian(f);

% Plot da função
hold on
fsurf(f, [-4 4 -4 4]);

%Ponto de partida (chute inicial)
r = [2, 0.25];
scatter3(r(1), r(2), double(subs(f, [x y], r)), 100, 'g', 'filled');
xlabel('x'); ylabel('y'); zlabel('f(x,y)')

% Variáveis da iteração
i = 0;
tol = 1e-5;
nmi = 150;
df = realmax;
grad_old = double(subs(G, [x y], r));
d = -grad_old;

while abs(df) > tol && i < nmi
    if i > 0
        grad = double(subs(G, [x y], r));
        beta = (grad' * grad) / (grad_old' * grad_old);
        d = -grad + beta * d;
        grad_old = grad;
    else
        d = -grad_old;
    end
    
    % Determinação do tamanho de passo para Gradiente Conjugado
    hessian = double(subs(H, [x y], r));
    a = (grad_old' * grad_old) / (d' * hessian * d);
    
    % Atualizar a posição
    r_new = r + a * d';
    
    % Calcular novo valor da função e a diferença
    f_old = double(subs(f, [x y], r));
    f_new = double(subs(f, [x y], r_new));
    df = f_new - f_old;
    
    % Atualizando o ponto de partida r
    r = r_new;
    
    % contador
    i = i + 1;
    
    % Plotagem
    scatter3(r(1), r(2), f_new, 15, 'y', 'filled');
    pause(0.05);
    
    % Display
    disp(['Iteração ' num2str(i) ': f = ' num2str(f_new) ', x = ' num2str(r(1)) ', y = ' num2str(r(2))]);
end

scatter3(r(1), r(2), f_new, 100, 'r', 'filled');
hold off

%Exibição dos resultados

disp(['Número total de iterações: ' num2str(i)]);
disp(['Valor final da função: ' num2str(f_new)]);
disp(['Coordenadas finais: x = ' num2str(r(1)) ', y = ' num2str(r(2))]);
