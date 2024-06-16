% Definição da função de Beale para plotagem
bealeFunctionPlot = @(x, y) (1.5 - x + x.*y).^2 + ...
                            (2.25 - x + x.*y.^2).^2 + ...
                            (2.625 - x + x.*y.^3).^2;

% Criação de uma grade de pontos
[xGrid, yGrid] = meshgrid(-4:0.1:4, -4:0.1:4);

% Avaliação da função de Beale na grade de pontos
zGrid = bealeFunctionPlot(xGrid, yGrid);

% Criação do gráfico de superfície
figure;
surf(xGrid, yGrid, zGrid);
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('Função de Beale');
colorbar;
shading interp;
colormap(gray); % Usar escala de cinza
