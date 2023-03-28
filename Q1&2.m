% Question 1
clc; clear; close all;
cvx_begin
cvx_precision high
    variables x_opt(3, 1)
    minimize (exp(x_opt(1) - 1) ...
    + exp(-x_opt(1) + 1) ...
    + exp(x_opt(2) - 2) ...
    + exp(-x_opt(2) + 2) ...
    + exp(x_opt(3) - 3) ...
    + exp(-x_opt(3) + 3) ...
    + sum(x_opt)^4);
cvx_end

disp(['x* = ', num2str([x_opt(1), x_opt(2), x_opt(3)])]);
disp(['f(x*) = ', num2str(cvx_optval)]);

% Function
f = @(x) exp(x(1) - 1) ...
    + exp(-x(1) + 1) ...
    + exp(x(2) - 2) ...
    + exp(-x(2) + 2) ...
    + exp(x(3) - 3) ...
    + exp(-x(3) + 3) ...
    + sum(x)^4;

% Gradient
g = @(x) [
    exp(x(1)-1) ...
        - exp(-x(1) + 1) ...
        + 4*(x(1) ...
        + x(2) ...
        + x(3))^3; 
    exp(x(2) - 2) ...
        - exp(-x(2) + 2) ...
        + 4*(x(1) + x(2) ...
        + x(3))^3; 
    exp(x(3) - 3) ...
        - exp(-x(3) + 3) ...
        + 4*(x(1) + x(2) ...
        + x(3))^3
        ];

h = @(x) 12*(x(1)+x(2)+x(3))^2 + diag([ ...
    exp(-x(1) + 1) + exp(x(1) - 1), ...
    exp(x(2) - 2) + exp(-x(2) + 2), ...
    exp(x(3) - 3) + exp(-x(3) + 3)...
    ]);

K=30;
alpha = 0.1;
beta = 0.5;
x0 = zeros(3,1);
X = zeros(3,K+1);
F = zeros(K+1,1);
xk = x0;
X(:,1) = xk;
F(1) = f(xk);

for k=1:K 
    dxk = -g(xk); 
    dxk = dxk / norm(dxk); 
    tk = alpha;
    while f(xk + tk * dxk) > f(xk)
        tk = beta * tk;
    end
    xk= xk + tk * dxk;
    X(:, k+1) = xk;
    F(k+1) = f(xk);
end

disp(xk);
disp(num2str(norm(xk - x_opt)));
disp(num2str(norm(x_opt)));
disp(['||x(K)-x*|| / ||x*|| = ', num2str(norm(xk - x_opt) / norm(x_opt))]);

fig1 = figure();
box on; hold on;
plot(0:K, abs(F-cvx_optval), 'linewidth', 1.5);
set(gca, 'yscale', 'log', 'ticklabelinterpreter', 'latex');
xlabel('Iteration number, $k$','interpreter', 'latex');
ylabel('Error, $|f(x^{(k)})-f(x^*)|$', 'interpreter', 'latex');
title(['Gradient method with backtracking line search, $\alpha=', ...
    num2str(alpha), '$, $\beta=', num2str(beta),'$'], 'interpreter', 'latex');
grid on
saveas(fig1, 'gradient_1.eps', 'epsc')

fig2 = figure();
box on; hold on;
plot3(X(1,:), X(2,:), X(3,:), 'linewidth', 1.5);
plot3(X(1,1), X(2,1), X(3,1), 'o','linewidth', 1.5);
plot3(x_opt(1), x_opt(2), x_opt(3), 'p', 'linewidth', 1.5);
view([-1, -1, 1]);
legend({'Trajectory of iterates, $x^{(k)}$',...
'Initial iterate, $x^{(0)}$',...
'CVX solution, $x^*$'},'interpreter','latex','location','northwest');
set(gca,'ticklabelinterpreter','latex');
title('Gradient method','interpreter','latex');
saveas(fig2, 'gradient_2.eps', 'epsc')

% Question 2

alpha = 0.1;
beta = 0.5;
xk=x0;
X(:,1) = xk;
F(1) = f(xk);
for k=1:K
    dxk = -h(xk) \ g(xk);
    dxk = dxk / norm(dxk);
    tk = alpha;
    while f(xk + tk * dxk) > f(xk)
        tk = beta * tk;
    end
    xk = xk + tk * dxk;
    X(:, k+1) = xk;
    F(k+1) = f(xk);
end 

disp(['||x(K)-x*|| / ||x*|| = ', num2str(norm(xk - x_opt) / norm(x_opt))]);

fig3 = figure();
box on; hold on;
plot(0:K, abs(F-cvx_optval), 'linewidth', 1.5);
set(gca, 'yscale', 'log', 'ticklabelinterpreter', 'latex');
xlabel('Iteration number, $k$','interpreter', 'latex');
ylabel('Error, $|f(x^{(k)})-f(x^*)|$','interpreter', 'latex');
title(['Newton''s method with backtracking line search, $\alpha=',...
    num2str(alpha), '$, $\beta=', num2str(beta), '$'], ...
    'interpreter', 'latex');
grid on
saveas(fig3, 'newton_1.eps', 'epsc')

fig4 = figure();
box on; hold on;
plot3(X(1,:), X(2,:), X(3,:), 'linewidth', 1.5);
plot3(X(1,1), X(2,1), X(3,1), 'o', 'linewidth', 1.5);
plot3(x_opt(1), x_opt(2), x_opt(3), 'p', 'linewidth', 1.5);
view([-1, -1, 1]);
legend({'Trajectory of iterates, $x^{(k)}$', 'Initial iterate, $x^{(0)}$', 'CVX solution, $x^*$'}, 'interpreter', 'latex', 'location', 'northwest'); 
set(gca, 'ticklabelinterpreter', 'latex');
title('Newton''s method', 'interpreter', 'latex');
saveas(fig4, 'newton_2.eps', 'epsc')











