clc; clear; close all;
% setup 
n = 1300;
A = rand(n, n);
A = 1/2 * (A+A.');
f = @(x) (-(ones(n,1).'))*(log(5 - x.^2) + log(1 + A*x));
g = @(x,A) ((-2*x)./(5-x.^2)) - ((ones(n,1).')*(A./((1+A*x)*ones(n,1).'))).';
h = @(x) (diag(10 + 2*x.^2)./((5-x.^2).^2)) - (-A' * diag(1./(1+A*x).^2) * A);
cvx_begin
cvx_precision high
    variables x_opt(n, 1)
    minimize (-(ones(n,1).')*(log(5 - x_opt.^2) + log(1 + A*x_opt)))
cvx_end

disp(['f(x*) = ', num2str(cvx_optval)]);

K = 80;
alpha = 0.1;
beta = 0.5;
x0 = zeros(n,1);
X = zeros(n,K+1);
F = zeros(K+1,1);
xk=x0;
X(:,1) = xk;
F(1) = f(xk);
for k=1:K
    dxk = -h(xk) \ g(xk,A);
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