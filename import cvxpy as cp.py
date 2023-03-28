import cvxpy as cp
import numpy as np

# Generate random data
n = 500
a = np.random.randn(n)

# Define the optimization variables
x = cp.Variable(n)

# Define the objective function
objective = cp.Minimize(-cp.sum(cp.log(5 - cp.square(x))) - cp.sum(cp.log(1 + a.T @ x)))

# Define the constraints
constraints = []

# Define the problem
problem = cp.Problem(objective, constraints)

# Solve the problem
problem.solve()

# Print the optimal value and the optimal solution
print("Optimal value:", problem.value)
print("Optimal solution:", x.value)



# # Define the function and its gradient
# def f(x):
#     return -np.sum(np.log(5 - np.square(x))) - np.sum(np.log(1 + a.T @ x))

# def grad_f(x):
#     return -2 * x / (5 - np.square(x)) - a / (1 + a.T @ x)

# # Define the Hessian
# def hessian_f(x):
#     n = len(x)
#     hessian = np.zeros((n, n))
#     for i in range(n):
#         hessian[i, i] = 2 / (np.square(5 - x[i]) + np.square(a[i] / (1 + a.T @ x)))
#     return hessian

# # Define the backtracking line search
# def backtracking_line_search(x, direction, alpha=0.1, beta=0.5):
#     t = 1.0
#     while f(x + t * direction) > f(x) + alpha * t * grad_f(x).T @ direction:
#         t *= beta
#     return t

# # Set the initial point and the tolerance
# x = np.zeros(n)
# tol = 1e-6

# # Set the maximum number of iterations
# max_iter = 1000

# # Implement the Newton's method
# for k in range(max_iter):
#     # Compute the direction and the step size
#     direction = -np.linalg.solve(hessian_f(x), grad_f(x))
#     step_size = backtracking_line_search(x, direction)
    
#     # Update the solution
#     x += step_size * direction
    
#     # Check the stopping criterion
#     if np.linalg.norm(grad_f(x)) < tol:
#         break

# # Print the optimal value and the optimal solution
# print("Optimal value:", f(x))
# print("Optimal solution:", x)








