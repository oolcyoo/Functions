import sympy as sp

def find_extrema_and_saddle_points(f, a, b):

    # Compute the first and second derivatives of f
    f_prime = sp.diff(f, x)
    f_double_prime = sp.diff(f_prime, x)

    # Find critical points by solving f'(x) = 0
    critical_points = sp.solve(f_prime, x)

    # Initialize variables to keep track of global min and max
    global_min = (None, float('inf'))
    global_max = (None, float('-inf'))

    # Analyze the critical points
    local_minima = []
    local_maxima = []
    saddle_points = []

    for point in critical_points:
        second_derivative = f_double_prime.subs(x, point)
        second_derivative_value = sp.N(second_derivative)  # Convert to numerical value

        # Check if the second_derivative is real (imaginary part is zero)
        if sp.im(second_derivative_value) == 0:
            if second_derivative_value > 0:
                local_minima.append(point)
                value = f.subs(x, point)
                if value < global_min[1]:
                    global_min = (point, value)
                if value > global_max[1]:
                    global_max = (point, value)

            elif second_derivative_value < 0:
                local_maxima.append(point)
                value = f.subs(x, point)
                if value < global_min[1]:
                    global_min = (point, value)
                if value > global_max[1]:
                    global_max = (point, value)

            else:
                n_derivative = 2
                f_n = f_double_prime
                derive_value = sp.N(f_n.subs(x, point))
                while derive_value == 0:
                    f_n = sp.diff(f_n, x)
                    derive_value = sp.N(f_n.subs(x, point))
                    n_derivative += 1
                if n_derivative % 2 != 0:
                    saddle_points.append(point)
                elif derive_value > 0:
                    local_minima.append(point)
                    value = f.subs(x, point)
                    if value < global_min[1]:
                        global_min = (point, value)
                    if value > global_max[1]:
                        global_max = (point, value)
                else:
                    local_maxima.append(point)
                    value = f.subs(x, point)
                    if value < global_min[1]:
                        global_min = (point, value)
                    if value > global_max[1]:
                        global_max = (point, value)
    
                

    # Evaluate the function at the endpoints
    if f_prime.subs(x, a) > 0:
        local_minima.append(a)

    elif f_prime.subs(x, a) < 0:
        local_maxima.append(a)


    if f_prime.subs(x, b) < 0:
        local_minima.append(b)

    elif f_prime.subs(x, b) > 0:
        local_maxima.append(b)

    # # Update the global min and max if needed
    # if f_a < global_min[1]:
    #     global_min = (a, f_a)
    #     local_minima.append(a)
    # if f_a > global_max[1]:
    #     global_max = (a, f_a)
    #     local_maxima.append(a)

    # if f_b < global_min[1]:
    #     global_min = (b, f_b)
    #     local_minima.append(b)
    # if f_b > global_max[1]:
    #     global_max = (b, f_b)
    #     local_minima.append(b)

    return local_minima, local_maxima, global_min, global_max, saddle_points

# Example usage
x = sp.Symbol('x')
f = x**6 - 2*x**3 - 1
a = -2
b = 2

local_minima, local_maxima, global_min, global_max, saddle_points = find_extrema_and_saddle_points(f, a, b)

print("Local Minima:", local_minima)
print("Local Maxima:", local_maxima)
print("Global Min:", global_min)
print("Global Max:", global_max)
print("Saddle Points:", saddle_points)