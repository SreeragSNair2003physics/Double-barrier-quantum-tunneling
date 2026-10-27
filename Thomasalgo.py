# Defining the Thomas algorithm for tridiagonal solving

def thomas_algo(lower, main, upper, rhs):


    n= len(main)     # Establishing grid size with number of main diagonal elements, saves in variable 'n'

    # Temporary copies to use outside function if needed 
    lower = lower.copy()
    main = main.copy()
    upper = upper.copy()
    rhs = rhs.copy()


# Forward elimination

    for i in range(1, n):     

        factor = lower[i - 1] / main[i - 1]  # Common factor to use for variable elimination

        main[i] = main[i] - factor * upper[i - 1] # Changing the main diagonal

        rhs[i] = rhs[i] - factor * rhs[i - 1]  # Update the right-hand side during elimination


# Back substitution     
   
    solution = [0j] * n   # List maker for complex numbers

    solution[n - 1] = rhs[n - 1] / main[n - 1]   # Solving the last grid point

    for i in range(n - 2, -1, -1):  # Countdown loop to solve grid backwards; n-2 bcoz last grid point is already solved by above line; visits indices 0 & 1

        solution[i] = ( rhs[i] - upper[i] * solution[i + 1] ) / main[i]   # Calculating answer for current position

    return solution


