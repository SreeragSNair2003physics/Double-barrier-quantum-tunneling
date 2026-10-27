import numpy as np
import matplotlib.pyplot as plt

###################################################################### DEFINING CONSTANTS  ##########################################################################

hbar= 1.054571817e-34 # Reduced Planck's constant (Js)
m_e= 9.1093837139e-31 # Mass of electron (Kg)
eV= 1.602176634e-19   # Electronvolt; energy equal to kinetic energy of electron when accelerated through potential difference of one volt in a vacuum (J)


###################################################################### CREATING THE GRID  ##########################################################################

x_min= -10e-9         # e-9 means 10^-9, so -10e-9 is -10 nano metres
x_max= 10e-9
N= 2000               # Specifying number of gridpoints; including the endpoints

x= np.linspace(x_min, x_max, N)     # Creates an array of 2000 evenly spaced positions
dx= x[1] - x[0]       # Defines space between each grid point


############################################################### DEFINING GUASSIAN WAVEPACKET ##########################################################################

x0= -4.0e-9           # Initial center of wavepacket; it is at -4 nano metres
sigma= 0.7e-9         # Packet width; smaller=> localized packet ; larger=> broader packet
k0= 4.0e9            # Central wave number; determines packet's central momentum; for a non-relativistic electron, it's central kinetic energy

E_mean = hbar**2 * k0**2 / (2 * m_e)           # Calculates electron's mean kinetic energy

delta_k = 1 / (2 * sigma)                      # Calculates wave number spread of the wavepacket
delta_E = hbar**2 * k0 * delta_k / m_e         # Calculates energy spread of the wavepacket

A= 1 / (2 * np.pi * sigma**2)**0.25     # Normalization constant
psi= A * np.exp(-(x - x0)**2 / (4 * sigma**2)) * np.exp(1j * k0 * x)     # Initial Guassian wavepacket with Guassian envelope and oscillating part


initial_width= np.sqrt( np.sum((x - x0)**2 * np.abs(psi)**2) * dx )      # Checking initial width of wavepacket
print("Initial spatial width (m):", initial_width)




############################################################### CREATING POTENTIAL BARRIER ##########################################################################

V0= 0.8 * eV          # Barrier height  in Joules
barrier_center = 0.0  # Barrier centre at x = 0 metres
barrier_width= 2e-9   # Barrier width; here it is 2 nanometres

left_edge = barrier_center - barrier_width / 2
right_edge = barrier_center + barrier_width / 2

V= np.zeros(N)        # Initially, the entire grid space has zero potential energy
#V[np.abs(x) <= barrier_width / 2] = V0   # Removed for free particle test


############################################################# CREATING DISCRETE HAMILTONIAN ##########################################################################

q= hbar**2 / (2 * m_e * dx**2)

main_H= 2 * q + V     # Main diagonal entries
off_H= -q * np.ones(N - 1)   # Off diagonal entries


############################################################# IMPOSING DIRICHLET CONDITIONS ##########################################################################

main_H_int = main_H[1:-1]    # Removes the fixed zero values at endpoints because of Dirichlet conditions; Only the interior points remain
off_H_int = off_H[1:-1]      # Removes the fixed zero values at endpoints because of Dirichlet conditions; Only the interior points remain


############################################################# CRANK NICOLSON CONSTANTS ###################################################################################

dt= 1e-18             # Time-step in seconds; here it is 1 x 10^-18 seconds
s= dt / (2 * hbar)    # Common factor from Crank Nicolson equation


############################################################# BUILDING MATRIXES A & B ################################################################################

main_A = 1 + 1j * s * main_H_int
off_A = 1j * s * off_H_int

main_B = 1 - 1j * s * main_H_int
off_B = -1j * s * off_H_int

psi_int = psi[1:-1].copy()    # Removes fixed zero values at endpoints because of Dirichlet conditions; Only interior points remain; Copies so psi isn't rewritten


############################################################# RHS MULTIPLIER ALGORITHM ##################################################################################

def rhs_multiplier(main, off, psi):
    result= main * psi
    result[1:] += off * psi[:-1]
    result[:-1] += off * psi[1:]
    return result

##################################################################### RHS CALCULATION ##################################################################################

rhs= rhs_multiplier(main_B, off_B, psi_int)

##################################################################### THOMAS ALGORITHM #################################################################################

def thomas_algo(main, off, rhs):
    n= len(rhs)
    b= main.astype(complex).copy()    # Makes a separate, complex-number copy of the main-diagonal array
    c= off.astype(complex).copy()     # Makes a separate, complex-number copy of the off-diagonal array
    d= rhs.astype(complex).copy()     # Makes a separate, complex-number copy of the rhs array

    for i in range(1, n):
        factor= off[i - 1] / b[i - 1]
        b[i]= b[i] - factor * off[i - 1]
        d[i]= d[i] - factor * d[i - 1]
    
    solution= np.zeros(n, dtype=complex)

    solution[-1]= d[-1] / b[-1]

    for i in range(n - 2, -1, -1):
        solution[i]= (d[i] - off[i] * solution[i + 1]) / b[i]

    return solution



psi_initial = psi_int.copy()          # Saves initial wavefunction for comparison (Used in sanity check)


##################################################################### SOLVING WITH TIME STEPS ##########################################################################

n_steps= 16000        # Time evolution

for step in range(n_steps):
    rhs = rhs_multiplier(main_B, off_B, psi_int)
    psi_next = thomas_algo(main_A, off_A, rhs)
    psi_int = psi_next.copy()



    if (step + 1) % 2000 == 0:                                                # Seperate probability checks at every checkpoint (2k,4k,6k,...)    
        psi_full = np.zeros(N, dtype=complex)
        psi_full[1:-1] = psi_int

        P_left = np.sum(np.abs(psi_full[x < left_edge])**2) * dx
        P_barrier = np.sum(np.abs(psi_full[(x >= left_edge) & (x <= right_edge)])**2) * dx
        P_right = np.sum(np.abs(psi_full[x > right_edge])**2) * dx

        P_total = np.sum(np.abs(psi_int)**2) * dx                              # Total probability checks across the grid at every checkpoint
        print(f"Step {step + 1}: Total probability = {P_total:.10f}")          # Sanity checked here to remove redundancy



psi_full= np.zeros(N, dtype=complex)                                          # Final width check of wavepacket
psi_full[1:-1] = psi_int

probability_density= np.abs(psi_full)**2
P_total= np.sum(probability_density) * dx

x_mean= np.sum(x * probability_density) * dx / P_total

final_width= np.sqrt( np.sum((x - x_mean)**2 * probability_density) * dx / P_total )

print("Final average position (m):", x_mean)
print("Final spatial width (m):", final_width)




















######################################################################### SANITY CHECKS ################################################################################

probability = np.sum(np.abs(psi_int)**2) * dx
change = np.max(np.abs(psi_int - psi_initial))

print("\nChange in wavefunction:", change)                          # If nonzero, wavefunction has evolved
print("Total probability after evolution:", probability)          # Has to be extremely close to 1; Probability conservation
print("Time steps completed:", n_steps)                           # Confirmation of number of time steps
print("Mean kinetic energy (eV):", E_mean / eV)                   # Prints the mean kinetic energy of the electron
print("Wave number spread (1/m):", delta_k)
print("Approximate energy spread (eV):", delta_E / eV)







probability_density = np.abs(psi_int)**2
x_int = x[1:-1]

average_position = np.sum(x_int * probability_density) * dx / probability      # Calculates average position of the wave packet after time evolution

print("\nAverage position after evolution (m):", average_position)
print("Initial packet centre (m):", x0)






velocity = hbar * k0 / m_e                                        # Calculates group velocity of initial wave packet

print("\nExpected group velocity (m/s):", velocity)
print("Expected travel distance (m):", velocity * dt * n_steps)   # Calculates how much distance the wave packet moved from initial position








barrier_left = barrier_center - barrier_width / 2
distance_to_barrier = barrier_left - x0                           # Estimates the distance for the packet centre to reach the barrier
time_to_barrier = distance_to_barrier / velocity                  # Estimates the time for the packet centre to reach the barrier
steps_to_barrier = time_to_barrier / dt                           # Estimate the steps for the packet centre to reach the barrier

print("\nDistance to barrier (m):", distance_to_barrier)
print("Estimated time to barrier (s):", time_to_barrier)
print("Estimated time steps to barrier:", steps_to_barrier)







psi_full = np.zeros(N, dtype=complex)
psi_full[1:-1] = psi_int

left_region = x < barrier_left
barrier_region = (x >= barrier_left) & (x <= right_edge)
right_region = x > right_edge

P_left = np.sum(np.abs(psi_full[left_region])**2) * dx             # Measuring probability in each region of space we've made
P_barrier = np.sum(np.abs(psi_full[barrier_region])**2) * dx
P_right = np.sum(np.abs(psi_full[right_region])**2) * dx

print("\nProbability left of barrier:", P_left)
print("Probability inside barrier:", P_barrier)
print("Probability right of barrier:", P_right)
print("Sum of regional probabilities:", P_left + P_barrier + P_right)





right_boundary_region = x > 8e-9                                   
P_near_right_boundary = ( np.sum(np.abs(psi_full[right_boundary_region])**2) * dx )           # checks the final 2 nm of the grid

print("\nProbability near right boundary:", P_near_right_boundary)











print("\nProbability distribution after time evolution:")
print(f"Left of barrier: {P_left:.4%}")                             # From the if loop inside for loop of time step definition
print(f"Inside barrier: {P_barrier:.4%}")
print(f"Right of barrier: {P_right:.4%}")






total_time = n_steps * dt
distance_to_right_boundary = x_max - right_edge                #calculating total simulated time and the approximate travel time from the barrier to the right boundary
time_to_boundary = distance_to_right_boundary / velocity            

print("\nTotal simulated time (s):", total_time)                                                  
print("Distance from barrier to right boundary (m):", distance_to_right_boundary)
print("Estimated travel time to right boundary (s):", time_to_boundary)








print("\nFinal-state validation:")                                                      # Validation checks
print(f"Total probability: {P_left + P_barrier + P_right:.10f}")
print(f"Probability accounting error: {abs(1.0 - (P_left + P_barrier + P_right)):.3e}")
































































