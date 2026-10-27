import numpy as np
import matplotlib.pyplot as plt


########################################################################### CONSTANTS ################################################################################

hbar= 1.054571817e-34                      # Reduced Planck constant (Js)
m_e= 9.1093837139e-31                      # Electron mass (Kg)
eV= 1.602176634e-19                        # One electron volt (J)

nm= 1.0e-9                                 # One nanometre in metres (m)



##################################################################### DEFINING SPATIAL GRID ###########################################################################

x_min= -30 * nm                            # Left edge of the grid
x_max= 30 * nm                             # Right edge of grid

N= 5001                                    # Establishing 2000 grid points

x= np.linspace(x_min, x_max, N)            # Creates the evenly spaced grid including endpoints

dx= x[1] - x[0]                            # Grid spacing (nm)


################################################################## ESTABLISHING DOUBLE BARRIER ########################################################################

V0= 0.8 * eV                               # Barrier height (J)
barrier_width= 0.5 * nm                    # Width of the potential barrier (0.5 nm)
well_width= 2.0 * nm                       # Width of the well in between barriers (2 nm)

edges= [-(well_width / 2 + barrier_width), -well_width / 2, well_width / 2, well_width / 2 + barrier_width]           # Defining 4 edges of 2 barriers

potential_values = [0.0, V0, 0.0, V0, 0.0]                                                                            # Defining potential values for regions


################################################################ CREATING POTENTIAL ENERGY ARRAY ######################################################################

V= np.zeros_like(x)                        # Creates an array of zeros with same size as 'x'

V[(x>= edges[0]) & (x< edges[1])]= V0      # Assigns V0= 0.8 eV potential to first barrier
V[(x>= edges[2]) & (x< edges[3])]= V0      # Assigns V0= 0.8 eV potential to second barrier



############################################################ DEFINING INITIAL GUASSIAN WAVE PACKET ####################################################################

x0= -4 * nm                                # Center of wave packet
sigma= 0.7 * nm                            # Width of wave packet
k0= 4e9                                    # Central wave number

psi= np.exp(-(x - x0)**2 / (4 * sigma**2)) * np.exp(1j * k0 * x)                # Initial wavefunction

psi= psi / np.sqrt(np.sum(np.abs(psi)**2) * dx)                                 # Normalizing wavefunction


########################################################################## CREATING DISCRETE HAMILTONIAN ##########################################################################

q= hbar**2 / (2 * m_e * dx**2)

main_H= 2 * q + V
off_H= -q * np.ones(N - 1)




########################################################################## IMPOSING DIRICHLET CONDITIONS ##########################################################################

main_H_int= main_H[1:-1]
off_H_int= off_H[1:-1]




########################################################################## CRANK-NICOLSON CONSTANTS ##########################################################################

dt= 5e-19                         # Time step in seconds
s= dt / (2 * hbar)





########################################################################## BUILDING MATRICES A AND B ##########################################################################

main_A= 1 + 1j * s * main_H_int
off_A= 1j * s * off_H_int

main_B= 1 - 1j * s * main_H_int
off_B= -1j * s * off_H_int

psi_int= psi[1:-1].copy()
psi_initial= psi_int.copy()



########################################################################## RHS MULTIPLIER ##########################################################################

def rhs_multiplier(main, off, psi):
    result= main * psi
    result[1:]+= off * psi[:-1]
    result[:-1]+= off * psi[1:]
    return result




########################################################################## THOMAS ALGORITHM ##########################################################################

def thomas_algo(main, off, rhs):
    n= len(rhs)

    b= main.astype(complex).copy()
    c= off.astype(complex).copy()
    d= rhs.astype(complex).copy()

    for i in range(1, n):
        factor= c[i - 1] / b[i - 1]
        b[i]= b[i] - factor * c[i - 1]
        d[i]= d[i] - factor * d[i - 1]

    solution= np.zeros(n, dtype=complex)
    solution[-1] = d[-1] / b[-1]

    for i in range(n - 2, -1, -1):
        solution[i]= (d[i] - c[i] * solution[i + 1]) / b[i]

    return solution

snapshots = {}

#################################################################### INTERPRETATION OF GRAPH ###########################################################################

def describe_snapshot(x, density, dx, nm):
    # Position of the largest probability-density peak
    peak_index = np.argmax(density)
    peak_position = x[peak_index] / nm
    peak_height = density[peak_index] * nm

    # Probability in the left, central, and right regions
    P_left = np.sum(density[x < edges[0]]) * dx
    P_barriers = np.sum(
        density[
            ((x >= edges[0]) & (x < edges[1])) |
            ((x >= edges[2]) & (x < edges[3]))
        ]
    ) * dx
    P_well = np.sum(
        density[(x >= edges[1]) & (x < edges[2])]
    ) * dx
    P_right = np.sum(density[x >= edges[3]]) * dx

    # Find prominent local peaks
    peak_indices = np.where(
        (density[1:-1] > density[:-2]) &
        (density[1:-1] > density[2:])
    )[0] + 1

    strongest = peak_indices[
        np.argsort(density[peak_indices])[-5:][::-1]
    ]

    print(f"\nPeak position: {peak_position:.3f} nm")
    print(f"Peak density: {peak_height:.5f} per nm")
    print(
        f"Probabilities: left={P_left:.4f}, "
        f"barriers={P_barriers:.4f}, "
        f"well={P_well:.4f}, right={P_right:.4f}"
    )
    print("Five strongest local peaks (position nm, density per nm):")

    for i in strongest:
        print(f"  ({x[i] / nm:.3f}, {density[i] * nm:.5f})")



########################################################################## TIME EVOLUTION ##########################################################################

n_steps= 60000

for step in range(n_steps):

    rhs= rhs_multiplier(main_B, off_B, psi_int)
    psi_int= thomas_algo(main_A, off_A, rhs)

    

    psi_full= np.zeros(N, dtype=complex)
    psi_full[1:-1]= psi_int

    probability_density= np.abs(psi_full)**2
    if (step + 1) % 12000 == 0:
        print(f"\n--- Snapshot at {(step + 1) * dt * 1e15:.2f} fs ---")
        describe_snapshot(x, probability_density, dx, nm)

    if (step + 1) % 12000 == 0:
        snapshots[(step + 1) * dt] = probability_density.copy()

    P_left= np.sum(probability_density[x < edges[0]]) * dx

    P_barriers= np.sum(probability_density[ ((x >= edges[0]) & (x < edges[1])) | ((x >= edges[2]) & (x < edges[3])) ] ) * dx

    P_well= np.sum(probability_density[ (x >= edges[1]) & (x < edges[2])] ) * dx

    P_right= np.sum(probability_density[x >= edges[3]]) * dx

    P_total= np.sum(probability_density) * dx

    x_mean= np.sum(x * probability_density) * dx / P_total

    final_width= np.sqrt(np.sum((x - x_mean)**2 * probability_density) * dx / P_total)


    
if (step + 1) % 12000 == 0:
    print(
        f"Time = {(step + 1) * dt * 1e15:.2f} fs, "
        f"Left = {P_left:.5f}, "
        f"Right = {P_right:.5f}, "
        f"Total = {P_total:.8f}"
    )


# Energy values for the transmission scan
energies_eV = np.linspace(0.05, 0.95, 300)
energies = energies_eV * eV

def transmission_probability(E):
    # Wave numbers in the five regions
    k = np.sqrt(2 * m_e * (E - np.array(potential_values) + 0j)) / hbar

    # Start with the identity transfer matrix
    M = np.eye(2, dtype=complex)

    # Match the wavefunction at each of the four boundaries
    for j, boundary in enumerate(edges):
        k_left = k[j]
        k_right = k[j + 1]

        C_left = np.array([
            [np.exp(1j * k_left * boundary),
             np.exp(-1j * k_left * boundary)],
            [1j * k_left * np.exp(1j * k_left * boundary),
             -1j * k_left * np.exp(-1j * k_left * boundary)]
        ])

        C_right = np.array([
            [np.exp(1j * k_right * boundary),
             np.exp(-1j * k_right * boundary)],
            [1j * k_right * np.exp(1j * k_right * boundary),
             -1j * k_right * np.exp(-1j * k_right * boundary)]
        ])

        M = M @ np.linalg.solve(C_left, C_right)

    # Transmission probability for identical left and right free regions
    transmission_amplitude = 1 / M[0, 0]
    reflection_amplitude = M[1, 0] / M[0, 0]

    T = abs(transmission_amplitude)**2
    R = abs(reflection_amplitude)**2

    return T, R

# Calculate transmission at each energy
# Calculate transmission and reflection at each energy
results = np.array([
    transmission_probability(E) for E in energies
])

T_values = results[:, 0]
R_values = results[:, 1]

print("Minimum transmission:", np.min(T_values))
print("Maximum transmission:", np.max(T_values))
print("All values finite:", np.all(np.isfinite(results)))

print("\nIndependent conservation check:")
for E_eV in [0.1, 0.3, 0.6, 0.9]:
    idx = np.argmin(np.abs(energies_eV - E_eV))
    print(
        f"E = {energies_eV[idx]:.3f} eV, "
        f"T = {T_values[idx]:.6f}, "
        f"R = {R_values[idx]:.6f}, "
        f"T+R = {T_values[idx] + R_values[idx]:.8f}"
    )
        

        













































########################################################################## SANITY CHECKS ##############################################################################

print("\nNumber of grid points:", len(x))                                      # Spatial grid checks
print("Grid spacing (nm):", dx / nm)


print("\nBarrier edges (nm):", np.array(edges) / nm)                           # Potential barrier checks


print("\nInitial total probability:", np.sum(np.abs(psi)**2) * dx)             # Guassian wave packet checks


print(
    f"Step {step + 1}: "
    f"Left={P_left:.4f}, "
    f"Barriers={P_barriers:.4f}, "
    f"Well={P_well:.4f}, "
    f"Right={P_right:.4f}, "
    f"Total={P_total:.8f}, "
    f"Mean position={x_mean/nm:.3f} nm, "
    f"Width={final_width/nm:.3f} nm"
    )


# Find sampled local maxima in the transmission curve
peak_indices = np.where(
    (T_values[1:-1] > T_values[:-2]) &
    (T_values[1:-1] > T_values[2:])
)[0] + 1

print("\nSampled transmission peaks:")
for i in peak_indices:
    print(f"E = {energies_eV[i]:.5f} eV, T = {T_values[i]:.6f}")





# Refine the energy grid around each sampled peak
print("\nRefined peak estimates:")

for i in peak_indices:
    center = energies_eV[i]

    fine_energies_eV = np.linspace(
        max(0.05, center - 0.01),
        min(0.95, center + 0.01),
        2001
    )

    fine_results = np.array([
        transmission_probability(E * eV)
        for E in fine_energies_eV
    ])

    fine_T = fine_results[:, 0]
    j = np.argmax(fine_T)

    print(
        f"E = {fine_energies_eV[j]:.6f} eV, "
        f"T = {fine_T[j]:.8f}"
    )


# Independent transmission and reflection conservation check
print("\nTransmission + reflection check:")

for E_eV in [0.1, 0.3, 0.6, 0.9]:
    E = E_eV * eV
    T, R = transmission_probability(E)

    print(
        f"E = {E_eV:.2f} eV, "
        f"T = {T:.6f}, "
        f"R = {R:.6f}, "
        f"T+R = {T + R:.8f}"
    )



######################################################################### PLOTTING ####################################################################################

plt.plot(x / nm, V / eV)                                                       # Plotting the double barriers      
plt.xlabel("Position (nm)")
plt.ylabel("Potential energy (eV)")
plt.title("Double-barrier potential")
plt.grid()
plt.show()


plt.figure(figsize=(10, 6))

for time, density in snapshots.items():
    plt.plot(
        x / nm,
        density * nm,
        label=f"{time * 1e15:.0f} fs"
    )

plt.plot(x / nm, V / eV, "k--", label="Potential (eV)")
plt.xlabel("Position (nm)")
plt.ylabel("Probability density (nm$^{-1}$)")
plt.title("Wave-packet evolution through a double barrier")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()


# Plot transmission probability versus energy
plt.figure(figsize=(9, 5))
plt.plot(energies_eV, T_values)
plt.xlabel("Incident energy (eV)")
plt.ylabel("Transmission probability")
plt.title("Quantum Tunnelling Through a Double Barrier")
plt.ylim(0, 1.05)
plt.grid(True)
plt.tight_layout()
plt.show()

