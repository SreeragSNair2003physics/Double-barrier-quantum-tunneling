import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

########################################################################### CONSTANTS ################################################################################

hbar=  1.054571817e-34               # Reduced Planck's constant (Js)
m_e= 9.1093837139e-31               # Mass of an electron (Kg)
eV= 1.602176634e-19                 # one eV; Joules per electronvolt

####################################################################### POTENTIAL BARRIER ############################################################################

V0= 0.8 * eV                        # Barrier height; Joules
W= 2e-9                             # Barrier width; Metres




################################################################### INITIAL GUASSIAN WAVEPACKET #####################################################################

sigma= 0.7e-9                        # Initial spatial width (m)
k0= 4.0e9                            # Central wave number (1/m




###########################################################  DEFINING EXACT TRANSMISSION COEFFICIENT ##################################################################

def T_exact(energies):

    energies= np.atleast_1d(energies).astype(float)         # Makes sure it is converted to a Numpy array
    T= np.zeros_like(energies)                              # Array of zeros with same size as energies

    for i, E in enumerate(energies):                        # Examines single list of pairs; (i,E) with one index number

        if E<= 0:
            T[i] = 0.0

        elif np.isclose(E, V0, rtol=1e-12, atol=0.0):
            T[i]= 1 / (1 + m_e * V0 * W**2 / (2 * hbar**2))

        elif E < V0:
            kappa= np.sqrt(2 * m_e * (V0 - E)) / hbar
            T[i]= 1 / (1 + V0**2 * np.sinh(kappa * W)**2 / (4 * E * (V0 - E)))

        else:
            k= np.sqrt(2 * m_e * (E - V0)) / hbar
            T[i]= 1 / (1 + V0**2 * np.sin(k * W)**2 / (4 * E * (E - V0)))


    return T

    



################################################################### CALCULATING ENERGY RANGE ###########################################################################

E0= hbar**2 * k0**2 / (2 * m_e)                             # Kinetic energy assosciated with packet's central wave number (k0)
delta_k= 1 / (2 * sigma)                                    # Wave number spread
delta_E= hbar**2 * k0 * delta_k / m_e                       # Approximate energy spread



############################################################## CALCULATING ENERGY-WEIGHTED TRANSMISSIONS ##########################################################

k_values= np.linspace(0, 5 * k0, 20000)                                      # Creates many possible positive wave numbers
dk= k_values[1] - k_values[0]

momentum_weights= np.exp(-2 * sigma**2 * (k_values - k0)**2)
momentum_weights/= np.sum(momentum_weights)

energy_values= hbar**2 * k_values**2 / (2 * m_e)                             # Converts each wave number into kinetic energy
transmission_values= T_exact(energy_values)                                  # Calculates the transmission probability for every one of those energies

T_packet= np.sum(momentum_weights * transmission_values)



################################################################ DEFINING PACKET TRANSMISSION ######################################################################

def packet_transmission(k0, sigma):
    delta_k = 1 / (2 * sigma)
    k_values = np.linspace(max(k0 - 8 * delta_k, 0), k0 + 8 * delta_k, 2000)
    weights = np.exp(-2 * sigma**2 * (k_values - k0)**2)
    weights /= np.sum(weights)
    energies = hbar**2 * k_values**2 / (2 * m_e)
    return np.sum(weights * T_exact(energies))

k0_scan = np.linspace(1.5e9, 6.5e9, 40)
E_scan_eV = hbar**2 * k0_scan**2 / (2 * m_e) / eV
T_packet_scan = np.array([packet_transmission(k, sigma) for k in k0_scan])




################################################################ CRANK-NICOLSON SIMULATION ##########################################################################

def simulate(k0, sigma, L=30e-9, dx_target=0.02e-9, dt=2e-18, x0=-6e-9):

    
    N = int(round(2 * L / dx_target)) + 1                                       # Grid
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    
    V = np.zeros(N)                                                             # Potential barrier
    V[np.abs(x) <= W / 2] = V0

     
    psi = (1 / (2 * np.pi * sigma**2) ** 0.25) * np.exp(-(x - x0)**2 / (4 * sigma**2)) * np.exp(1j * k0 * x)      # Initial Gaussian packet
    psi[0] = 0
    psi[-1] = 0
    psi_int = psi[1:-1].copy()

    
    q = hbar**2 / (2 * m_e * dx**2)                                             # Hamiltonian and Crank-Nicolson matrices A and B
    main_H = (2 * q + V)[1:-1]
    s = dt / (2 * hbar)

    main_A = 1 + 1j * s * main_H
    off_A = 1j * s * (-q)
    main_B = 1 - 1j * s * main_H
    off_B = -1j * s * (-q)

    
    ab = np.zeros((3, N - 2), dtype=complex)                                    # Matrix A in the layout solve_banded expects
    ab[0, 1:] = off_A
    ab[1, :] = main_A
    ab[2, :-1] = off_A

                             
    v = hbar * k0 / m_e                                                         # Number of time steps: enough for the packet to travel 14 nm past the barrier
    n_steps = int(round((-x0 + W / 2 + 14e-9) / v / dt))

    
    for step in range(n_steps):                                                 # Time evolution
        rhs = main_B * psi_int
        rhs[1:] += off_B * psi_int[:-1]
        rhs[:-1] += off_B * psi_int[1:]
        psi_int = solve_banded((1, 1), ab, rhs)

    
    psi_full = np.zeros(N, dtype=complex)                                       # Probability to the right of the barrier
    psi_full[1:-1] = psi_int
    density = np.abs(psi_full)**2
    total = np.sum(density) * dx
    P_right = np.sum(density[x > 0]) * dx / total

    return P_right, total


################################################################ RUNNING THE SIMULATION AT SEVERAL ENERGIES ##########################################################

sim_k0 = np.linspace(2.0e9, 6.0e9, 8)
sim_E_eV = hbar**2 * sim_k0**2 / (2 * m_e) / eV

sim_results = [simulate(k, sigma) for k in sim_k0]
sim_T = np.array([r[0] for r in sim_results])
sim_norm = np.array([r[1] for r in sim_results])

for E_i, Ts, Ta, nrm in zip(sim_E_eV, sim_T, [packet_transmission(k, sigma) for k in sim_k0], sim_norm):
    print(f"E = {E_i:.3f} eV   simulation = {Ts:.4f}   exact = {Ta:.4f}   total probability = {nrm:.8f}")





############################################################################ PLOTTING ###############################################################################

energies_eV = np.linspace(0.01, 1.5, 1000)                                   # Energy range for the plot
energies_J = energies_eV * eV

T_values = T_exact(energies_J)                                               # Calculating exact transmission probabilities

plt.figure(figsize=(8, 5))
plt.plot(energies_eV, T_values, "--", color="gray", label="Exact T(E), single energy")
plt.plot(E_scan_eV, T_packet_scan, label="Exact T averaged over the packet")
plt.plot(sim_E_eV, sim_T, "o", label="Crank-Nicolson simulation")
plt.axvline(V0 / eV, color="orange", alpha=0.5, label="Barrier height")
plt.xlabel("Central energy of packet (eV)")
plt.ylabel("Transmitted probability")
plt.title("Wave-packet transmission vs analytic transmission")
plt.ylim(0, 1)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("packet_vs_single_energy.png", dpi=300)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(energies_eV, T_values, label="Exact transmission")
plt.axvline(V0 / eV, linestyle="--", label="Barrier height")

plt.xlabel("Particle energy (eV)")
plt.ylabel("Transmission probability")
plt.title("Transmission Through a Rectangular Potential Barrier")
plt.ylim(0, 1)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("transmission_vs_energy.png", dpi=300)
plt.show()



########################################################################### SANITY CHECKS ##############################################################################





energies= np.array([0.2, 0.5, 0.79, 0.81, 1.0]) * eV
transmission= T_exact(energies)

print("\nTransmission values:", transmission)
print("All finite:", np.all(np.isfinite(transmission)))
print("All between 0 and 1:", np.all((transmission >= 0) & (transmission <= 1)))

print("\nAt threshold:", T_exact(np.array([V0])))

print(f"\nCentral energy: {E0 / eV:.4f} eV")
print(f"Approximate energy spread: {delta_E / eV:.4f} eV")
print(f"Barrier height: {V0 / eV:.4f} eV")


print(f"\nEnergy-weighted packet transmission: {T_packet:.6f}")

print("\nSum of weights:", np.sum(momentum_weights))
print("Largest weight:", np.max(momentum_weights))