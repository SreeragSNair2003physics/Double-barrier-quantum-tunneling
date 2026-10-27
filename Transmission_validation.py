
import numpy as np

########################################################################### CONSTANTS ################################################################################

hbar=  1.054571817e-34               # Reduced Planck's constant (Js)
m_e= 9.1093837139e-31               # Mass of an electron (Kg)
eV= 1.602176634e-19                 # one eV; Joules per electronvolt

####################################################################### POTENTIAL BARRIER ############################################################################

V0= 0.8 * eV                        # Barrier height; Joules
W= 2e-9                             # Barrier width; Metres


###########################################################  DEFINING EXACT TRANSMISSION COEFFICIENT ##################################################################


def T_exact(E):                                                                    # Transmission probability for a particle with energy E
    
    E = np.atleast_1d(E).astype(float)
    T = np.zeros_like(E)

    
    T[E <= 0] = 0.0                                                                # Zero or negative energy

    
    threshold = np.isclose(E, V0, rtol=1e-12, atol=0.0)                            # Energy exactly at the barrier threshold
    T[threshold] = 1 / (1 + m_e * V0 * W**2 / (2 * hbar**2))

    
    below = (E > 0) & (E < V0) & ~threshold                                        # Energies below the barrier
    kappa = np.sqrt(2 * m_e * (V0 - E[below])) / hbar
    T[below] = 1 / (1 + V0**2 * np.sinh(kappa * W)**2 / (4 * E[below] * (V0 - E[below])))

    
    above = E > V0                                                                 # Energies above the barrier
    k = np.sqrt(2 * m_e * (E[above] - V0)) / hbar
    T[above] = 1 / (1 + V0**2 * np.sin(k * W)**2 / (4 * E[above] * (E[above] - V0)))

    return T



########################################################### ESTIMATING TRANSMISSION AVERAGED OVER A GAUSSIAN PACKET ###################################################


def T_packet(k0, sigma):
    
    delta_k = 1 / (2 * sigma)                                     # Momentum-space width of the packet

    k_values = np.linspace(max(k0 - 8 * delta_k, 1e6), k0 + 8 * delta_k, 4001)                  # Sample positive wave numbers around the central wave number

    weights = np.exp(-(k_values - k0)**2 / (2 * delta_k**2))      # Gaussian momentum probability distribution

    energies = hbar**2 * k_values**2 / (2 * m_e)                  # Convert wave numbers to energies

    transmissions = T_exact(energies)                             # Average the exact transmission over the distribution
    return np.trapezoid(weights * transmissions, k_values) / np.trapezoid(weights, k_values)



sigma = 0.7e-9
k0 = 4.0e9

packet_transmission = T_packet(k0, sigma)



################################################################## COMPARING NUMERICAL & THEORETICAL RESULTS #########################################################

simulation_right_probability = 0.12134019703908003                                # Finite-time probability found to the right of the barrier in the existing simulation

absolute_difference = abs(packet_transmission - simulation_right_probability)

difference_percentage = (absolute_difference / simulation_right_probability) * 100







########################################################################### SANITY CHECKS ##############################################################################


test_energies = np.array([0.2, 0.5, 0.79, 0.8, 0.81, 1.0]) * eV
test_transmission = T_exact(test_energies)

print("\nTest energies (eV):", test_energies / eV)
print("Transmission:", test_transmission)
print("All values finite:", np.all(np.isfinite(test_transmission)))
print("All values between 0 and 1:",
      np.all((test_transmission >= 0) & (test_transmission <= 1)))
print("Threshold transmission:", T_exact(np.array([V0]))[0])



print("\nGaussian wave-packet test")
print("Central wave number (1/m):", k0)
print("Packet width (m):", sigma)
print("Estimated packet transmission:", packet_transmission)
print("Result is finite:", np.isfinite(packet_transmission))
print("Result is between 0 and 1:", 0 <= packet_transmission <= 1)


print("\nPreliminary theory-simulation comparison")
print(f"Theoretical Gaussian-packet transmission:{packet_transmission:.6f}")
print(f"Finite-time probability to the right:{simulation_right_probability:.6f}")
print(f"Absolute difference:{absolute_difference:.6f}")
print(f"Difference (percentage points):{difference_percentage:.3f}")
print("\nNote: These quantities are not strictly equivalent.")
print("The simulated value is measured at a finite time and may include")
print("effects from wave-packet spreading and incomplete separation.")