import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.backends.backend_pdf import PdfPages

# Define global variables
m = None
g = 32.2  # Gravity (ft/s²)
Ixx = Iyy = Izz = Ixz = None

FAX = FAY = FAZ = FTX = FTY = FTZ = La = Ma = Na = Lt = Mt = Nt = 0

# Load flight data from file
data = np.loadtxt("data.txt")

# Aircraft mass and inertia properties
Weight = 564000  # lbs
m = Weight / g  # Compute mass
Ixx = 13700000  # slug.ft²
Iyy = 30500000  # slug.ft²
Izz = 49700000  # slug.ft²
Ixz = 830000  # slug.ft²

# Initial flight conditions
Altitude = 1000  # ft
U0, V0, W0 = 262.6411, 0, 22.2674  # Initial velocities (ft/s)
P0, Q0, R0 = 0, 0, 0  # Initial angular rates (rad/s)
phi0, psi0 = 0, 0  # Initial roll & yaw angles (rad)
theta0 = np.radians(1.81)  # Convert pitch angle from degrees to radians
X0, Y0, Z0 = 0, 0, -Altitude  # Initial position (ft)

# State vector initialization
initial_conditions = np.array([U0, V0, W0, P0, Q0, R0, phi0, theta0, psi0, X0, Y0, Z0])

# Extract time data
t = data[:, 0]
dt = t[1] - t[0]

# Define 6DOF Equations of Motion
def EOM_6DoF_Body(t, states):
    global m, g, Ixx, Iyy, Izz, Ixz
    global FAX, FAY, FAZ, FTX, FTY, FTZ, La, Ma, Na, Lt, Mt, Nt

    # Extract state variables
    U, V, W, P, Q, R, phi, theta, psi, X, Y, Z = states

    # Compute total forces
    F_x = FAX + FTX
    F_y = FAY + FTY
    F_z = FAZ + FTZ

    # Linear motion equations
    Udot = (1/m) * (-m*g*np.sin(theta) + F_x) - Q*W + R*V
    Vdot = (1/m) * (m*g*np.cos(theta)*np.sin(phi) + F_y) - R*U + P*W
    Wdot = (1/m) * (m*g*np.cos(theta)*np.cos(phi) + F_z) - P*V + U*Q

    # Angular motion equations
    L, M, N = La + Lt, Ma + Mt, Na + Nt

    Pdot = (1 / (Ixx * Izz - Ixz**2)) * ((Ixx - Iyy + Izz) * Ixz * P * Q + (Iyy * Izz - Izz**2 - Ixz**2) * Q * R + Izz * L + Ixz * N)
    Qdot = (1 / Iyy) * ((Izz - Ixx) * P * R + (R**2 - P**2) * Ixz + M)
    Rdot = (1 / (Ixx * Izz - Ixz**2)) * ((Ixx**2 - Ixx * Iyy + Ixz**2) * P * Q + (-Ixx + Iyy - Izz) * Ixz * Q * R + Ixz * L + Ixx * N)

    # Euler angle rate equations
    phidot = P + Q*np.sin(phi)*np.tan(theta) + R*np.cos(phi)*np.tan(theta)
    thetadot = Q*np.cos(phi) - R*np.sin(phi)
    psidot = (Q*np.sin(phi) + R*np.cos(phi)) / np.cos(theta)

    # Rotation matrix (Body to Inertial)
    Cb2i = np.array([
        [np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
        [np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
        [np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi) - np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]
    ])

    XYZdot = Cb2i @ np.array([U, V, W])

    return np.array([Udot, Vdot, Wdot, Pdot, Qdot, Rdot, phidot, thetadot, psidot, XYZdot[0], XYZdot[1], XYZdot[2]])

# Solve using Euler method
state_history = np.zeros((len(t), len(initial_conditions)))
state_history[0, :] = initial_conditions

for k in range(1, len(t)):
    print(f"t = {t[k]}")
    
    # Update forces from data file
    FAX, FAY, FAZ = data[k, 1], data[k, 2], data[k, 3]
    FTX, FTY, FTZ = data[k, 4], data[k, 5], data[k, 6]
    La, Ma, Na = data[k, 7], data[k, 8], data[k, 9]
    Lt, Mt, Nt = data[k, 10], data[k, 11], data[k, 12]
    
    # Compute next state
    state_history[k, :] = state_history[k-1, :] + dt * EOM_6DoF_Body(t[k-1], state_history[k-1, :])

# Extract states
U, V, W = state_history[:, 0], state_history[:, 1], state_history[:, 2]
P, Q, R = state_history[:, 3], state_history[:, 4], state_history[:, 5]
phi, theta, psi = state_history[:, 6], state_history[:, 7], state_history[:, 8]
X, Y, Z = state_history[:, 9], state_history[:, 10], state_history[:, 11]

# Compute Angle of Attack (α) and Sideslip Angle (β)
alpha = np.arctan(W / U) * 180 / np.pi  # Angle of Attack (deg)
beta = np.arcsin(V / np.sqrt(U**2 + V**2 + W**2)) * 180 / np.pi  # Sideslip Angle (deg)

# Save plots in PDF
with PdfPages("flight_simulation_results.pdf") as pdf:
    variables = [(U, "U Velocity"), (V, "V Velocity"), (W, "W Velocity"),
                 (P, "Roll Rate (P)"), (Q, "Pitch Rate (Q)"), (R, "Yaw Rate (R)"),
                 (phi, "Roll Angle (Phi)"), (theta, "Pitch Angle (Theta)"), (psi, "Yaw Angle (Psi)"),
                 (X, "X Position"), (Y, "Y Position"), (Z, "Z Position"),
                 (alpha, "Angle of Attack (Alpha)"), (beta, "Sideslip Angle (Beta)")]

    for var, title in variables:
        plt.figure()
        plt.plot(t, var)
        plt.title(title)
        plt.xlabel("Time (s)")
        plt.ylabel(title)
        plt.grid()
        pdf.savefig()
        plt.close()

    # 3D Flight Path
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(X, Y, Z)
    ax.set_xlabel("X (ft)")
    ax.set_ylabel("Y (ft)")
    ax.set_zlabel("Z (ft)")
    ax.set_title("3D Flight Path")
    pdf.savefig()
    plt.close()


