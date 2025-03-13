#plot Z ba matlab fargh dare va male matlab doroste
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def EOM_6DoF_Body(t, states, forces, moments, params):
    m, g, Ixx, Iyy, Izz, Ixz = params
    FAX, FAY, FAZ, FTX, FTY, FTZ = forces
    La, Ma, Na, Lt, Mt, Nt = moments

    U, V, W, P, Q, R, phi, theta, psi, X, Y, Z = states

    Fx = FAX + FTX
    Fy = FAY + FTY
    Fz = FAZ + FTZ

    Udot = (1/m)*(-m*g*np.sin(theta) + Fx) - Q*W + R*V
    Vdot = (1/m)*(m*g*np.cos(theta)*np.sin(phi) + Fy) - R*U + P*W
    Wdot = (1/m)*(m*g*np.cos(theta)*np.cos(phi) + Fz) - P*V + U*Q

    L = La + Lt
    M = Ma + Mt
    N = Na + Nt

    Gamma = Ixx * Izz - Ixz**2

    Pdot = (1/Gamma)*((Izz*L + Ixz*N - (Ixz*(Iyy - Ixx - Izz))*P*Q - (Ixz**2 + Izz*(Izz - Iyy))*Q*R))
    Qdot = (M - (Ixx - Izz)*P*R - Ixz*(P**2 - R**2))/Iyy
    Rdot = (Ixz*L + Ixx*N - (Ixz*(Iyy - Ixx))*P*Q + (Ixz*(Ixx - Iyy + Izz))*Q*R)/Gamma

    phidot = P + Q*np.sin(phi)*np.tan(theta) + R*np.cos(phi)*np.tan(theta)
    thetadot = Q*np.cos(phi) - R*np.sin(phi)
    psidot = (Q*np.sin(phi) + R*np.cos(phi))/np.cos(theta)

    Cb2i = np.array([
        [np.cos(theta)*np.cos(psi), np.cos(theta)*np.sin(psi), -np.sin(theta)],
        [np.sin(phi)*np.sin(theta)*np.cos(psi)-np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(theta)*np.sin(psi)+np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(theta)],
        [np.cos(phi)*np.sin(theta)*np.cos(psi)+np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.sin(psi)-np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(theta)]
    ])

    XYZdot = Cb2i @ np.array([U, V, W])
    Xdot, Ydot, Zdot = XYZdot

    return np.array([Udot, Vdot, Wdot, Pdot, Qdot, Rdot, phidot, thetadot, psidot, Xdot, Ydot, Zdot])

# Load data
data = np.loadtxt('data.txt')
t = data[:, 0]
dt = t[1] - t[0]

# Aircraft parameters
Weight = 564000
g = 32.2
m = Weight / g
Ixx, Iyy, Izz, Ixz = 13700000, 30500000, 49700000, 830000
params = (m, g, Ixx, Iyy, Izz, Ixz)

# Initial conditions
Altitude = 1000
U0, V0, W0 = 262.6411, 0, 22.2674
P0, Q0, R0 = 0, 0, 0
phi0, psi0 = 0, 0
theta0 = np.deg2rad(1.81)
X0, Y0, Z0 = 0, 0, -Altitude

IC = np.array([U0, V0, W0, P0, Q0, R0, phi0, theta0, psi0, X0, Y0, Z0])
states = np.zeros((len(t), len(IC)))
states[0, :] = IC

# Euler simulation
for k in range(1, len(t)):
    forces = data[k, 1:7]
    moments = data[k, 7:13]
    states[k, :] = states[k-1, :] + dt * EOM_6DoF_Body(
        t[k-1], states[k-1, :], forces, moments, params
    )

U, V, W, P, Q, R, phi, theta, psi, X, Y, Z = states.T

alpha = np.arctan2(W, U)
beta = np.arcsin(V / np.sqrt(U**2 + V**2 + W**2))

with PdfPages('All_Plots.pdf') as pdf:
    plot_params = [
        (t, np.rad2deg(phi), 'Phi (deg)'),
        (t, np.rad2deg(theta), 'Theta (deg)'),
        (t, np.rad2deg(psi), 'Psi (deg)'),
        (t, np.rad2deg(P), 'P (deg/s)'),
        (t, np.rad2deg(Q), 'Q (deg/s)'),
        (t, np.rad2deg(R), 'R (deg/s)'),
        (t, np.rad2deg(alpha), 'Alpha (deg)'),
        (t, np.rad2deg(beta), 'Beta (deg)'),
        (t, U, 'U (ft/s)'),
        (t, V, 'V (ft/s)'),
        (t, W, 'W (ft/s)')

    ]

with PdfPages('All_Plots.pdf') as pdf:
    for x, y, label in plot_params:
        plt.figure(figsize=(10, 4))
        plt.plot(x, y)
        plt.xlabel('Time (sec)')
        plt.ylabel(label)
        plt.grid()
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    fig, axs = plt.subplots(3, 1, figsize=(10, 9))
    axs[0].plot(t, X); axs[0].set(ylabel='X (ft)'); axs[0].grid()
    axs[1].plot(t, Y); axs[1].set(ylabel='Y (ft)'); axs[1].grid()
    axs[2].plot(t, Z); axs[2].set(ylabel='Z (ft)', xlabel='Time (sec)'); axs[2].grid()
    plt.tight_layout()
    pdf.savefig(); plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(X, Y, Z); ax.set(xlabel='X (ft)', ylabel='Y (ft)', zlabel='Z (ft)')
    plt.tight_layout()
    pdf.savefig(); plt.close()

