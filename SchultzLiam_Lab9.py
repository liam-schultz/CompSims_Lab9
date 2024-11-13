import numpy as np
import matplotlib.pyplot as plt


def solve():
    # * Select numerical parameters (time step, grid spacing, etc.).
    N = 600 #number of grid points
    L = 1200.  # System size (meters)
    h = 2  # Grid spacing for periodic boundary conditions
    v_max = 25.  # Maximum car speed (m/s)
    tau = h/v_max
    nstep = 1500
    coeff = tau / (2 * h)  # Coefficient used by all schemes

    # * Set initial and boundary conditions
    rho_max = 1.0  # Maximum density
    Flow_max = 0.25 * rho_max * v_max  # Maximum Flow
    Flow = np.empty(N)
    cp = np.empty(N)
    cm = np.empty(N)

    """# Initial condition is a square pulse from x = -L/4 to x = 0
    rho = np.zeros(N)
    for i in range(int(N / 4), int(N / 2)):
        rho[i] = rho_max  # Max density in the square pulse

    rho[int(N / 2)] = rho_max / 2  # Try running without this line"""

    #initial condition given in lab assignment description
    rho = np.zeros(N)
    rho[:N//2] = rho_max

    # Use periodic boundary conditions
    ip = np.arange(N) + 1
    ip[N - 1] = 0  # ip = i+1 with periodic b.c.
    im = np.arange(N) - 1
    im[0] = N - 1  # im = i-1 with periodic b.c.

    # * Initialize plotting variables.
    iplot = 1
    xplot = (np.arange(N) - 1 / 2.) * h - L / 2.  # Record x scale for plot
    rplot = np.empty((N, nstep + 1))
    tplot = np.empty(nstep + 1)
    rplot[:, 0] = np.copy(rho)  # Record the initial state
    tplot[0] = 0  # Record the initial time (t=0)

    # * Loop over desired number of steps.
    for istep in range(nstep):

        # * Compute the flow = (Density)*(Velocity)
        Flow[:] = rho[:] * (v_max * (1 - rho[:] / rho_max))

        # * Compute new values of density using Lax method
        rho[:] = .5 * (rho[ip] + rho[im]) - coeff * (Flow[ip] - Flow[im])

        # * Record density for plotting.
        rplot[:, iplot] = np.copy(rho)
        tplot[iplot] = tau * (istep + 1)
        iplot += 1

    return rplot, xplot, tplot

rplot, xplot, tplot = solve()

#Figure 2
snapshots = np.arange(0, 1501, 1500/6, dtype=int)
snapshots[-1] = 1500

for i in snapshots:
    plt.plot(xplot, rplot[:, i])
plt.legend([f"$\\rho(x,t)$ at $t={tplot[i]}$s" for i in snapshots])
plt.show()