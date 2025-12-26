import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid
import matplotlib.pyplot as plt
from scipy.special import gamma, gammaincc  # gamma and upper incomplete gamma

phi_g = 0.7134 
g0 = 0.0838 
phi_a = 0.749 
a0 = 0.396 
E0 = 0.0001
myu0 = 0.035 
alpha = 0.1

x = np.exp(np.linspace(np.log(1), np.log(1e6), 50))   
z = np.exp(np.linspace(np.log(1), np.log(1e6), 200))                           

# --- Functions ---
def g(x,Ez): 
    return g0 * x**phi_g

def mu(x,Ez):
    return myu0

def betax(x,Ez):
    return (alpha/(1-alpha)) * g(x,Ez) * (1 - Ez[0])  

# Weight func for env calc
def Ac(x): 
    return a0 * x**phi_a

def E_func(u, x, z):
    Ez = np.zeros_like(z)
    for j, zj in enumerate(z):
        a = x >= zj   # integrate only over diameters >= z, an array a is formed, and it contains all values over the zj.
        if np.any(a):
            Ez[j] = simpson(Ac(x[a]) * u[a], x[a])
    return E0 * Ez
    
# --- Solver ---
def int_1_over_g(x, Ez):
    return cumulative_trapezoid(1.0/g(x, Ez), x, initial=0)

def solve_u():
    u = x**(-phi_g)
    u = u / simpson(u, x)
    Ez = E_func(u, x, z) 
    cnt = 0
    for it in range(2000):
        Ez = E_func(u, x, z)          # E(z) for all heights
        cnt = cnt + 1
        I = int_1_over_g(x, Ez)
        u_safe = np.clip(u, 1e-12, 1e3)  # avoid under/overflow
        Bx = betax(x, Ez)
        B = simpson(Bx * u_safe, x)

        u_new = (B / g(x, Ez)) * np.exp(-mu(x,Ez) * I)
        

        if np.max(np.abs(u_new - u)) < 1e-6:
            print(f"Converged after {cnt} iterations")
            return u_new, Ez
        u = u_new
   
    return u, Ez


u_final, Ez_final = solve_u()

def Ueq(x):
    # Constants
    a0 = 0.396
    phiA = 0.749
    g0 = 0.0838
    phiG = 0.7134
    m0 = 1
    mort = 0.035
    mu0 = mort * m0 / g0
    alpha = 0.10

    temp = mu0 / (1 - phiG)
    
    # Coverage term (note: gammaincc = upper incomplete gamma normalized by gamma(a))
    coverage = 1 - (1 - alpha) / alpha * mu0 / (
        (temp ** (phiG / (phiG - 1))) *
        np.exp(mu0 / (1 - phiG)) *
        gammaincc(phiG / (1 - phiG) + 1, mu0 / (1 - phiG)) * gamma(phiG / (1 - phiG) + 1)
    )

    Neq = coverage / a0 / (
        (temp ** (phiA / (phiG - 1))) *
        np.exp(temp) *
        gammaincc(phiA / (1 - phiG) + 1, temp) * gamma(phiA / (1 - phiG) + 1)
    )

    n0 = Neq * mort / g0

    # Vectorized output
    return n0 * (x / m0) * (-phiG) * np.exp(mu0 / (1 - phiG) * (1 - (x / m0) * (1 - phiG))) * 10000


plt.figure()
plt.plot(np.log10(x), np.log10(u_final), 'o')
plt.plot(np.log10(x), np.log10(Ueq(x)))
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Converged u(x)")

plt.figure()
Ez_plot = np.clip(Ez_final, 1e-12, None)
plt.plot(np.log10(z), np.log10(Ez_plot))
plt.xlabel("Height z")
plt.ylabel("E(z)")
plt.title("E(z) vs height")