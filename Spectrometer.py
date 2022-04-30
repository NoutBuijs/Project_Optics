import KrakenOS as Kos
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

# Optics params

D = 200
D1 = 300
D1_in = D1*0.3
D2 = D1*0.3
B = D*1.2
F = 5*D1
F2 = 65
angle = 30

def beta_finder_der2(b, W, D, alpha):
    return W/D * W * ( (np.sin(b)/(-np.sin(b - np.pi/3) - np.sin(b)) -
                        np.cos(b)*(-np.cos(b - np.pi/3) - np.cos(b))/(-np.sin(b - np.pi/3) - np.sin(b))**2))

def beta_finder2(b, W, R, D, alpha):
    return W / D * (W * np.cos(b) / (np.sin(np.pi / 3 - b) - np.sin(b))) - W / R

def beta_finder_der(b, W, D, alpha):
    return W/D * W * (np.cos(b)**2 / (np.sin(alpha) - np.sin(b))**2 - np.sin(b)/(np.sin(alpha) - np.sin(b)))

def beta_finder(b, W, R, D, alpha):
    return W/D * (W*np.cos(b)/(np.sin(alpha) - np.sin(b))) - W/R

def grating(W, R, D, alpha, order = 3, start = 0.1):
    a = np.radians(alpha)
    b = newton(lambda x: beta_finder(x, W, R, D, a), start, lambda x: beta_finder_der(x, W, D, a))%(np.pi*2)
    d = W*order/(np.sin(a) - np.sin(b))
    dbdl = (np.sin(a) - np.sin(b))/(W*np.cos(b))
    return np.degrees(b), d, dbdl

def grating2(W, R, D, alpha, order = 3, start = 0.1):
    a = np.radians(alpha)
    b = newton(lambda x: beta_finder2(x, W, R, D, a), start, lambda x: beta_finder_der2(x, W, D, a))%(np.pi*2)
    d = W*order/(np.sin(a) - np.sin(b))
    dbdl = (np.sin(a) - np.sin(b))/(W*np.cos(b))
    return np.degrees(b), d, dbdl

def Ritchey(D,B,F):
    M = (F-B)/D
    R1 = -2*F/M
    R2 = -2*B/(M-1)
    K1 = -1 - 2/M**3 * B/D
    K2 = -1 - 2/(M-1)**3 * (M*(2*M-1) + B/D)
    return R1, R2, K1, K2

def Hyperbola(F, R):
    return -(F/R - 1)**2

def Parabola(F):
    return -F*2

def thick_lens_finder_der(R, d, n):
    return (n-1)*(-2/R**2 - 2*(n-1)*d/(n*R**3))

def thick_lens_finder(R, d, f, n):
    return (n-1)*(2/R + (n-1)*d/(n*R**2)) - 1/f


def Spherical_lens(f, d, n):
    return newton(lambda x: thick_lens_finder(x, d, f, n), 1, lambda x: thick_lens_finder_der(x, d, n))

R1, R2, K1, K2 = Ritchey(D, B, F)

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Glass = "AIR"
P_Obj.Thickness = D*2
P_Obj.Diameter = D1

# ______________________________________#

M1 = Kos.surf()
M1.Rc = R1
M1.k = K1
M1.Thickness = -D
M1.Diameter = D1
M1.InDiameter = D1_in
M1.Glass = "MIRROR"

# ______________________________________#

M2 = Kos.surf()
M2.Rc = R2
M2.Thickness = B
M2.k = K2
M2.Diameter = D2
M2.Glass = "MIRROR"

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Diameter = 200
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"
A = [P_Obj, M1, M2, P_Ima]

# ______________________________________#

config = Kos.Setup()
Telescope = Kos.system(A, config)
Rays = Kos.raykeeper(Telescope)

# ______________________________________#

W = 10.0
sup = 1
AperType = "STOP"
Pup = Kos.PupilCalc(Telescope, sup, W, AperType)
Pup.Samp = 7
Pup.FieldType = "angle"

# ______________________________________#

Pup.FieldX = 0.0
Pup.FieldY = 0.0
Pup.Ptype = "hexapolar"
xa, ya, za, La, Ma, Na = Pup.Pattern2Field()

# ______________________________________#

for i in range(0, len(xa)):
    pSource_0 = [xa[i], ya[i], za[i]]
    dCos = [La[i], Ma[i], Na[i]]
    Telescope.Trace(pSource_0, dCos, W)
    Rays.push()

X, Y, Z, L, M, N = Rays.pick(-1)

# ______________________________________#

focus = False
while not focus:

    # ______________________________________#

    v = Kos.BestFocus(X, Y, Z, L, M, N)
    M2.Thickness += v+F2
    print(f"Image plane shift: {v} [mm]")

    # ______________________________________#

    M3 = Kos.surf()
    M3.Rc = Parabola(F2)
    M3.Thickness = -100
    M3.k = -1
    M3.Diameter = 50
    M3.Glass = "MIRROR"
    M3.TiltX = angle
    M3.AxisMove = 2.0

    # ______________________________________#

    b, d, dbdl = grating2(W, 2557, 40E3, 51.40292275686964, order=3, start=0.1)

    Dif_Obj = Kos.surf()
    Dif_Obj.Rc = 0.0
    Dif_Obj.Thickness = 70
    Dif_Obj.Glass = "MIRROR"
    Dif_Obj.Diameter = 50
    Dif_Obj.Grating_D = d
    Dif_Obj.Diff_Ord = 3
    Dif_Obj.TiltX = -(60-b)
    Dif_Obj.AxisMove = 1.0
    Dif_Obj.Surface_type = 1
    Dif_Obj.Grating_Angle = 180

    # ______________________________________#

    P_Ima.TiltX = -b
    P_Ima.DespY = np.sin(np.radians(b))*Dif_Obj.Thickness

    # ______________________________________#

    A = [P_Obj, M1, M2, M3, Dif_Obj, P_Ima]
    Telescope = Kos.system(A, config)
    Rays = Kos.raykeeper(Telescope)

    # ______________________________________#

    W = 10.0
    sup = 1
    AperType = "STOP"
    Pup = Kos.PupilCalc(Telescope, sup, W, AperType)
    Pup.Samp = 4
    Pup.FieldType = "angle"

    # ______________________________________#

    Pup.FieldX = 0.0
    Pup.FieldY = 0.0
    Pup.Ptype = "hexapolar"
    xa, ya, za, La, Ma, Na = Pup.Pattern2Field()

    # ______________________________________#

    for i in range(0, len(xa)):
        pSource_0 = [xa[i], ya[i], za[i]]
        dCos = [La[i], Ma[i], Na[i]]
        Telescope.Trace(pSource_0, dCos, W)
        Rays.push()
        Telescope.Trace(pSource_0, dCos, W-.1)
        Rays.push()
        Telescope.Trace(pSource_0, dCos, W+.1)
        Rays.push()
        Telescope.Trace(pSource_0, dCos, 0.45)
        Rays.push()
        Telescope.Trace(pSource_0, dCos, 0.5)
        Rays.push()
        Telescope.Trace(pSource_0, dCos, 0.7)
        Rays.push()

    # ______________________________________#

    focus = True

RMS = np.mean(np.abs(Rays.pick(-1)[0:2]), axis=1) - np.mean(Rays.pick(-1)[0:2], axis=1)
X, Y, Z, L, M, N = Rays.pick(-1)

# ______________________________________#

Kos.display2d(Telescope, Rays, 0)
Kos.display3d(Telescope, Rays, 2)

# ______________________________________#

X_p, Y_p, Z_p, P2V = Kos.Phase(Pup)
EFFL = Telescope.Parax(W)[-6]
NC = 18
A = np.ones(38)
Zcoef, Mat, w_rms = Kos.Zernike_Fitting(X_p, Y_p, Z_p, A)
S = np.exp(-(2*np.pi*w_rms/(W*10**6))**2)
ima = Kos.WavefrontData2Image(Zcoef, 400)
AB = Kos.Seidel(Pup)

for i in range(0, NC):
    print("z", i + 1, "  ", "{0:.6f}".format(float(Zcoef[i])), ":", Mat[i])

# ______________________________________#

fig = plt.figure()

ax1 = plt.subplot2grid((1,2),(0,0), rowspan=1, colspan=1)
ax1.scatter(X-np.mean(X), Y-np.mean(Y), marker = 'x')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Spot Diagram')
ax1.axis('equal')

ax2 = plt.subplot2grid((1,2),(0,1), rowspan=1, colspan=1, projection = "3d")
ax2.plot_trisurf(X_p, Y_p, Z_p, cmap="coolwarm")
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_title('PSF')

# ______________________________________#
