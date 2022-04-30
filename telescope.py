import KrakenOS as Kos
import numpy as np
import matplotlib.pyplot as plt

# Ritchey params

D = 400
D1 = 300
D1_in = D1*0.3
D2 = D1*0.3
B = D*1.2
F = 5*D

def Ritchey(D,B,F):
    M = (F-B)/D
    R1 = -2*F/M
    R2 = -2*B/(M-1)
    K1 = -1 - 2/M**3 * B/D
    K2 = -1 - 2/(M-1)**3 * (M*(2*M-1) + B/D)
    return R1, R2, K1, K2

def Hyperbola(F, R):
    return -(F/R - 1)**2

R1, R2, K1, K2 = Ritchey(D, B, F)
# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
# P_Obj.Thickness = 200.0
# P_Obj.Diameter = 100
P_Obj.Glass = "AIR"
P_Obj.Thickness = D*2
P_Obj.Diameter = D1

# ______________________________________#

M1 = Kos.surf()
# M1.Rc = -269.8364227993029
# M1.Thickness = -100.0
# M1.k = -1.056604646839576
# M1.Diameter = 100
# M1.InDiameter = 30.0
M1.Rc = R1
M1.k = K1
M1.Thickness = -D
M1.Diameter = D1
M1.InDiameter = D1_in
M1.Glass = "MIRROR"

# ______________________________________#

M2 = Kos.surf()
# M2.Rc = -97.90599864690856
# M2.Thickness = 97.0
# M2.k = -3.803238575153494
# M2.Diameter = 30.0
M2.Rc = R2
M2.Thickness = B
M2.k = K2
M2.Diameter = D2
M2.Glass = "MIRROR"
M2.AxisMove = 0.0

# # ______________________________________#
#
# L1a = Kos.surf()
# L1a.Rc = 45.27999130048428
# L1a.Thickness = 5.0
# L1a.Glass = "CAF2SCHOTT"
# L1a.Diameter = 16.0
#
# L1b = Kos.surf()
# L1b.Rc = -122.7256512429732
# L1b.Thickness = 112.3 - 97 - 5.0
# L1b.Glass = "AIR"
# L1b.Diameter = 16.0
#
# # ______________________________________#
#
# L2a = Kos.surf()
# L2a.Rc = -30.54229218826073
# L2a.Thickness = 5.5
# L2a.Glass = "CAF2SCHOTT"
# L2a.Diameter = 12.0
#
# L2b = Kos.surf()
# L2b.Rc = 48.11681314249844
# L2b.Thickness = 120 - 112.3 - 5.5 + 1.4
# L2b.Glass = "AIR"
# L2b.Diameter = 12.0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Diameter = 20.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"
# A = [P_Obj, M1, M2, L1a, L1b, L2a, L2b, P_Ima]
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
margin = 10
focus = False
while not focus:

    # ______________________________________#

    v = Kos.BestFocus(X, Y, Z, L, M, N)
    # L2b.Thickness += v
    M2.Thickness += v
    print(f"Image plane shift: {v} [mm]")

    # ______________________________________#

    # A = [P_Obj, M1, M2, L1a, L1b, L2a, L2b, P_Ima]
    A = [P_Obj, M1, M2, P_Ima]
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

    # ______________________________________#

    focus = True

RMS = np.mean(np.abs(Rays.pick(-1)[0:2]), axis=1) - np.mean(Rays.pick(-1)[0:2], axis=1)
X, Y, Z, L, M, N = Rays.pick(-1)

# ______________________________________#

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
ax1.scatter(X, Y, marker = 'x')
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
