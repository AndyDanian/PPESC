# LRESC CONSTANTS
M: float = 1000000.0
C: float = 137.0359998
ALPHA: float = 1.0 / C
ALPHA2: float = ALPHA * ALPHA
#
lresc_scale: dict[str, float] = {
    # NR
    "paranr": M * ALPHA2 / 2.0,
    # paramagnetic
    "fclap": -0.5 * M * ALPHA2 * ALPHA2 / 3.0,
    "sddxx": -M
    * ALPHA2
    * ALPHA2
    / 4.0
    / 4.0,  # one 1./4. is because old lresc didn't eliminate this constant in the calculate the dptovl in the dalton
    "sddyy": -M * ALPHA2 * ALPHA2 / 4.0 / 4.0,
    "sddzz": -M * ALPHA2 * ALPHA2 / 4.0 / 4.0,
    "sddxy": -M * ALPHA2 * ALPHA2 / 4.0 / 4.0,
    "sddxz": -M * ALPHA2 * ALPHA2 / 4.0 / 4.0,
    "sddyx": -M * ALPHA2 * ALPHA2 / 4.0 / 4.0,
    "sddyz": -M * ALPHA2 * ALPHA2 / 4.0 / 4.0,
    "sddzx": -M * ALPHA2 * ALPHA2 / 4.0 / 4.0,
    "sddzy": -M * ALPHA2 * ALPHA2 / 4.0 / 4.0,
    "fcbso": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsoxx": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsoyy": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsozz": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsoxy": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsoxz": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsoyx": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsoyz": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsozx": -M * ALPHA2 * ALPHA2 / 4.0,
    "sdbsozy": -M * ALPHA2 * ALPHA2 / 4.0,
    "lpsokin": M * ALPHA2 * ALPHA2 / 2.0,
    "lkinpso": M * ALPHA2 * ALPHA2,
    "lfcso": -M * ALPHA2 / 4.0,
    "lsdsoxx": -M * ALPHA2 / 4.0,
    "lsdsoyy": -M * ALPHA2 / 4.0,
    "lsdsozz": -M * ALPHA2 / 4.0,
    "lsdsoxy": -M * ALPHA2 / 4.0,
    "lsdsoxz": -M * ALPHA2 / 4.0,
    "lsdsoyx": -M * ALPHA2 / 4.0,
    "lsdsoyz": -M * ALPHA2 / 4.0,
    "lsdsozx": -M * ALPHA2 / 4.0,
    "lsdsozy": -M * ALPHA2 / 4.0,
    "lpsomv": -M * ALPHA2 / 2.0,
    "lpsodw": -M * ALPHA2 / 2.0,
    # NR
    "dianr": M * ALPHA2,
    # diamagnetic
    "fc": -7.0 * M * ALPHA2 * ALPHA2 / 16.0,
    "psooz": -M * ALPHA2 * ALPHA2 * 0.5,
    "dnske": 10.0 * M * ALPHA2 * ALPHA2 * ALPHA,
    "a2dw": -M * ALPHA2,
    "a2mv": -M * ALPHA2,
}
#
lresc: dict[str, float] = {
    # NR
    "paranr": M * ALPHA2 / 2.0,
    # paramagnetic
    "fclap": -M * ALPHA2 * ALPHA2 * (-0.5),
    "sddxx": M * ALPHA2 * ALPHA2 / 2.0,
    "sdkinyy": M * ALPHA2 * ALPHA2 / 2.0,
    "sddzz": M * ALPHA2 * ALPHA2 / 2.0,
    "sdkinxy": M * ALPHA2 * ALPHA2 / 2.0,
    "sddxz": M * ALPHA2 * ALPHA2 / 2.0,
    "sdkinyx": M * ALPHA2 * ALPHA2 / 2.0,
    "sddyz": M * ALPHA2 * ALPHA2 / 2.0,
    "sdkinzx": M * ALPHA2 * ALPHA2 / 2.0,
    "sddzy": M * ALPHA2 * ALPHA2 / 2.0,
    "fcbso": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsoxx": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsoyy": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsozz": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsoxy": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsoxz": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsoyx": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsoyz": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsozx": 2.0 * M * ALPHA2 * ALPHA2,
    "sdbsozy": 2.0 * M * ALPHA2 * ALPHA2,
    "lpsokin": M * ALPHA2 * ALPHA2 / 8.0,
    "lkinpso": M * ALPHA2 * ALPHA2 / 8.0,
    "lfcso": M * ALPHA2 / 2.0,
    "lsdsoxx": M * ALPHA2 / 2.0,
    "lsdsoyy": M * ALPHA2 / 2.0,
    "lsdsozz": M * ALPHA2 / 2.0,
    "lsdsoxy": M * ALPHA2 / 2.0,
    "lsdsoxz": M * ALPHA2 / 2.0,
    "lsdsoyx": M * ALPHA2 / 2.0,
    "lsdsoyz": M * ALPHA2 / 2.0,
    "lsdsoxz": M * ALPHA2 / 2.0,
    "lsdsozy": M * ALPHA2 / 2.0,
    "lpsomv": M * ALPHA2 / 2.0,
    "lpsodw": M * ALPHA2 / 2.0,
    # NR
    "dianr": M * ALPHA2,
    # diamagnetic
    "fc": -M * ALPHA2 * ALPHA2 / 8.0,
    "psooz": -M * ALPHA2 * ALPHA2 * 0.5,
    "dnske": 2.0 * M * ALPHA2 * ALPHA2 / 3.0,
    "a2dw": M * ALPHA2,
    "a2mv": M * ALPHA2,
}
#
name_order_responses: list[str] = [
    "paranr",
    "fclap",
    "sdlap",
    "fcbso",
    "sdbso",
    "lpsokin",
    "lkinpso",
    "lfcso",
    "lsdso",
    "lpsomv",
    "lpsodw",
    "a2dw",
    "a2mv",
]
#
lresc_constants: dict[str, dict[str, float]] = {
    "lresc_scale": lresc_scale,
    "lresc": lresc,
}
###### Paramagnetic
# - NR
paranr = [
    lambda a: ["pso " + str(1 + 3 * a), "angmom x"],
    lambda a: ["pso " + str(1 + 3 * a), "angmom y"],
    lambda a: ["pso " + str(1 + 3 * a), "angmom z"],
    lambda a: ["pso " + str(2 + 3 * a), "angmom x"],
    lambda a: ["pso " + str(2 + 3 * a), "angmom y"],
    lambda a: ["pso " + str(2 + 3 * a), "angmom z"],
    lambda a: ["pso " + str(3 + 3 * a), "angmom x"],
    lambda a: ["pso " + str(3 + 3 * a), "angmom y"],
    lambda a: ["pso " + str(3 + 3 * a), "angmom z"],
]
paramagnetic_nr = {"paranr": paranr}
# -- Lineal Response
# - Triplet
# FcKin
fclap = [
    lambda a: ["laplacian xx", "fc " + str(1 + a)],
    lambda a: ["laplacian yy", "fc " + str(1 + a)],
    lambda a: ["laplacian zz", "fc " + str(1 + a)],
]
# Sddxx
## In old LRESC used DPTOVL, which according of DALTON
## source is -1/4*ddi, !falta las cruzadas
sddxx = [
    lambda a: ["laplacian xx", "sd " + str(1 + 3 * a) + " x"],
    lambda a: ["laplacian yy", "sd " + str(1 + 3 * a) + " x"],
    lambda a: ["laplacian zz", "sd " + str(1 + 3 * a) + " x"],
]
# Sddxy
sddxy = [
    lambda a: ["laplacian xx", "sd " + str(1 + 3 * a) + " y"],
    lambda a: ["laplacian yy", "sd " + str(1 + 3 * a) + " y"],
    lambda a: ["laplacian zz", "sd " + str(1 + 3 * a) + " y"],
]
# Sddxz
sddxz = [
    lambda a: ["laplacian xx", "sd " + str(1 + 3 * a) + " z"],
    lambda a: ["laplacian yy", "sd " + str(1 + 3 * a) + " z"],
    lambda a: ["laplacian zz", "sd " + str(1 + 3 * a) + " z"],
]
# Sddyx
sddyx = [
    lambda a: ["laplacian xx", "sd " + str(2 + 3 * a) + " x"],
    lambda a: ["laplacian yy", "sd " + str(2 + 3 * a) + " x"],
    lambda a: ["laplacian zz", "sd " + str(2 + 3 * a) + " x"],
]
# Sddyy
sddyy = [
    lambda a: ["laplacian xx", "sd " + str(2 + 3 * a) + " y"],
    lambda a: ["laplacian yy", "sd " + str(2 + 3 * a) + " y"],
    lambda a: ["laplacian zz", "sd " + str(2 + 3 * a) + " y"],
]
# Sddyz
sddyz = [
    lambda a: ["laplacian xx", "sd " + str(2 + 3 * a) + " z"],
    lambda a: ["laplacian yy", "sd " + str(2 + 3 * a) + " z"],
    lambda a: ["laplacian zz", "sd " + str(2 + 3 * a) + " z"],
]
# Sddzx
sddzx = [
    lambda a: ["laplacian xx", "sd " + str(3 + 3 * a) + " x"],
    lambda a: ["laplacian yy", "sd " + str(3 + 3 * a) + " x"],
    lambda a: ["laplacian zz", "sd " + str(3 + 3 * a) + " x"],
]
# Sddzy
sddzy = [
    lambda a: ["laplacian xx", "sd " + str(3 + 3 * a) + " y"],
    lambda a: ["laplacian yy", "sd " + str(3 + 3 * a) + " y"],
    lambda a: ["laplacian zz", "sd " + str(3 + 3 * a) + " y"],
]
# Sddzz
sddzz = [
    lambda a: ["laplacian xx", "sd " + str(3 + 3 * a) + " z"],
    lambda a: ["laplacian yy", "sd " + str(3 + 3 * a) + " z"],
    lambda a: ["laplacian zz", "sd " + str(3 + 3 * a) + " z"],
]
# FcBso
fcbso = [
    lambda a: ["fc " + str(1 + a), "sofiel xx"],
    lambda a: ["fc " + str(1 + a), "sofiel xy"],
    lambda a: ["fc " + str(1 + a), "sofiel xz"],
    lambda a: ["fc " + str(1 + a), "sofiel yx"],
    lambda a: ["fc " + str(1 + a), "sofiel yy"],
    lambda a: ["fc " + str(1 + a), "sofiel yz"],
    lambda a: ["fc " + str(1 + a), "sofiel zx"],
    lambda a: ["fc " + str(1 + a), "sofiel zy"],
    lambda a: ["fc " + str(1 + a), "sofiel zz"],
]
# SdBsoxx
sdbsoxx = [
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel xx"],
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel xy"],
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel xz"],
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel yx"],
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel yy"],
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel yz"],
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel zx"],
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel zy"],
    lambda a: ["sd " + str(1 + 3 * a) + " x", "sofiel zz"],
]
# SdBsoxy
sdbsoxy = [
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel xx"],
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel xy"],
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel xz"],
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel yx"],
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel yy"],
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel yz"],
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel zx"],
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel zy"],
    lambda a: ["sd " + str(1 + 3 * a) + " y", "sofiel zz"],
]
# SdBsoxz
sdbsoxz = [
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel xx"],
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel xy"],
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel xz"],
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel yx"],
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel yy"],
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel yz"],
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel zx"],
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel zy"],
    lambda a: ["sd " + str(1 + 3 * a) + " z", "sofiel zz"],
]
# SdBsoyx
sdbsoyx = [
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel xx"],
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel xy"],
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel xz"],
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel yx"],
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel yy"],
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel yz"],
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel zx"],
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel zy"],
    lambda a: ["sd " + str(2 + 3 * a) + " x", "sofiel zz"],
]
# SdBsoyy
sdbsoyy = [
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel xy"],
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel yy"],
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel zy"],
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel xy"],
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel yy"],
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel zy"],
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel xy"],
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel yy"],
    lambda a: ["sd " + str(2 + 3 * a) + " y", "sofiel zy"],
]
# SdBsoyz
sdbsoyz = [
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel xx"],
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel xy"],
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel xz"],
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel yx"],
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel yy"],
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel yz"],
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel zx"],
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel zy"],
    lambda a: ["sd " + str(2 + 3 * a) + " z", "sofiel zz"],
]
# SdBsozx
sdbsozx = [
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel xx"],
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel xy"],
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel xz"],
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel yx"],
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel yy"],
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel yz"],
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel zx"],
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel zy"],
    lambda a: ["sd " + str(3 + 3 * a) + " x", "sofiel zz"],
]
# SdBsozy
sdbsozy = [
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel xx"],
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel xy"],
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel xz"],
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel yx"],
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel yy"],
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel yz"],
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel zx"],
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel zy"],
    lambda a: ["sd " + str(3 + 3 * a) + " y", "sofiel zz"],
]
# SdBsozz
sdbsozz = [
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel xz"],
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel yz"],
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel zz"],
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel xz"],
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel yz"],
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel zz"],
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel xz"],
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel yz"],
    lambda a: ["sd " + str(3 + 3 * a) + " z", "sofiel zz"],
]
triplet_lineal_responses = {
    "fclap": fclap,
    "sddxx": sddxx,
    "sddyy": sddyy,
    "sddzz": sddzz,
    "sddxy": sddxy,
    "sddxz": sddxz,
    "sddyx": sddyx,
    "sddyz": sddyz,
    "sddzx": sddzx,
    "sddzy": sddzy,
    "fcbso": fcbso,
    "sdbsoxx": sdbsoxx,
    "sdbsoyy": sdbsoyy,
    "sdbsozz": sdbsozz,
    "sdbsoxy": sdbsoxy,
    "sdbsoxz": sdbsoxz,
    "sdbsoyx": sdbsoyx,
    "sdbsoyz": sdbsoyz,
    "sdbsozx": sdbsozx,
    "sdbsozy": sdbsozy,
}
# - Singlet
# L-PsoKin
lpsokin = [
    lambda a: ["psoke " + str(1 + a * 3), "angmom x"],
    lambda a: ["psoke " + str(1 + a * 3), "angmom y"],
    lambda a: ["psoke " + str(1 + a * 3), "angmom z"],
    lambda a: ["psoke " + str(2 + a * 3), "angmom x"],
    lambda a: ["psoke " + str(2 + a * 3), "angmom y"],
    lambda a: ["psoke " + str(2 + a * 3), "angmom z"],
    lambda a: ["psoke " + str(3 + a * 3), "angmom x"],
    lambda a: ["psoke " + str(3 + a * 3), "angmom y"],
    lambda a: ["psoke " + str(3 + a * 3), "angmom z"],
]
# Lkin-Pso
lkinpso = [
    lambda a: ["pso " + str(1 + 3 * a), "ozke x"],
    lambda a: ["pso " + str(1 + 3 * a), "ozke y"],
    lambda a: ["pso " + str(1 + 3 * a), "ozke z"],
    lambda a: ["pso " + str(2 + 3 * a), "ozke x"],
    lambda a: ["pso " + str(2 + 3 * a), "ozke y"],
    lambda a: ["pso " + str(2 + 3 * a), "ozke z"],
    lambda a: ["pso " + str(3 + 3 * a), "ozke x"],
    lambda a: ["pso " + str(3 + 3 * a), "ozke y"],
    lambda a: ["pso " + str(3 + 3 * a), "ozke z"],
]
singlet_lineal_responses = {"lpsokin": lpsokin, "lkinpso": lkinpso}
# -- Quadratic Responses
# - Triplet
# LFcSO
lfcso = [
    lambda a: ["angmom x", "fc " + str(1 + a), "spinorbit x"],
    lambda a: ["angmom x", "fc " + str(1 + a), "spinorbit y"],
    lambda a: ["angmom x", "fc " + str(1 + a), "spinorbit z"],
    lambda a: ["angmom y", "fc " + str(1 + a), "spinorbit x"],
    lambda a: ["angmom y", "fc " + str(1 + a), "spinorbit y"],
    lambda a: ["angmom y", "fc " + str(1 + a), "spinorbit z"],
    lambda a: ["angmom z", "fc " + str(1 + a), "spinorbit x"],
    lambda a: ["angmom z", "fc " + str(1 + a), "spinorbit y"],
    lambda a: ["angmom z", "fc " + str(1 + a), "spinorbit z"],
]
# LSdSOxx
lsdsoxx = [
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " x", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " x", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " x", "spinorbit z"],
]
# LSdSOxy
lsdsoxy = [
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " y", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " y", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " y", "spinorbit z"],
]
# LSdSOxz
lsdsoxz = [
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(1 + 3 * a) + " z", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(1 + 3 * a) + " z", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(1 + 3 * a) + " z", "spinorbit z"],
]
# LSdSOyx
lsdsoyx = [
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " x", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " x", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " x", "spinorbit z"],
]
# LSdSOyy
lsdsoyy = [
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " y", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " y", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " y", "spinorbit z"],
]
# LSdSOyz
lsdsoyz = [
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(2 + 3 * a) + " z", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(2 + 3 * a) + " z", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(2 + 3 * a) + " z", "spinorbit z"],
]
# LSdSOzx
lsdsozx = [
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " x", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " x", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " x", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " x", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " x", "spinorbit z"],
]
# LSdSOzy
lsdsozy = [
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " y", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " y", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " y", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " y", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " y", "spinorbit z"],
]
# LSdSOzz
lsdsozz = [
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom x", "sd " + str(3 + 3 * a) + " z", "spinorbit z"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom y", "sd " + str(3 + 3 * a) + " z", "spinorbit z"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " z", "spinorbit x"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " z", "spinorbit y"],
    lambda a: ["angmom z", "sd " + str(3 + 3 * a) + " z", "spinorbit z"],
]
triplet_quadratic_responses = {
    "lfcso": lfcso,
    "lsdsoxx": lsdsoxx,
    "lsdsoyy": lsdsoyy,
    "lsdsozz": lsdsozz,
    "lsdsoxy": lsdsoxy,
    "lsdsoxz": lsdsoxz,
    "lsdsoyx": lsdsoyx,
    "lsdsoyz": lsdsoyz,
    "lsdsozx": lsdsozx,
    "lsdsozy": lsdsozy,
}
# - Singlet
# LpsoMv
lpsomv = [
    lambda a: ["angmom x", "pso " + str(1 + 3 * a), "massvelo"],
    lambda a: ["angmom x", "pso " + str(2 + 3 * a), "massvelo"],
    lambda a: ["angmom x", "pso " + str(3 + 3 * a), "massvelo"],
    lambda a: ["angmom y", "pso " + str(1 + 3 * a), "massvelo"],
    lambda a: ["angmom y", "pso " + str(2 + 3 * a), "massvelo"],
    lambda a: ["angmom y", "pso " + str(3 + 3 * a), "massvelo"],
    lambda a: ["angmom z", "pso " + str(1 + 3 * a), "massvelo"],
    lambda a: ["angmom z", "pso " + str(2 + 3 * a), "massvelo"],
    lambda a: ["angmom z", "pso " + str(3 + 3 * a), "massvelo"],
]
# LpsoDw
lpsodw = [
    lambda a: ["angmom x", "pso " + str(1 + 3 * a), "darwin"],
    lambda a: ["angmom x", "pso " + str(2 + 3 * a), "darwin"],
    lambda a: ["angmom x", "pso " + str(3 + 3 * a), "darwin"],
    lambda a: ["angmom y", "pso " + str(1 + 3 * a), "darwin"],
    lambda a: ["angmom y", "pso " + str(2 + 3 * a), "darwin"],
    lambda a: ["angmom y", "pso " + str(3 + 3 * a), "darwin"],
    lambda a: ["angmom z", "pso " + str(1 + 3 * a), "darwin"],
    lambda a: ["angmom z", "pso " + str(2 + 3 * a), "darwin"],
    lambda a: ["angmom z", "pso " + str(3 + 3 * a), "darwin"],
]
singlet_quadratic_responses = {"lpsomv": lpsomv, "lpsodw": lpsodw}

###### Diamagnetic part
## - NR
dianr = [
    lambda a: ["nstcgo " + str(1 + 3 * a) + " x"],
    lambda a: ["nstcgo " + str(1 + 3 * a) + " y"],
    lambda a: ["nstcgo " + str(1 + 3 * a) + " z"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " x"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " y"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " z"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " x"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " y"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " z"],
]
diamagnetic_nr = {"dianr": dianr}
## - Averages
# - Fc
fc = [lambda a: ["fc " + str(1 + a)]]
# - PsoOZ
psooz = [
    lambda a: ["psooz " + str(1 + 3 * a) + " x"],
    lambda a: ["psooz " + str(1 + 3 * a) + " y"],
    lambda a: ["psooz " + str(1 + 3 * a) + " z"],
    lambda a: ["psooz " + str(2 + 3 * a) + " x"],
    lambda a: ["psooz " + str(2 + 3 * a) + " y"],
    lambda a: ["psooz " + str(2 + 3 * a) + " z"],
    lambda a: ["psooz " + str(3 + 3 * a) + " x"],
    lambda a: ["psooz " + str(3 + 3 * a) + " y"],
    lambda a: ["psooz " + str(3 + 3 * a) + " z"],
]
# - DNSKE
dnske = [
    lambda a: ["dnske " + str(1 + 3 * a) + " x"],
    lambda a: ["dnske " + str(1 + 3 * a) + " y"],
    lambda a: ["dnske " + str(1 + 3 * a) + " z"],
    lambda a: ["dnske " + str(2 + 3 * a) + " x"],
    lambda a: ["dnske " + str(2 + 3 * a) + " y"],
    lambda a: ["dnske " + str(2 + 3 * a) + " z"],
    lambda a: ["dnske " + str(3 + 3 * a) + " x"],
    lambda a: ["dnske " + str(3 + 3 * a) + " y"],
    lambda a: ["dnske " + str(3 + 3 * a) + " z"],
]
dia_averages = {"fc": fc, "psooz": psooz, "dnske": dnske}
## - Lineal Responses
# - a2dw
a2dw = [
    lambda a: ["nstcgo " + str(1 + 3 * a) + " x", "darwin"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " x", "darwin"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " x", "darwin"],
    lambda a: ["nstcgo " + str(1 + 3 * a) + " y", "darwin"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " y", "darwin"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " y", "darwin"],
    lambda a: ["nstcgo " + str(1 + 3 * a) + " z", "darwin"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " z", "darwin"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " z", "darwin"],
]
# - a2mv
a2mv = [
    lambda a: ["nstcgo " + str(1 + 3 * a) + " x", "massvelo"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " x", "massvelo"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " x", "massvelo"],
    lambda a: ["nstcgo " + str(1 + 3 * a) + " y", "massvelo"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " y", "massvelo"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " y", "massvelo"],
    lambda a: ["nstcgo " + str(1 + 3 * a) + " z", "massvelo"],
    lambda a: ["nstcgo " + str(2 + 3 * a) + " z", "massvelo"],
    lambda a: ["nstcgo " + str(3 + 3 * a) + " z", "massvelo"],
]
dia_singlet_lineal_responses = {"a2dw": a2dw, "a2mv": a2mv}
# if __name__ == "__main__":
#     for sd in sdkinxx.values():
#         print(sd(1))

##NAME
lresc_label: dict = {
    "paranr": "<<PSO;ANG>>",
    "fclap": "<<FC;P²>>",
    "sdlap": "<<SD;P²>>",
    "fcbso": "<<FC;SOFIEL>>",
    "sdbso": "<<SD;SOFIEL>>",
    "lpsokin": "<<ANG;{PSO,P²}>>",
    "lkinpso": "<<PSO;{ANG,P²}>>",
    "lfcso": "<<ANG;FC,SO>>",
    "lsdso": "<<ANG;SD,SO>>",
    "lpsomv": "<<ANG;PSO,Mv>>",
    "lpsodw": "<<ANG;PSO,Dw>>",
    "dianr": "<A²>",
    "fc": "<FC>",
    "sd": "<SD>",
    "psooz": "<PSO-OZ>",
    "dnske": "<{A², P²}>",
    "pnstcgop": "<PA²P>",
    "a2mv": "<<A²;Mv>>",
    "a2dw": "<<A²;Dw>>",
}
