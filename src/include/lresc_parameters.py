# LRESC CONSTANTS
M: float = 1000000.0
C: float = 137.0359998
ALPHA: float = 1.0/C
ALPHA2: float = ALPHA*ALPHA
lresc_scale = {"fckin": M*ALPHA2*ALPHA2, "sdkinxx": -M*ALPHA2*ALPHA2/4.0/4.0, #one 1./4. is because old lresc didn't eliminate this constant in the calculate the dptovl in the dalton
        "sdkinyy": -M*ALPHA2*ALPHA2/4.0/4.0, "sdkinzz": -M*ALPHA2*ALPHA2/4.0/4.0,
        "fcbso": -M*ALPHA2*ALPHA2/4.0, "sdbsoxx": -M*ALPHA2*ALPHA2/4.0,
        "sdbsoyy": -M*ALPHA2*ALPHA2/4.0, "sdbsozz": -M*ALPHA2*ALPHA2/4.0,
        "lpsokin": M*ALPHA2*ALPHA2/2.0, "lkinpso": M*ALPHA2*ALPHA2,
        "lfcso": -M*ALPHA2/4.0, "lsdsoxx": -M*ALPHA2/4.0,
        "lsdsoyy": -M*ALPHA2/4.0, "lsdsozz": -M*ALPHA2/4.0,
        "lpsomv": -M*ALPHA2/2.0, "lpsodw": -M*ALPHA2/2.0}
lresc = {"fckin": -M*ALPHA2*ALPHA2, "sdkinxx": M*ALPHA2*ALPHA2/2.0,
        "sdkinyy": M*ALPHA2*ALPHA2/2.0, "sdkinzz": M*ALPHA2*ALPHA2/2.0,
        "fcbso": 2.0*M*ALPHA2*ALPHA2, "sdbsoxx": 2.0*M*ALPHA2*ALPHA2,
        "sdbsoyy": 2.0*M*ALPHA2*ALPHA2, "sdbsozz": 2.0*M*ALPHA2*ALPHA2,
        "lpsokin": M*ALPHA2*ALPHA2/8.0, "lkinpso": M*ALPHA2*ALPHA2/8.0,
        "lfcso": M*ALPHA2/2.0, "lsdsoxx": M*ALPHA2/2.0,
        "lsdsoyy": M*ALPHA2/2.0, "lsdsozz": M*ALPHA2/2.0,
        "lpsomv": M*ALPHA2/2.0, "lpsodw": M*ALPHA2/2.0}
lresc_constants = {"lresc_scale": lresc_scale, "lresc": lresc}
###### Paramagnetic Corrections
# -- Lineal Response
# - Triplet
# FcKin
fckin = {"fckin" : lambda a: ["kinetic", "fc " + str(1 + a)]}
# SdKinxx
sdkinxx = {
"sd1xddxx" : lambda a: ["laplacian xx", "sd " + str(1 + 3*a) + " x"], #In old LRESC used DPTOVL, which according of DALTON
"sd1xddxy" : lambda a: ["laplacian xy", "sd " + str(1 + 3*a) + " x"], #source is -1/4*ddi, !falta las cruzadas
"sd1xddxz" : lambda a: ["laplacian xz", "sd " + str(1 + 3*a) + " x"],
"sd1yddxx" : lambda a: ["laplacian xx", "sd " + str(1 + 3*a) + " y"],
"sd1yddxy" : lambda a: ["laplacian xy", "sd " + str(1 + 3*a) + " y"],
"sd1yddxz" : lambda a: ["laplacian xz", "sd " + str(1 + 3*a) + " y"],
"sd1zddxx" : lambda a: ["laplacian xx", "sd " + str(1 + 3*a) + " z"],
"sd1zddxy" : lambda a: ["laplacian xy", "sd " + str(1 + 3*a) + " z"],
"sd1zddxz" : lambda a: ["laplacian xz", "sd " + str(1 + 3*a) + " z"],
}
# SdKinyy
sdkinyy = {
"sd2xddyx" : lambda a: ["laplacian xy", "sd " + str(2 + 3*a) + " x"], # dxy
"sd2xddyy" : lambda a: ["laplacian yy", "sd " + str(2 + 3*a) + " x"], # dyy
"sd2xddyz" : lambda a: ["laplacian yz", "sd " + str(2 + 3*a) + " x"], # dyz
"sd2yddyx" : lambda a: ["laplacian xy", "sd " + str(2 + 3*a) + " y"],
"sd2yddyy" : lambda a: ["laplacian yy", "sd " + str(2 + 3*a) + " y"],
"sd2yddyz" : lambda a: ["laplacian yz", "sd " + str(2 + 3*a) + " y"],
"sd2zddyx" : lambda a: ["laplacian xy", "sd " + str(2 + 3*a) + " z"],
"sd2zddyy" : lambda a: ["laplacian yy", "sd " + str(2 + 3*a) + " z"],
"sd2zddyz" : lambda a: ["laplacian yz", "sd " + str(2 + 3*a) + " z"],
}
# SdKinzz
sdkinzz = {
"sd3xddzx" : lambda a: ["laplacian xz", "sd " + str(3 + 3*a) + " x"], # dxz
"sd3xddzy" : lambda a: ["laplacian yz", "sd " + str(3 + 3*a) + " x"], # dyz
"sd3xddzz" : lambda a: ["laplacian zz", "sd " + str(3 + 3*a) + " x"], # dzz
"sd3yddzx" : lambda a: ["laplacian xz", "sd " + str(3 + 3*a) + " y"],
"sd3yddzy" : lambda a: ["laplacian yz", "sd " + str(3 + 3*a) + " y"],
"sd3yddzz" : lambda a: ["laplacian zz", "sd " + str(3 + 3*a) + " y"],
"sd3zddzx" : lambda a: ["laplacian xz", "sd " + str(3 + 3*a) + " z"],
"sd3zddzy" : lambda a: ["laplacian yz", "sd " + str(3 + 3*a) + " z"],
"sd3zddzz" : lambda a: ["laplacian zz", "sd " + str(3 + 3*a) + " z"],
}
# FcBso
fcbso = {
"fcsofielxx" : lambda a: ["fc " + str(1 + a), "sofiel xx"],
"fcsofielyy" : lambda a: ["fc " + str(1 + a), "sofiel yy"],
"fcsofielzz" : lambda a: ["fc " + str(1 + a), "sofiel zz"],
}
# SdBsoxx
sdbsoxx = {
"sd1xsofielxx" : lambda a: ["sd " + str(1 + 3*a) + " x", "sofiel xx"],
"sd1xsofielxy" : lambda a: ["sd " + str(1 + 3*a) + " x", "sofiel xy"],
"sd1xsofielxz" : lambda a: ["sd " + str(1 + 3*a) + " x", "sofiel xz"],
"sd1ysofielxx" : lambda a: ["sd " + str(1 + 3*a) + " y", "sofiel xx"],
"sd1ysofielxy" : lambda a: ["sd " + str(1 + 3*a) + " y", "sofiel xy"],
"sd1ysofielxz" : lambda a: ["sd " + str(1 + 3*a) + " y", "sofiel xz"],
"sd1zsofielxx" : lambda a: ["sd " + str(1 + 3*a) + " z", "sofiel xx"],
"sd1zsofielxy" : lambda a: ["sd " + str(1 + 3*a) + " z", "sofiel xy"],
"sd1zsofielxz" : lambda a: ["sd " + str(1 + 3*a) + " z", "sofiel xz"],
}
# SdBsoyy
sdbsoyy = {
"sd2xsofielyx" : lambda a: ["sd " + str(2 + 3*a) + " x", "sofiel yx"],
"sd2xsofielyy" : lambda a: ["sd " + str(2 + 3*a) + " x", "sofiel yy"],
"sd2xsofielyz" : lambda a: ["sd " + str(2 + 3*a) + " x", "sofiel yz"],
"sd2ysofielyx" : lambda a: ["sd " + str(2 + 3*a) + " y", "sofiel yx"],
"sd2ysofielyy" : lambda a: ["sd " + str(2 + 3*a) + " y", "sofiel yy"],
"sd2ysofielyz" : lambda a: ["sd " + str(2 + 3*a) + " y", "sofiel yz"],
"sd2zsofielyx" : lambda a: ["sd " + str(2 + 3*a) + " z", "sofiel yx"],
"sd2zsofielyy" : lambda a: ["sd " + str(2 + 3*a) + " z", "sofiel yy"],
"sd2zsofielyz" : lambda a: ["sd " + str(2 + 3*a) + " z", "sofiel yz"],
}
# SdBsozz
sdbsozz = {
"sd3xsofielzx" : lambda a: ["sd " + str(3 + 3*a) + " x", "sofiel zx"],
"sd3xsofielzy" : lambda a: ["sd " + str(3 + 3*a) + " x", "sofiel zy"],
"sd3xsofielzz" : lambda a: ["sd " + str(3 + 3*a) + " x", "sofiel zz"],
"sd3ysofielzx" : lambda a: ["sd " + str(3 + 3*a) + " y", "sofiel zx"],
"sd3ysofielzy" : lambda a: ["sd " + str(3 + 3*a) + " y", "sofiel zy"],
"sd3ysofielzz" : lambda a: ["sd " + str(3 + 3*a) + " y", "sofiel zz"],
"sd3zsofielzx" : lambda a: ["sd " + str(3 + 3*a) + " z", "sofiel zx"],
"sd3zsofielzy" : lambda a: ["sd " + str(3 + 3*a) + " z", "sofiel zy"],
"sd3zsofielzz" : lambda a: ["sd " + str(3 + 3*a) + " z", "sofiel zz"],
}
triplet_lineal_responses = {"fckin": fckin, "sdkinxx": sdkinxx, "sdkinyy": sdkinyy, "sdkinzz": sdkinzz,
                        "fcbso": fcbso, "sdbsoxx":sdbsoxx, "sdbsoyy": sdbsoyy, "sdbsozz": sdbsozz}
# - Singlet
# L-PsoKin
lpsokin = {
"angmomxpsoke1" : lambda a: ["angmom x", "psoke " + str(1 + a*3)],
"angmomypsoke2" : lambda a: ["angmom y", "psoke " + str(2 + a*3)],
"angmomzpsoke3" : lambda a: ["angmom z", "psoke " + str(3 + a*3)],
}
# Lkin-Pso
lkinpso = {
"psoxozke1" : lambda a: ["pso " + str(1 + 3*a), "ozke x"],
"psoyozke2" : lambda a: ["pso " + str(2 + 3*a), "ozke y"],
"psozozke3" : lambda a: ["pso " + str(3 + 3*a), "ozke z"],
}
singlet_lineal_responses = {"lpsokin": lpsokin, "lkinpso": lkinpso}
# -- Quadratic Responses
# - Triplet
# LFcSO
lfcso = {
"angmomxfcspinox" : lambda a: ["angmom x", "fc " + str(1 + a), "spinorbit x"],
"angmomyfcspinoy" : lambda a: ["angmom y", "fc " + str(1 + a), "spinorbit y"],
"angmomzfcspinoz" : lambda a: ["angmom z", "fc " + str(1 + a), "spinorbit z"],
}
# LSdSOxx
lsdsoxx = {
"angmomxsd1xspinox" : lambda a: ["angmom x", "sd " + str(1 + 3*a) + " x", "spinorbit x"],
"angmomysd1yspinoy" : lambda a: ["angmom x", "sd " + str(1 + 3*a) + " y", "spinorbit y"],
"angmomzsd1zspinoz" : lambda a: ["angmom x", "sd " + str(1 + 3*a) + " z", "spinorbit z"],
}
# LSdSOyy
lsdsoyy = {
"angmomxsd2xspinox" : lambda a: ["angmom y", "sd " + str(2 + 3*a) + " x", "spinorbit x"],
"angmomysd2yspinoy" : lambda a: ["angmom y", "sd " + str(2 + 3*a) + " y", "spinorbit y"],
"angmomzsd2zspinoz" : lambda a: ["angmom y", "sd " + str(2 + 3*a) + " z", "spinorbit z"],
}
# LSdSOzz
lsdsozz = {
"angmomxsd3xspinox" : lambda a: ["angmom z", "sd " + str(3 + 3*a) + " x", "spinorbit x"],
"angmomysd3yspinoy" : lambda a: ["angmom z", "sd " + str(3 + 3*a) + " y", "spinorbit y"],
"angmomzsd3zspinoz" : lambda a: ["angmom z", "sd " + str(3 + 3*a) + " z", "spinorbit z"],
}
triplet_quadratic_responses = {"lfcso": lfcso, "lsdsoxx": lsdsoxx, "lsdsoyy": lsdsoyy, "lsdsozz": lsdsozz}
# - Singlet
# LpsoMv
lpsomv = {
"angmomxpso1massvelo" : lambda a: ["angmom x", "pso " + str(1 + 3*a), "massvelo"],
"angmomypso2massvelo" : lambda a: ["angmom y", "pso " + str(2 + 3*a), "massvelo"],
"angmomzpso3massvelo" : lambda a: ["angmom z", "pso " + str(3 + 3*a), "massvelo"],
}
# LpsoDw
lpsodw = {
"angmomxpso1darwin" : lambda a: ["angmom x", "pso " + str(1 + 3*a), "darwin"],
"angmomypso2darwin" : lambda a: ["angmom y", "pso " + str(2 + 3*a), "darwin"],
"angmomzpso3darwin" : lambda a: ["angmom z", "pso " + str(3 + 3*a), "darwin"],
}
singlet_quadratic_responses = {"lpsomv" : lpsomv, "lpsodw": lpsodw}

# if __name__ == "__main__":
#     for sd in sdkinxx.values():
#         print(sd(1))