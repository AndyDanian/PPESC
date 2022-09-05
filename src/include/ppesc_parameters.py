# LRESC CONSTANTS
M: float = 1000000.0
C: float = 137.0359998
ALPHA: float = 1.0/C
ALPHA2: float = ALPHA*ALPHA
#
ppesc_constants = {
        # NR
        "paranr": M*ALPHA2/2.0,
        # paramagnetic
        "fclap": -M*ALPHA2*ALPHA2*0.5, "sddxx": M*ALPHA2*ALPHA2/2.0,
        "sddyy": M*ALPHA2*ALPHA2/2.0, "sddzz": M*ALPHA2*ALPHA2/2.0,
        "sddxy": M*ALPHA2*ALPHA2/2.0, "sddxz": M*ALPHA2*ALPHA2/2.0,
        "sddyx": M*ALPHA2*ALPHA2/2.0, "sddyz": M*ALPHA2*ALPHA2/2.0,
        "sddzx": M*ALPHA2*ALPHA2/2.0, "sddzy": M*ALPHA2*ALPHA2/2.0,
        "lpsokin": M*ALPHA2*ALPHA2/8.0, "lkinpso": M*ALPHA2*ALPHA2/8.0,
        #"lfcso": M*ALPHA2/2.0, "lsdsoxx": M*ALPHA2/2.0,
        #"lsdsoyy": M*ALPHA2/2.0, "lsdsozz": M*ALPHA2/2.0,
        #"lpsomv": M*ALPHA2/2.0, "lpsodw": M*ALPHA2/2.0,
        # NR
        "dianr": M*ALPHA2,
        # diamagnetic
        "fc": -M*ALPHA2*ALPHA2*0.25, "psooz": -M*ALPHA2*ALPHA2*0.25,
        "sd": -M*ALPHA2*ALPHA2*0.25, "dnske": M*ALPHA2*ALPHA2/12.0,
        "pnstcgop": -M*ALPHA2*ALPHA2*0.25,
        #"a2dw": M*ALPHA2, "a2mv": M*ALPHA2
        }
###### Paramagnetic
# - NR
paranr = [
lambda a: ["pso " + str(1 + 3*a), "angmom x"],
lambda a: ["pso " + str(1 + 3*a), "angmom y"],
lambda a: ["pso " + str(1 + 3*a), "angmom z"],
lambda a: ["pso " + str(2 + 3*a), "angmom x"],
lambda a: ["pso " + str(2 + 3*a), "angmom y"],
lambda a: ["pso " + str(2 + 3*a), "angmom z"],
lambda a: ["pso " + str(3 + 3*a), "angmom x"],
lambda a: ["pso " + str(3 + 3*a), "angmom y"],
lambda a: ["pso " + str(3 + 3*a), "angmom z"],
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
sddxx = [
lambda a: ["laplacian xx", "sd " + str(1 + 3*a) + " x"],
lambda a: ["laplacian yy", "sd " + str(1 + 3*a) + " x"],
lambda a: ["laplacian zz", "sd " + str(1 + 3*a) + " x"],
]
# Sddxy
sddxy = [
lambda a: ["laplacian xx", "sd " + str(1 + 3*a) + " y"],
lambda a: ["laplacian yy", "sd " + str(1 + 3*a) + " y"],
lambda a: ["laplacian zz", "sd " + str(1 + 3*a) + " y"],
]
# Sddxz
sddxz = [
lambda a: ["laplacian xx", "sd " + str(1 + 3*a) + " z"],
lambda a: ["laplacian yy", "sd " + str(1 + 3*a) + " z"],
lambda a: ["laplacian zz", "sd " + str(1 + 3*a) + " z"],
]
# Sddyx
sddyx = [
lambda a: ["laplacian xx", "sd " + str(2 + 3*a) + " x"],
lambda a: ["laplacian yy", "sd " + str(2 + 3*a) + " x"],
lambda a: ["laplacian zz", "sd " + str(2 + 3*a) + " x"],
]
# Sddyy
sddyy = [
lambda a: ["laplacian xx", "sd " + str(2 + 3*a) + " y"],
lambda a: ["laplacian yy", "sd " + str(2 + 3*a) + " y"],
lambda a: ["laplacian zz", "sd " + str(2 + 3*a) + " y"],
]
# Sddyz
sddyz = [
lambda a: ["laplacian xx", "sd " + str(2 + 3*a) + " z"],
lambda a: ["laplacian yy", "sd " + str(2 + 3*a) + " z"],
lambda a: ["laplacian zz", "sd " + str(2 + 3*a) + " z"],
]
# Sddzx
sddzx = [
lambda a: ["laplacian xx", "sd " + str(3 + 3*a) + " x"],
lambda a: ["laplacian yy", "sd " + str(3 + 3*a) + " x"],
lambda a: ["laplacian zz", "sd " + str(3 + 3*a) + " x"],
]
# Sddzy
sddzy = [
lambda a: ["laplacian xx", "sd " + str(3 + 3*a) + " y"],
lambda a: ["laplacian yy", "sd " + str(3 + 3*a) + " y"],
lambda a: ["laplacian zz", "sd " + str(3 + 3*a) + " y"],
]
# Sddzz
sddzz = [
lambda a: ["laplacian xx", "sd " + str(3 + 3*a) + " z"],
lambda a: ["laplacian yy", "sd " + str(3 + 3*a) + " z"],
lambda a: ["laplacian zz", "sd " + str(3 + 3*a) + " z"],
]
triplet_lineal_responses = {"fclap": fclap,
                        "sddxx": sddxx, "sddyy": sddyy, "sddzz": sddzz,
                        "sddxy": sddxy, "sddxz": sddxz, "sddyx": sddyx,
                        "sddyz": sddyz, "sddzx": sddzx, "sddzy": sddzy}

# - Singlet
# L-PsoKin
lpsokin = [
lambda a: ["psoke " + str(1 + a*3), "angmom x"],
lambda a: ["psoke " + str(1 + a*3), "angmom y"],
lambda a: ["psoke " + str(1 + a*3), "angmom z"],
lambda a: ["psoke " + str(2 + a*3), "angmom x"],
lambda a: ["psoke " + str(2 + a*3), "angmom y"],
lambda a: ["psoke " + str(2 + a*3), "angmom z"],
lambda a: ["psoke " + str(3 + a*3), "angmom x"],
lambda a: ["psoke " + str(3 + a*3), "angmom y"],
lambda a: ["psoke " + str(3 + a*3), "angmom z"],
]
# Lkin-Pso
lkinpso = [
lambda a: ["pso " + str(1 + 3*a), "ozke x"],
lambda a: ["pso " + str(1 + 3*a), "ozke y"],
lambda a: ["pso " + str(1 + 3*a), "ozke z"],
lambda a: ["pso " + str(2 + 3*a), "ozke x"],
lambda a: ["pso " + str(2 + 3*a), "ozke y"],
lambda a: ["pso " + str(2 + 3*a), "ozke z"],
lambda a: ["pso " + str(3 + 3*a), "ozke x"],
lambda a: ["pso " + str(3 + 3*a), "ozke y"],
lambda a: ["pso " + str(3 + 3*a), "ozke z"],
]
singlet_lineal_responses = {"lpsokin": lpsokin, "lkinpso": lkinpso}
# # -- Quadratic Responses
# # - Triplet
# # LFcSO
# lfcso = {
# "angmomxfcspinox" : lambda a: ["angmom x", "fc " + str(1 + a), "spinorbit x"],
# "angmomyfcspinoy" : lambda a: ["angmom y", "fc " + str(1 + a), "spinorbit y"],
# "angmomzfcspinoz" : lambda a: ["angmom z", "fc " + str(1 + a), "spinorbit z"],
# }
# # LSdSOxx
# lsdsoxx = {
# "angmomxsd1xspinox" : lambda a: ["angmom x", "sd " + str(1 + 3*a) + " x", "spinorbit x"],
# "angmomysd1yspinoy" : lambda a: ["angmom x", "sd " + str(1 + 3*a) + " y", "spinorbit y"],
# "angmomzsd1zspinoz" : lambda a: ["angmom x", "sd " + str(1 + 3*a) + " z", "spinorbit z"],
# }
# # LSdSOyy
# lsdsoyy = {
# "angmomxsd2xspinox" : lambda a: ["angmom y", "sd " + str(2 + 3*a) + " x", "spinorbit x"],
# "angmomysd2yspinoy" : lambda a: ["angmom y", "sd " + str(2 + 3*a) + " y", "spinorbit y"],
# "angmomzsd2zspinoz" : lambda a: ["angmom y", "sd " + str(2 + 3*a) + " z", "spinorbit z"],
# }
# # LSdSOzz
# lsdsozz = {
# "angmomxsd3xspinox" : lambda a: ["angmom z", "sd " + str(3 + 3*a) + " x", "spinorbit x"],
# "angmomysd3yspinoy" : lambda a: ["angmom z", "sd " + str(3 + 3*a) + " y", "spinorbit y"],
# "angmomzsd3zspinoz" : lambda a: ["angmom z", "sd " + str(3 + 3*a) + " z", "spinorbit z"],
# }
# triplet_quadratic_responses = {"lfcso": lfcso, "lsdsoxx": lsdsoxx, "lsdsoyy": lsdsoyy, "lsdsozz": lsdsozz}
# # - Singlet
# # LpsoMv
# lpsomv = {
# "angmomxpso1massvelo" : lambda a: ["angmom x", "pso " + str(1 + 3*a), "massvelo"],
# "angmomypso2massvelo" : lambda a: ["angmom y", "pso " + str(2 + 3*a), "massvelo"],
# "angmomzpso3massvelo" : lambda a: ["angmom z", "pso " + str(3 + 3*a), "massvelo"],
# }
# # LpsoDw
# lpsodw = {
# "angmomxpso1darwin" : lambda a: ["angmom x", "pso " + str(1 + 3*a), "darwin"],
# "angmomypso2darwin" : lambda a: ["angmom y", "pso " + str(2 + 3*a), "darwin"],
# "angmomzpso3darwin" : lambda a: ["angmom z", "pso " + str(3 + 3*a), "darwin"],
# }
# singlet_quadratic_responses = {"lpsomv" : lpsomv, "lpsodw": lpsodw}

###### Diamagnetic part
## - NR
dianr = [
lambda a: ["nstcgo " + str(1 + 3*a) + " x"],
lambda a: ["nstcgo " + str(1 + 3*a) + " y"],
lambda a: ["nstcgo " + str(1 + 3*a) + " z"],
lambda a: ["nstcgo " + str(2 + 3*a) + " x"],
lambda a: ["nstcgo " + str(2 + 3*a) + " y"],
lambda a: ["nstcgo " + str(2 + 3*a) + " z"],
lambda a: ["nstcgo " + str(3 + 3*a) + " x"],
lambda a: ["nstcgo " + str(3 + 3*a) + " y"],
lambda a: ["nstcgo " + str(3 + 3*a) + " z"],
]
diamagnetic_nr = {"dianr": dianr}
## - Averages
# - FC
fc = [lambda a: ["fc " + str(a)]]
# - SD
sd = [
lambda a : ["sd " + str(1 + 3*a) + " x"],
lambda a : ["sd " + str(1 + 3*a) + " y"],
lambda a : ["sd " + str(1 + 3*a) + " z"],
lambda a : ["sd " + str(2 + 3*a) + " x"],
lambda a : ["sd " + str(2 + 3*a) + " y"],
lambda a : ["sd " + str(2 + 3*a) + " z"],
lambda a : ["sd " + str(3 + 3*a) + " x"],
lambda a : ["sd " + str(3 + 3*a) + " y"],
lambda a : ["sd " + str(3 + 3*a) + " z"],
]
# - PSO-OZ
psooz = [
lambda a : ["psooz " + str(1 + 3*a) + " x"],
lambda a : ["psooz " + str(1 + 3*a) + " y"],
lambda a : ["psooz " + str(1 + 3*a) + " z"],
lambda a : ["psooz " + str(2 + 3*a) + " x"],
lambda a : ["psooz " + str(2 + 3*a) + " y"],
lambda a : ["psooz " + str(2 + 3*a) + " z"],
lambda a : ["psooz " + str(3 + 3*a) + " x"],
lambda a : ["psooz " + str(3 + 3*a) + " y"],
lambda a : ["psooz " + str(3 + 3*a) + " z"],
]
# - DNSKE
dnske =[
lambda a : ["dnske " + str(1 + 3*a) + " x"],
lambda a : ["dnske " + str(1 + 3*a) + " y"],
lambda a : ["dnske " + str(1 + 3*a) + " z"],
lambda a : ["dnske " + str(2 + 3*a) + " x"],
lambda a : ["dnske " + str(2 + 3*a) + " y"],
lambda a : ["dnske " + str(2 + 3*a) + " z"],
lambda a : ["dnske " + str(3 + 3*a) + " x"],
lambda a : ["dnske " + str(3 + 3*a) + " y"],
lambda a : ["dnske " + str(3 + 3*a) + " z"],
]
# - P NSTCGO P
pnstcgop =[
lambda a : ["pnstcgop " + str(1 + 3*a) + " x"],
lambda a : ["pnstcgop " + str(1 + 3*a) + " y"],
lambda a : ["pnstcgop " + str(1 + 3*a) + " z"],
lambda a : ["pnstcgop " + str(2 + 3*a) + " x"],
lambda a : ["pnstcgop " + str(2 + 3*a) + " y"],
lambda a : ["pnstcgop " + str(2 + 3*a) + " z"],
lambda a : ["pnstcgop " + str(3 + 3*a) + " x"],
lambda a : ["pnstcgop " + str(3 + 3*a) + " y"],
lambda a : ["pnstcgop " + str(3 + 3*a) + " z"],
]
dia_averages ={"fc": fc, "sd": sd, "psooz": psooz, "dnske": dnske, "pnstcgop": pnstcgop}
## - Lineal Responses
# - a2dw
# a2dw = {
# "nstcgo1xdw": lambda a : ["nstcgo " + str(1 + 3*a) + " x", "darwin"],
# "nstcgo2ydw": lambda a : ["nstcgo " + str(2 + 3*a) + " y", "darwin"],
# "nstcgo3zdw": lambda a : ["nstcgo " + str(3 + 3*a) + " z", "darwin"],
# }
# # - a2mv
# a2mv = {
# "nstcgo1xmv": lambda a : ["nstcgo " + str(1 + 3*a) + " x", "massvelo"],
# "nstcgo2ymv": lambda a : ["nstcgo " + str(2 + 3*a) + " y", "massvelo"],
# "nstcgo3zmv": lambda a : ["nstcgo " + str(3 + 3*a) + " z", "massvelo"],
# }
# dia_singlet_lineal_responses = {"a2dw": a2dw, "a2mv": a2mv}
# if __name__ == "__main__":
#     for sd in sdkinxx.values():
#         print(sd(1))

##NAME
ppesc_label: dict = {
        "paranr": "<<PSO;ANG>>",
        "fclap": "<<FC;P²>>",
        "sdlap": "<<SD;P²>>",
        "lpsokin": "<<ANG;{PSO,P²}>>",
        "lkinpso": "<<PSO;{ANG,P²}>>",
        "dianr": "<A²>",
        "fc": "<FC>",
        "sd": "<SD>",
        "psooz": "<PSO-OZ>",
        "dnske": "<{A², P²}>",
        "pnstcgop": "<PA²P>"
}