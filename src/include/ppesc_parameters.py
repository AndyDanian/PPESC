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
paranr = {
"lxpso1": lambda a: ["angmom x", "pso " + str(1 + 3*a)],
"lypso2": lambda a: ["angmom y", "pso " + str(2 + 3*a)],
"lzpso3": lambda a: ["angmom z", "pso " + str(3 + 3*a)],
}
paramagnetic_nr = {"paranr": paranr}
# -- Lineal Response
# - Triplet
# FcKin
fclap = {
"fcdxx" : lambda a: ["laplacian xx", "fc " + str(1 + a)],
"fcdyy" : lambda a: ["laplacian yy", "fc " + str(1 + a)],
"fcdzz" : lambda a: ["laplacian zz", "fc " + str(1 + a)],
        }
# Sddxx
sddxx = {
"sd1xdxx" : lambda a: ["laplacian xx", "sd " + str(1 + 3*a) + " x"], #In old LRESC used DPTOVL, which according of DALTON
"sd1xdyy" : lambda a: ["laplacian yy", "sd " + str(1 + 3*a) + " x"], #source is -1/4*ddi, !falta las cruzadas
"sd1xdzz" : lambda a: ["laplacian zz", "sd " + str(1 + 3*a) + " x"],
}
# Sddyy
sddyy = {
"sd2ydyx" : lambda a: ["laplacian xx", "sd " + str(2 + 3*a) + " y"],
"sd2ydyy" : lambda a: ["laplacian yy", "sd " + str(2 + 3*a) + " y"],
"sd2ydyz" : lambda a: ["laplacian zz", "sd " + str(2 + 3*a) + " y"],
}
# Sddzz
sddzz = {
"sd3zdzx" : lambda a: ["laplacian xx", "sd " + str(3 + 3*a) + " z"],
"sd3zdzy" : lambda a: ["laplacian yy", "sd " + str(3 + 3*a) + " z"],
"sd3zdzz" : lambda a: ["laplacian zz", "sd " + str(3 + 3*a) + " z"],
}
triplet_lineal_responses = {"fclap": fclap, "sddxx": sddxx, "sddyy": sddyy, "sddzz": sddzz}
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
dianr = {
"nstcgo1x": lambda a: ["nstcgo " + str(1 + 3*a) + " x"],
"nstcgo2y": lambda a: ["nstcgo " + str(2 + 3*a) + " y"],
"nstcgo3z": lambda a: ["nstcgo " + str(3 + 3*a) + " z"],
}
diamagnetic_nr = {"dianr": dianr}
## - Averages
# - FC
fc = {"fc": lambda a: ["fc " + str(a)]}
# - SD
sd = {
"sd1x": lambda a : ["sd " + str(1 + 3*a) + " x"],
"sd2y": lambda a : ["sd " + str(2 + 3*a) + " y"],
"sd3z": lambda a : ["sd " + str(3 + 3*a) + " z"],
}
# - PSO-OZ
psooz = {
"psooz1x": lambda a : ["psooz " + str(1 + 3*a) + " x"],
"psooz2y": lambda a : ["psooz " + str(2 + 3*a) + " y"],
"psooz3z": lambda a : ["psooz " + str(3 + 3*a) + " z"],
}
# - DNSKE
dnske ={
"dnske1x": lambda a : ["dnske " + str(1 + 3*a) + " x"],
"dnske2y": lambda a : ["dnske " + str(2 + 3*a) + " y"],
"dnske3z": lambda a : ["dnske " + str(3 + 3*a) + " z"],
}
# - P NSTCGO P
pnstcgop ={
"pnstcgop1x": lambda a : ["pnstcgop " + str(1 + 3*a) + " x"],
"pnstcgop2y": lambda a : ["pnstcgop " + str(2 + 3*a) + " y"],
"pnstcgop3z": lambda a : ["pnstcgop " + str(3 + 3*a) + " z"],
}
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