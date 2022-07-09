###### Paramagnetic Corrections
# -- Lineal Response
# - Triplet
# FcKin
fckin = ["fc " + str(1 + a), "kinetic"]
# SdKinxx
sd1xddxx = ["sd " + str(1 + 3*a) + " x", "laplacian 1"] #In old LRESC used DPTOVL, which according of DALTON
sd1xddxy = ["sd " + str(1 + 3*a) + " x", "laplacian 4"] #source is -1/4*ddi, !falta las cruzadas
sd1xddxz = ["sd " + str(1 + 3*a) + " x", "laplacian 5"]
sd1yddxx = ["sd " + str(1 + 3*a) + " y", "laplacian 1"]
sd1yddxy = ["sd " + str(1 + 3*a) + " y", "laplacian 4"]
sd1yddxz = ["sd " + str(1 + 3*a) + " y", "laplacian 5"]
sd1zddxx = ["sd " + str(1 + 3*a) + " z", "laplacian 1"]
sd1zddxy = ["sd " + str(1 + 3*a) + " z", "laplacian 4"]
sd1zddxz = ["sd " + str(1 + 3*a) + " z", "laplacian 5"]
# SdKinyy
sd2xddyx = ["sd " + str(2 + 3*a) + " x", "laplacian 4"] # dxy
sd2xddyy = ["sd " + str(2 + 3*a) + " x", "laplacian 2"] # dyy
sd2xddyz = ["sd " + str(2 + 3*a) + " x", "laplacian 6"] # dyz
sd2yddyx = ["sd " + str(2 + 3*a) + " y", "laplacian 4"]
sd2yddyy = ["sd " + str(2 + 3*a) + " y", "laplacian 2"]
sd2yddyz = ["sd " + str(2 + 3*a) + " y", "laplacian 6"]
sd2zddyx = ["sd " + str(2 + 3*a) + " z", "laplacian 4"]
sd2zddyy = ["sd " + str(2 + 3*a) + " z", "laplacian 2"]
sd2zddyz = ["sd " + str(2 + 3*a) + " z", "laplacian 6"]
# SdKinzz
sd3xddzx = ["sd " + str(3 + 3*a) + " x", "laplacian 5"] # dxz
sd3xddzy = ["sd " + str(3 + 3*a) + " x", "laplacian 6"] # dyz
sd3xddzz = ["sd " + str(3 + 3*a) + " x", "laplacian 3"] # dzz
sd3yddzx = ["sd " + str(3 + 3*a) + " y", "laplacian 5"]
sd3yddzy = ["sd " + str(3 + 3*a) + " y", "laplacian 6"]
sd3yddzz = ["sd " + str(3 + 3*a) + " y", "laplacian 3"]
sd3zddzx = ["sd " + str(3 + 3*a) + " z", "laplacian 5"]
sd3zddzy = ["sd " + str(3 + 3*a) + " z", "laplacian 6"]
sd3zddzz = ["sd " + str(3 + 3*a) + " z", "laplacian 3"]
# FcBso
fcsofielxx = ["fc " + str(1 + a), "sofiel xx"]
fcsofielyy = ["fc " + str(1 + a), "sofiel yy"]
fcsofielzz = ["fc " + str(1 + a), "sofiel zz"]
# SdBsoxx
sd1xsofielxx = ["sd " + str(1 + 3*a) + " x", "sofiel xx"]
sd1xsofielxy = ["sd " + str(1 + 3*a) + " x", "sofiel xy"]
sd1xsofielxz = ["sd " + str(1 + 3*a) + " x", "sofiel xz"]
sd1ysofielxx = ["sd " + str(1 + 3*a) + " y", "sofiel xx"]
sd1ysofielxy = ["sd " + str(1 + 3*a) + " y", "sofiel xy"]
sd1ysofielxz = ["sd " + str(1 + 3*a) + " y", "sofiel xz"]
sd1zsofielxx = ["sd " + str(1 + 3*a) + " z", "sofiel xx"]
sd1zsofielxy = ["sd " + str(1 + 3*a) + " z", "sofiel xy"]
sd1zsofielxz = ["sd " + str(1 + 3*a) + " z", "sofiel xz"]
# SdBsoyy
sd2xsofielyx = ["sd " + str(2 + 3*a) + " x", "sofiel yx"]
sd2xsofielyy = ["sd " + str(2 + 3*a) + " x", "sofiel yy"]
sd2xsofielyz = ["sd " + str(2 + 3*a) + " x", "sofiel yz"]
sd2ysofielyx = ["sd " + str(2 + 3*a) + " y", "sofiel yx"]
sd2ysofielyy = ["sd " + str(2 + 3*a) + " y", "sofiel yy"]
sd2ysofielyz = ["sd " + str(2 + 3*a) + " y", "sofiel yz"]
sd2zsofielyx = ["sd " + str(2 + 3*a) + " z", "sofiel yx"]
sd2zsofielyy = ["sd " + str(2 + 3*a) + " z", "sofiel yy"]
sd2zsofielyz = ["sd " + str(2 + 3*a) + " z", "sofiel yz"]
# SdBsozz
sd3xsofielzx = ["sd " + str(3 + 3*a) + " x", "sofiel zx"]
sd3xsofielzy = ["sd " + str(3 + 3*a) + " x", "sofiel zy"]
sd3xsofielzz = ["sd " + str(3 + 3*a) + " x", "sofiel zz"]
sd3ysofielzx = ["sd " + str(3 + 3*a) + " y", "sofiel zx"]
sd3ysofielzy = ["sd " + str(3 + 3*a) + " y", "sofiel zy"]
sd3ysofielzz = ["sd " + str(3 + 3*a) + " y", "sofiel zz"]
sd3zsofielzx = ["sd " + str(3 + 3*a) + " z", "sofiel zx"]
sd3zsofielzy = ["sd " + str(3 + 3*a) + " z", "sofiel zy"]
sd3zsofielzz = ["sd " + str(3 + 3*a) + " z", "sofiel zz"]
# - Singlet
# L-PsoKin
angmomxpsoke1 = ["angmom 1, psoke " + str(1 + a*3)]
angmomypsoke2 = ["angmom 2, psoke " + str(2 + a*3)]
angmomzpsoke3 = ["angmom 3, psoke " + str(3 + a*3)]
# Lkin-Pso
psoxozke1 = ["pso " + str(1 + 3*a) + ", ozke 1"]
psoyozke2 = ["pso " + str(2 + 3*a) + ", ozke 2"]
psozozke3 = ["pso " + str(3 + 3*a) + ", ozke 3"]
# -- Quadratic Responses
# - Triplet
# LFcSO
angmomxfcspinox = ["angmom 1, fc " + str(1 + a) + ", spinorbit 1"]
angmomyfcspinoy = ["angmom 2, fc " + str(1 + a) + ", spinorbit 2"]
angmomzfcspinoz = ["angmom 3, fc " + str(1 + a) + ", spinorbit 3"]
# LSdSOxx
angmomxsd1xspinox = ["angmom 1, sd " + str(1 + 3*a) + "1, spinorbit 1"]
angmomysd1yspinoy = ["angmom 1, sd " + str(1 + 3*a) + "2, spinorbit 2"]
angmomzsd1zspinoz = ["angmom 1, sd " + str(1 + 3*a) + "3, spinorbit 3"]
# LSdSOyy
angmomxsd2xspinox = ["angmom 2, sd " + str(2 + 3*a) + "1, spinorbit 1"]
angmomysd2yspinoy = ["angmom 2, sd " + str(2 + 3*a) + "2, spinorbit 2"]
angmomzsd2zspinoz = ["angmom 2, sd " + str(2 + 3*a) + "3, spinorbit 3"]
# LSdSOzz
angmomxsd3xspinox = ["angmom 3, sd " + str(3 + 3*a) + "1, spinorbit 1"]
angmomysd3yspinoy = ["angmom 3, sd " + str(3 + 3*a) + "2, spinorbit 2"]
angmomzsd3zspinoz = ["angmom 3, sd " + str(3 + 3*a) + "3, spinorbit 3"]
# - Singlet
# LpsoMv
angmomxpso1massvelo = ["angmom 1, pso " + str(1 + 3*a) + ", massvelo"]
angmomypso2massvelo = ["angmom 2, pso " + str(2 + 3*a) + ", massvelo"]
angmomzpso3massvelo = ["angmom 3, pso " + str(3 + 3*a) + ", massvelo"]
# LpsoDw
angmomxpso1darwin = ["angmom 1, pso " + str(1 + 3*a) + ", darwin"]
angmomypso2darwin = ["angmom 2, pso " + str(2 + 3*a) + ", darwin"]
angmomzpso3darwin = ["angmom 3, pso " + str(3 + 3*a) + ", darwin"]