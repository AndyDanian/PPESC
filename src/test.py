from wave_function import wave_function as wf
from e_integral import eint as h

wfn = wf("input/H2.molden")

s = h(wfn._wfn_array)

s.integration([0, 1], ["pot"], 12)
# s.molecular_matrix(integral, sym)
