"""
A library of functions (for now just one function) to create a thermal
(blackbody) spectrum.
"""

def planck_spec_lambda(T, lammin, lammax, showfig=True, returnspec=False):
   """
   Generates a Planck spectrum in terms of wavelength, i.e., B_lambda(T).

   Inputs:
      T          - temperature in K
      lammin     - minimum wavelength in m
      lammax     - maximum wavelength in m
      showfig    - show the figure (default=True)
      returnspec - return the Planck spectrum (default=False)
   """
