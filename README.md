# spectroscopy
DADOS spectroscopy processing library for Walter-Hohmann-Sternwarte, Essen, Germany

the following setup is used:

Newton 30cm
DADOS
ASI 1600 MM PRO
TheSky X for Mac


processing_standardstar.ipynb
used to process standard star measurement
- additional FITS Header data
- dark current removal (currently, flat and bias are ignored)
- wavelength calibration
- instrument response function (sensitivity) including air effects
and delivers:

wavelength_solution.dat (polynomial fit for dispersion)
waves.dat: list of wavelengths
instr.dat: list of sensitivity values for waves.dat

processing_object.ipynb allows processing of object measurement
- additional FITS Header data
- dark current removal
- apply wavelength calibration
- apply instrument response
- apply athmospheric extinction

(C) 2020-2024 Dr. Michael Werger
