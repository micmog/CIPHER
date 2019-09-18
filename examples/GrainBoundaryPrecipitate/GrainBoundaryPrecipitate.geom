<header>
  dimension 3
  grid a 64 b 64 c 64
  size a 64 b 64 c 64
  n_phases 3
  n_materials 2
  n_components 4
  componentnames al cu mg zn
  interpolation_type cubic
</header>
<solution_parameters>
  finaltime 72000.0 
  timestep0 1e-9
  timestepmin 1e-9
  timestepmax 1e+2 
  interfacewidth 1 
  initrefine 5
  maxnrefine 7
  initcoarsen 4
  amrinterval 200 
  reltol 1e-4 
  abstol 1e-4 
  outputfreq 5000
  outfile GBP
  petscoptions -ts_adapt_monitor
</solution_parameters>
<material 1>
  chemicalenergy_type calphaddis
  RT 3267.402
  molarvolume 1e-6
  mobilityc 7.3451110109984390e-8 5.3486841665396531e-8 1.4696210067717206e-6 2.2897371183228479e-6 
  c0 0.90 0.02 0.04 0.04
  calphad_refenthalpy 0.0
  calphad_unaryenthalpy -11469.906443934595 -13367.567247075782 -10232.482939321606 -14359.062357670766  
  calphad_nbinaryenthalpy 1 2 3
  calphad_binaryenthalpy 1 2 -52734.0 37804 1170
  calphad_nbinaryenthalpy 1 3 3
  calphad_binaryenthalpy 1 3 3595.5 1066.239 950.0
  calphad_nbinaryenthalpy 1 4 3
  calphad_binaryenthalpy 1 4 7484.20216 4808.5777 -1797.79445
  calphad_nbinaryenthalpy 2 3 1
  calphad_binaryenthalpy 2 3 -19846.10824 
  calphad_nbinaryenthalpy 2 4 3
  calphad_binaryenthalpy 2 4 -38864.87606 1736.4706099999999 6914.69598
  calphad_nbinaryenthalpy 3 4 2
  calphad_binaryenthalpy 3 4 -841.0820699999999 -904.5974099999999
  calphad_nternaryenthalpy 1 2 4 3
  calphad_iternaryenthalpy 1 2 4 1 2 4
  calphad_ternaryenthalpy 1 2 4 7.2990788e+04 1.6799927e+04 5.9580949e+04
</material 1>
<material 2>
  chemicalenergy_type quadratic
  molarvolume 1e-6
  mobilityc 7.3451110109984390e-8 5.3486841665396531e-8 1.4696210067717206e-6 2.2897371183228479e-6 
  c0 0.16757 0.04471 0.33333 0.45439   
  quad_ceq 0.16757 0.04471 0.33333 0.45439   
  quad_refenthalpy -2.6605e+04
  quad_unaryenthalpy  0.0 0.0 0.0 0.0
  quad_binaryenthalpy 1.0e+6 1.0e+6 1.0e+6 1.0e+6 
</material 2>
<phase_material_mapping>
  2 of 1
  2
</phase_material_mapping>
<voxel_phase_mapping>
  116510 of 1
  4 of 3
  59 of 1
  6 of 3
  57 of 1
  8 of 3
  56 of 1
  8 of 3
  56 of 1
  8 of 3
  56 of 1
  8 of 3
  57 of 1
  6 of 3
  59 of 1
  4 of 3
  3515 of 1
  6 of 3
  57 of 1
  8 of 3
  55 of 1
  10 of 3
  53 of 1
  12 of 3
  52 of 1
  12 of 3
  52 of 1
  12 of 3
  52 of 1
  12 of 3
  52 of 1
  12 of 3
  52 of 1
  12 of 3
  53 of 1
  10 of 3
  55 of 1
  8 of 3
  57 of 1
  6 of 3
  3321 of 1
  8 of 3
  55 of 1
  10 of 3
  53 of 1
  12 of 3
  51 of 1
  14 of 3
  50 of 1
  14 of 3
  50 of 1
  14 of 3
  50 of 1
  14 of 3
  50 of 1
  14 of 3
  50 of 1
  14 of 3
  50 of 1
  14 of 3
  50 of 1
  14 of 3
  51 of 1
  12 of 3
  53 of 1
  10 of 3
  55 of 1
  8 of 3
  3193 of 1
  6 of 3
  56 of 1
  10 of 3
  53 of 1
  12 of 3
  51 of 1
  14 of 3
  50 of 1
  14 of 3
  49 of 1
  16 of 3
  48 of 1
  16 of 3
  48 of 1
  16 of 3
  48 of 1
  16 of 3
  48 of 1
  16 of 3
  48 of 1
  16 of 3
  49 of 1
  14 of 3
  50 of 1
  14 of 3
  51 of 1
  12 of 3
  53 of 1
  10 of 3
  56 of 1
  6 of 3
  1565 of 1
  1565 of 2
  6 of 3
  56 of 2
  10 of 3
  53 of 2
  12 of 3
  51 of 2
  14 of 3
  50 of 2
  14 of 3
  49 of 2
  16 of 3
  48 of 2
  16 of 3
  48 of 2
  16 of 3
  48 of 2
  16 of 3
  48 of 2
  16 of 3
  48 of 2
  16 of 3
  49 of 2
  14 of 3
  50 of 2
  14 of 3
  51 of 2
  12 of 3
  53 of 2
  10 of 3
  56 of 2
  6 of 3
  3193 of 2
  8 of 3
  55 of 2
  10 of 3
  53 of 2
  12 of 3
  51 of 2
  14 of 3
  50 of 2
  14 of 3
  50 of 2
  14 of 3
  50 of 2
  14 of 3
  50 of 2
  14 of 3
  50 of 2
  14 of 3
  50 of 2
  14 of 3
  50 of 2
  14 of 3
  51 of 2
  12 of 3
  53 of 2
  10 of 3
  55 of 2
  8 of 3
  3321 of 2
  6 of 3
  57 of 2
  8 of 3
  55 of 2
  10 of 3
  53 of 2
  12 of 3
  52 of 2
  12 of 3
  52 of 2
  12 of 3
  52 of 2
  12 of 3
  52 of 2
  12 of 3
  52 of 2
  12 of 3
  53 of 2
  10 of 3
  55 of 2
  8 of 3
  57 of 2
  6 of 3
  3515 of 2
  4 of 3
  59 of 2
  6 of 3
  57 of 2
  8 of 3
  56 of 2
  8 of 3
  56 of 2
  8 of 3
  56 of 2
  8 of 3
  57 of 2
  6 of 3
  59 of 2
  4 of 3
  116510 of 2
</voxel_phase_mapping>
