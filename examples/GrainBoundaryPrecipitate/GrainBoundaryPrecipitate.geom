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
  temperature 393.0
  interfacewidth 4
  initblocksize 2 2 2
  initrefine 5
  maxnrefine 5
  minnrefine 3
  initcoarsen 3
  amrinterval 200 
  reltol 1e-6
  abstol 1e-6 
  outputfreq 5000
  outfile GBP
  petscoptions -ts_adapt_monitor
</solution_parameters>
<material 1>
  chemicalenergy_type calphaddis
  molarvolume 1e-6
  statekineticcoeff 1e+0
  mobilityc_al 7.3451110109984390e-8
  mobilityc_cu 5.3486841665396531e-8
  mobilityc_mg 1.4696210067717206e-6
  mobilityc_zn 2.2897371183228479e-6 
  c0 0.90 0.02 0.04 0.04
  calphad_refenthalpy_coeff 0.0
  calphad_refenthalpy_exp 0
  calphad_unaryenthalpy_coeff_al -7976.15 137.093038 -1.884662e-3 -8.77664e-6 74092.0
  calphad_unaryenthalpy_exp_al 0 1 2 3 -1
  calphad_unaryenthalpy_logCoeff_al -24.3671976
  calphad_unaryenthalpy_coeff_cu -7770.458 130.485235 -2.65684e-3 1.29223e-6 52478.0
  calphad_unaryenthalpy_exp_cu 0 1 2 3 -1
  calphad_unaryenthalpy_logCoeff_cu -24.112392
  calphad_unaryenthalpy_coeff_mg -5767.34 142.775547 4.858e-3 -1.393669e-6 78950.0
  calphad_unaryenthalpy_exp_mg 0 1 2 3 -1
  calphad_unaryenthalpy_logCoeff_mg -26.1849782
  calphad_unaryenthalpy_coeff_zn -4315.967 116.900389 -1.712034e-3 1.264963e-6
  calphad_unaryenthalpy_exp_zn 0 1 2 3
  calphad_unaryenthalpy_logCoeff_zn -23.701314
  calphad_nbinaryenthalpy_al_cu 3
  calphad_binaryenthalpy_coeff_al_cu_1 -53520.0 2.0
  calphad_binaryenthalpy_exp_al_cu_1 0 1
  calphad_binaryenthalpy_coeff_al_cu_2 38590.0 -2.0
  calphad_binaryenthalpy_exp_al_cu_2 0 1
  calphad_binaryenthalpy_coeff_al_cu_3 1170.0
  calphad_binaryenthalpy_exp_al_cu_3 0
  calphad_nbinaryenthalpy_al_mg 3
  calphad_binaryenthalpy_coeff_al_mg_1 4971.0 -3.500025
  calphad_binaryenthalpy_exp_al_mg_1 0 1
  calphad_binaryenthalpy_coeff_al_mg_2 225.0 0.1
  calphad_binaryenthalpy_exp_al_mg_2 0 1
  calphad_binaryenthalpy_coeff_al_mg_3 950.0 2.5e-5
  calphad_binaryenthalpy_exp_al_mg_3 0 1
  calphad_nbinaryenthalpy_al_zn 3
  calphad_binaryenthalpy_coeff_al_zn_1 7297.48 0.47512
  calphad_binaryenthalpy_exp_al_cu_1 0 1
  calphad_binaryenthalpy_coeff_al_zn_2 6612.88 -4.59110
  calphad_binaryenthalpy_exp_al_cu_2 0 1
  calphad_binaryenthalpy_coeff_al_zn_3 -3097.19 3.30635
  calphad_binaryenthalpy_exp_al_cu_3 0 1
  calphad_nbinaryenthalpy_cu_mg 1
  calphad_binaryenthalpy_coeff_cu_mg_1 -22279.28 5.868
  calphad_binaryenthalpy_exp_cu_mg_1 0 1
  calphad_nbinaryenthalpy_cu_zn 3
  calphad_binaryenthalpy_coeff_cu_zn_1 -42803.75 10.02258
  calphad_binaryenthalpy_exp_cu_zn_1 0 1
  calphad_binaryenthalpy_coeff_cu_zn_2 2936.39 -3.05323
  calphad_binaryenthalpy_exp_cu_zn_2 0 1
  calphad_binaryenthalpy_coeff_cu_zn_3 9034.20 -5.39314
  calphad_binaryenthalpy_exp_cu_zn_3 0 1
  calphad_nbinaryenthalpy_mg_zn 2
  calphad_binaryenthalpy_coeff_mg_zn_1 -3056.82 5.63801
  calphad_binaryenthalpy_exp_mg_zn_1 0 1
  calphad_binaryenthalpy_coeff_mg_zn_2 -3127.26 5.65563
  calphad_binaryenthalpy_exp_mg_zn_2 0 1
  calphad_nternaryenthalpy_al_cu_mg 1
  calphad_ternaryenthalpy_coeff_al_cu_mg_1 60000.0
  calphad_ternaryenthalpy_exp_al_cu_mg_1 0
  calphad_nternaryenthalpy_al_cu_zn 3
  calphad_ternaryenthalpy_coeff_al_cu_zn_1 7.2990788e+4
  calphad_ternaryenthalpy_exp_al_cu_zn_1 0
  calphad_ternaryenthalpy_coeff_al_cu_zn_2 1.6799927e+4
  calphad_ternaryenthalpy_exp_al_cu_zn_2 0
  calphad_ternaryenthalpy_coeff_al_cu_zn_3 5.9580949e+4
  calphad_ternaryenthalpy_exp_al_cu_zn_3 0
</material 1>
<material 2>
  chemicalenergy_type quadratic
  molarvolume 1e-6
  statekineticcoeff 1e+0
  mobilityc_al 7.3451110109984390e-8
  mobilityc_cu 5.3486841665396531e-8
  mobilityc_mg 1.4696210067717206e-6
  mobilityc_zn 2.2897371183228479e-6 
  c0 0.16757 0.04471 0.33333 0.45439   
  quad_ceq_coeff_al 0.16757
  quad_ceq_exp_al 0
  quad_ceq_coeff_cu 0.04471
  quad_ceq_exp_cu 0
  quad_ceq_coeff_mg 0.33333
  quad_ceq_exp_mg 0
  quad_ceq_coeff_zn 0.45439   
  quad_ceq_exp_zn 0
  quad_refenthalpy_coeff -2.6605e+04
  quad_refenthalpy_exp 0
  quad_binaryenthalpy_coeff_al 1.0e+6
  quad_binaryenthalpy_exp_al 0
  quad_binaryenthalpy_coeff_cu 1.0e+6
  quad_binaryenthalpy_exp_cu 0
  quad_binaryenthalpy_coeff_mg 1.0e+6
  quad_binaryenthalpy_exp_mg 0
  quad_binaryenthalpy_coeff_zn 1.0e+6 
  quad_binaryenthalpy_exp_zn 0
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
