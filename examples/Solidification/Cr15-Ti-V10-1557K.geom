<header>
dimension 3
grid a 2 b 2 c 512
size a 2 b 2 c 512
n_phases 2
n_materials 2
n_components 3
componentnames cr ti v
interpolation_type cubic
</header>
<solution_parameters>
finaltime 1e+1
timestep0 1e-16
timestepmin 1e-16
timestepmax 1e-3
temperature 1557.0
interfacewidth 4
initblocksize 2 2 512
initrefine 0
maxnrefine 0
initcoarsen 0
amrinterval 1000000000
reltol 1e-6
abstol 1e-6
outputfreq 10000
petscoptions -ts_adapt_monitor
outfile CrTiV
</solution_parameters>
<material 1>
chemicalenergy_type calphaddis
molarvolume 1.0e-5
mobilityc_cr 0.0
mobilityc_ti 0.0
mobilityc_v 0.0
c0 0.174 0.77 0.056
statekineticcoeff 1.0e+2
calphad_refenthalpy 0.0
calphad_unaryenthalpy_coeff_cr -8856.94 +157.48 +1.89435E-3 -1.47721E-6 +139250.0
calphad_unaryenthalpy_exp_cr 0.0 1.0 2.0 3.0 -1.0
calphad_unaryenthalpy_logCoeff_cr -26.908
calphad_unaryenthalpy_coeff_ti 6667.385 105.438379 +1.21707E-3 -0.84534E-6 -2002750.0
calphad_unaryenthalpy_exp_ti 0.0 1.0 2.0 3.0 -1.0
calphad_unaryenthalpy_logCoeff_ti -22.3771
calphad_unaryenthalpy_coeff_v -7967.842 +143.291093 +6.25E-5 -6.8E-7
calphad_unaryenthalpy_exp_v 0.0 1.0 2.0 3.0
calphad_unaryenthalpy_logCoeff_v -25.9
calphad_nbinaryenthalpy_cr_ti 2
calphad_binaryenthalpy_coeff_cr_ti_1 -2247.87 +9.14144
calphad_binaryenthalpy_exp_cr_ti_1 0.0 1.0
calphad_binaryenthalpy_coeff_cr_ti_2 198.73
calphad_binaryenthalpy_exp_cr_ti_2 0.0
calphad_nbinaryenthalpy_cr_v 4
calphad_binaryenthalpy_coeff_cr_v_1 -8253.85 -3.61592
calphad_binaryenthalpy_exp_cr_v_1 0.0 1.0
calphad_binaryenthalpy_coeff_cr_v_2 7494.82 -8.69424
calphad_binaryenthalpy_exp_cr_v_2 0.0 1.0
calphad_binaryenthalpy_coeff_cr_v_3 -17599.07 +10.13142
calphad_binaryenthalpy_exp_cr_v_3 0.0 1.0
calphad_binaryenthalpy_coeff_cr_v_4 1042.58
calphad_binaryenthalpy_exp_cr_v_4 0.0
calphad_nbinaryenthalpy_ti_v 2
calphad_binaryenthalpy_coeff_ti_v_1 6523.17
calphad_binaryenthalpy_exp_ti_v_1 0.0
calphad_binaryenthalpy_coeff_ti_v_2 2025.39
calphad_binaryenthalpy_exp_ti_v_2 0.0
calphad_nternaryenthalpy_cr_ti_v 3
calphad_iternaryenthalpy_cr_ti_v 1 2 3
calphad_ternaryenthalpy_coeff_cr_ti_v_1 -31683.57 +77.85716
calphad_ternaryenthalpy_exp_cr_ti_v_1 0.0 1.0
calphad_ternaryenthalpy_coeff_cr_ti_v_2 -29472.15 +40.11438
calphad_ternaryenthalpy_exp_cr_ti_v_2 0.0 1.0
calphad_ternaryenthalpy_coeff_cr_ti_v_3 96150.43 -50.14293
calphad_ternaryenthalpy_exp_cr_ti_v_3 0.0 1.0
</material 1>
<material 2>
chemicalenergy_type calphaddis
molarvolume 1.0e-5
mobilityc_cr 1.0
mobilityc_ti 1.0
mobilityc_v 1.0
c0 0.15 0.75 0.1
statekineticcoeff 1.0e+2
calphad_refenthalpy 0.0
calphad_unaryenthalpy_coeff_cr 15483.015 +146.059775 +1.89435E-3 -1.47721E-6 +2.37615E-21 +139250
calphad_unaryenthalpy_exp_cr 0.0 1.0 2.0 3.0 7.0 -1.0
calphad_unaryenthalpy_logCoeff_cr -26.908
calphad_unaryenthalpy_coeff_ti 369519.198 -2553.9505 -163.409355E-3 +12.457117E-6 -67034516
calphad_unaryenthalpy_exp_ti 0.0 1.0 2.0 3.0 -1.0
calphad_unaryenthalpy_logCoeff_ti +342.059267
calphad_unaryenthalpy_coeff_v 12796.275 +133.835541 +6.25E-5 -6.8E-7 -519.136E-24
calphad_unaryenthalpy_exp_v 0.0 1.0 2.0 3.0 7.0
calphad_unaryenthalpy_logCoeff_v -25.9
calphad_nbinaryenthalpy_cr_ti 3
calphad_binaryenthalpy_coeff_cr_ti_1 -365.81
calphad_binaryenthalpy_exp_cr_ti_1 0.0
calphad_binaryenthalpy_coeff_cr_ti_2 -3030.23
calphad_binaryenthalpy_exp_cr_ti_2 0.0
calphad_binaryenthalpy_coeff_cr_ti_3 1549.08
calphad_binaryenthalpy_exp_cr_ti_3 0.0
calphad_nbinaryenthalpy_cr_v 2
calphad_binaryenthalpy_coeff_cr_v_1 -20224.86
calphad_binaryenthalpy_exp_cr_v_1 0.0
calphad_binaryenthalpy_coeff_cr_v_2 -13469.90
calphad_binaryenthalpy_exp_cr_v_2 0.0
calphad_nbinaryenthalpy_ti_v 2
calphad_binaryenthalpy_coeff_ti_v_1 368.55
calphad_binaryenthalpy_exp_ti_v_1 0.0
calphad_binaryenthalpy_coeff_ti_v_2 2838.63
calphad_binaryenthalpy_exp_ti_v_2 0.0
</material 2>
<phase_material_mapping>
1
2
</phase_material_mapping>
<voxel_phase_mapping>
80 of 1
1968 of 2
</voxel_phase_mapping>
