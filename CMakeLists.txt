project(BayesX)
cmake_minimum_required(VERSION 2.8.9)

option(WithoutDefaultOs "Compiles BayesX without adding -Os to definitions" OFF)
option(WithoutGSL "Compiles BayesX without gsl" OFF)

set(CMAKE_MODULE_PATH ${BayesX_SOURCE_DIR}/share/cmake ${CMAKE_MODULE_PATH})
set(ANDREA_SRC
	andrea/baseline.cpp
	andrea/cox.cpp
	andrea/IWLS_baseline.cpp
	andrea/multibaseline.cpp
	andrea/multistate.cpp
)
set(BIB_SRC
	bib/Random.cpp
	bib/adminparse_gnu.cpp
	bib/bandmat.cpp
	bib/bandmat_penalty.cpp
	bib/bayesreg.cpp
	bib/bayesreg2.cpp
	bib/bayesreg3.cpp
	bib/clstring.cpp
	bib/command.cpp
	bib/data.cpp
	bib/dataobj.cpp
	bib/envmatrix.cpp
	bib/envmatrix_penalty.cpp
	bib/graph.cpp
	bib/map.cpp
	bib/mapobject.cpp
	bib/model.cpp
	bib/model_remlreg.cpp
	bib/model_stepwise.cpp
	bib/option.cpp
	bib/realobs.cpp
	bib/realvar.cpp
	bib/remlreg.cpp
	bib/sparsemat.cpp
	bib/statmat.cpp
	bib/statmat_penalty.cpp
	bib/statobj.cpp
	bib/stepwisereg.cpp
	bib/use.cpp
	bib/vectorn.cpp
)
set(DAG_SRC
	dag/adjacency.cpp
	dag/dagobject.cpp
	dag/fullcond_dag.cpp
	dag/fullcond_dag_d.cpp
	dag/fullcond_dag_ia.cpp
	dag/fullcond_dag_ia_mixed.cpp
	dag/fullcond_rj.cpp
	dag/fullcond_rj_int.cpp
	dag/fullcond_rj_mix.cpp
	dag/func_dag.cpp
	dag/ia.cpp
	dag/ia_mixed.cpp
)
set(LEYRE_SRC
	leyre/nbinomial.cpp
	leyre/zip.cpp
)
set(MCMC_SRC
	mcmc/baseline_reml.cpp
	mcmc/distribution.cpp
	mcmc/fullcond.cpp
	mcmc/fullcond_merror.cpp
	mcmc/fullcond_mult.cpp
	mcmc/fullcond_nonp_gaussian.cpp
	mcmc/fullcond_nonp_gaussian_stepwise.cpp
	mcmc/gaussian_heteroskedastic.cpp
	mcmc/hrandom.cpp
	mcmc/kriging2.cpp
	mcmc/kriging.cpp
	mcmc/mcmc.cpp
	mcmc/mcmc_const.cpp
	mcmc/mcmc_const_stepwise.cpp
	mcmc/mcmc_nonp.cpp
	mcmc/mcmc_nonpbasis.cpp
	mcmc/mcmcsimul.cpp
	mcmc/mcmcsimul2.cpp
	mcmc/mcmcsimul2_multi.cpp
	mcmc/randomeffect.cpp
	mcmc/randomeffect_stepwise.cpp
	mcmc/remlest.cpp
	mcmc/remlest_multi.cpp
	mcmc/remlest_multi2.cpp
	mcmc/remlest_multi3.cpp
	mcmc/tvariance.cpp
	mcmc/tvariance2dim.cpp
	mcmc/variance_nonp.cpp
	mcmc/variance_nonp_vector.cpp
	mcmc/variance_nonp_vector_nigmix.cpp
)
set(PSPLINES_SRC
	psplines/IWLS_pspline.cpp
	psplines/bsplinemat.cpp
	psplines/fullcond_pspline_gaussian.cpp
	psplines/fullcond_pspline_stepwise.cpp
	psplines/fullcond_pspline_surf_gaussian.cpp
	psplines/fullcond_pspline_surf_stepwise.cpp
	psplines/mcmc_pspline.cpp
	psplines/mcmc_pspline_surf.cpp
	psplines/spline_basis.cpp
	psplines/spline_basis_surf.cpp
)
set(STRUCTADD_SRC
	structadd/FC.cpp
	structadd/FC_hrandom.cpp
	structadd/FC_hrandom_variance.cpp
	structadd/FC_hrandom_variance_vec.cpp
	structadd/FC_hrandom_variance_vec_nmig.cpp
	structadd/FC_linear.cpp
	structadd/FC_merror.cpp
	structadd/FC_mult.cpp
	structadd/FC_mult.cpp
	structadd/FC_nonp.cpp
	structadd/FC_nonp_variance.cpp
	structadd/FC_nonp_variance_vec.cpp
	structadd/FC_predict.cpp
	structadd/FC_predict_mult.cpp
	structadd/FC_cv.cpp
	structadd/FC_variance_pen_vector.cpp
	structadd/GENERAL_OPTIONS.cpp
	structadd/MASTER_obj.cpp
	structadd/design.cpp
	structadd/design_hrandom.cpp
	structadd/design_kriging.cpp
	structadd/design_mrf.cpp
	structadd/design_pspline.cpp
	structadd/design_userdefined.cpp
	structadd/distr.cpp
	structadd/distr_categorical.cpp
	structadd/distr_categorical_mult.cpp
	structadd/distr_mixture.cpp
	structadd/distr_gamlss.cpp
	structadd/distr_gamlss_copula.cpp
	structadd/distr_gamlss_nadja.cpp
	structadd/distr_zeroadjusted.cpp
        structadd/mcmcsim.cpp
	structadd/model_parameters.cpp
	structadd/superbayesreg.cpp
	structadd/FC_predictive_check.cpp
        structadd/FC_predict_predictor.cpp
)
set(SRC
	${ANDREA_SRC}
	${BIB_SRC}
	${DAG_SRC}
	${LEYRE_SRC}
	${MCMC_SRC}
	${PSPLINES_SRC}
	${STRUCTADD_SRC}
	main.cpp
        samson/multgaussian.cpp
	adaptiv/fullcond_adaptiv.cpp
	alex/mixture.cpp
)

if(CMAKE_SIZEOF_VOID_P MATCHES 8)
  set(64BIT 1)
else()
  set(32BIT 1)
endif()

include_directories(bib alex adaptiv andrea dag graph mcmc psplines samson leyre structadd)
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
link_libraries( ${GSL_LIBRARIES}  ${GSLCBLAS_LIBRARIES})
add_definitions(${GSL_CFLAGS})
add_executable(bayesx ${SRC})

add_definitions(-D__BUILDING_GNU -DTEMPL_INCL_DEF -D_MSC_VER2 -DNO_TEMPLATE_FRIENDS -DINCLUDE_REML -DINCLUDE_MCMC)
if(NOT WithoutDefaultOs)
  add_definitions(-Os)
endif()
if(NOT WithoutGSL)
  add_definitions(-DBayesX_gsl_included -DBAYESX_GSL_INCLUDED)
endif()

if(32BIT)
  add_definitions(-m32)
else()
  add_definitions(-m64)
endif()


if(NOT WIN32)
add_definitions(-D__BUILDING_LINUX)
find_package(Readline REQUIRED)
target_link_libraries(bayesx ${READLINE_LIBRARY})
endif()
