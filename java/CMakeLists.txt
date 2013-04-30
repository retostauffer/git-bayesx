set(ANDREA_SRC 
	../andrea/baseline.cpp
	../andrea/cox.cpp
	../andrea/IWLS_baseline.cpp
	../andrea/multibaseline.cpp
	../andrea/multistate.cpp
)
set(JAVA_SRC
        ../bib/adminparse.cpp
        ../bib/adminparse_basic.cpp
        ../bib/adminparse_pointers.cpp
        ../graph/graphobj.cpp
)
set(BIB_SRC 
	../bib/Random.cpp
	# ../bib/adminparse_gnu.cpp
	../bib/bandmat.cpp
	../bib/bandmat_penalty.cpp
	../bib/bayesreg.cpp
	../bib/bayesreg2.cpp
	../bib/bayesreg3.cpp
	../bib/clstring.cpp
	../bib/command.cpp
	../bib/data.cpp
	../bib/dataobj.cpp
	../bib/envmatrix.cpp
	../bib/envmatrix_penalty.cpp
	../bib/graph.cpp
	../bib/map.cpp
	../bib/mapobject.cpp
	../bib/model.cpp
	../bib/model_remlreg.cpp
	../bib/model_stepwise.cpp
	../bib/option.cpp
	../bib/realobs.cpp
	../bib/realvar.cpp
	../bib/remlreg.cpp
	../bib/sparsemat.cpp
	../bib/statmat.cpp
	../bib/statmat_penalty.cpp
	../bib/statobj.cpp
	../bib/stepwisereg.cpp
	../bib/use.cpp
	../bib/vectorn.cpp
)
set(DAG_SRC
	../dag/adjacency.cpp
	../dag/dagobject.cpp
	../dag/fullcond_dag.cpp
	../dag/fullcond_dag_d.cpp
	../dag/fullcond_dag_ia.cpp
	../dag/fullcond_dag_ia_mixed.cpp
	../dag/fullcond_rj.cpp
	../dag/fullcond_rj_int.cpp
	../dag/fullcond_rj_mix.cpp
	../dag/func_dag.cpp
	../dag/ia.cpp
	../dag/ia_mixed.cpp
)
set(LEYRE_SRC 
	../leyre/nbinomial.cpp
	../leyre/zip.cpp
)
set(MCMC_SRC 
	../mcmc/baseline_reml.cpp
	../mcmc/distribution.cpp
	../mcmc/fullcond.cpp
	../mcmc/fullcond_merror.cpp
	../mcmc/fullcond_mult.cpp
	../mcmc/fullcond_nonp_gaussian.cpp
	../mcmc/fullcond_nonp_gaussian_stepwise.cpp
	../mcmc/gaussian_heteroskedastic.cpp
	../mcmc/hrandom.cpp
	../mcmc/kriging2.cpp
	../mcmc/kriging.cpp
	../mcmc/mcmc.cpp
	../mcmc/mcmc_const.cpp
	../mcmc/mcmc_const_stepwise.cpp
	../mcmc/mcmc_nonp.cpp
	../mcmc/mcmc_nonpbasis.cpp
	../mcmc/mcmcsimul.cpp
	../mcmc/mcmcsimul2.cpp
	../mcmc/mcmcsimul2_multi.cpp
	../mcmc/randomeffect.cpp
	../mcmc/randomeffect_stepwise.cpp
	../mcmc/remlest.cpp
	../mcmc/remlest_multi.cpp
	../mcmc/remlest_multi2.cpp
	../mcmc/remlest_multi3.cpp
	../mcmc/tvariance.cpp
	../mcmc/tvariance2dim.cpp
	../mcmc/variance_nonp.cpp
	../mcmc/variance_nonp_vector.cpp
	../mcmc/variance_nonp_vector_nigmix.cpp
)
set(PSPLINES_SRC 
	../psplines/IWLS_pspline.cpp
	../psplines/bsplinemat.cpp
	../psplines/fullcond_pspline_gaussian.cpp
	../psplines/fullcond_pspline_stepwise.cpp
	../psplines/fullcond_pspline_surf_gaussian.cpp
	../psplines/fullcond_pspline_surf_stepwise.cpp
	../psplines/mcmc_pspline.cpp
	../psplines/mcmc_pspline_surf.cpp
	../psplines/spline_basis.cpp
	../psplines/spline_basis_surf.cpp
)
set(STRUCTADD_SRC 
	../structadd/FC.cpp
	../structadd/FC_hrandom.cpp
	../structadd/FC_hrandom_variance.cpp
	../structadd/FC_hrandom_variance_vec.cpp
	../structadd/FC_hrandom_variance_vec_nmig.cpp
	../structadd/FC_linear.cpp
	../structadd/FC_mult.cpp
	../structadd/FC_nonp.cpp
	../structadd/FC_nonp_variance.cpp
	../structadd/FC_predict.cpp
	../structadd/FC_predict_mult.cpp
	../structadd/FC_cv.cpp
	../structadd/FC_variance_pen_vector.cpp
	../structadd/GENERAL_OPTIONS.cpp
	../structadd/MASTER_obj.cpp
	../structadd/design.cpp
	../structadd/design_hrandom.cpp
	../structadd/design_mrf.cpp
	../structadd/design_pspline.cpp
	../structadd/distr.cpp
	../structadd/distr_categorical.cpp
	../structadd/distr_categorical_mult.cpp
	../structadd/distr_mixture.cpp
	../structadd/distr_gamlss.cpp
	../structadd/distr_gamlss_nadja.cpp
	../structadd/distr_zeroadjusted.cpp
        ../structadd/mcmcsim.cpp
	../structadd/model_parameters.cpp
	../structadd/superbayesreg.cpp
	../structadd/design_kriging.cpp
	../structadd/FC_predictive_check.cpp
        ../structadd/FC_predict_predictor.cpp
        ../structadd/FC_nonp_variance_vec.cpp
)
set(SRC 
	${ANDREA_SRC}
	${BIB_SRC}
	${DAG_SRC}
	${LEYRE_SRC}
	${MCMC_SRC}
	${PSPLINES_SRC}
	${STRUCTADD_SRC}
	# ../main.cpp
       	../samson/multgaussian.cpp
	../adaptiv/fullcond_adaptiv.cpp
	../alex/mixture.cpp
)

if(CMAKE_SIZEOF_VOID_P MATCHES 8)
  set(64BIT 1)
else()
  set(32BIT 1)
endif()

add_definitions(-DJAVA_OUTPUT_WINDOW)
add_library(BayesXdll SHARED ${SRC} ${JAVA_SRC} ../java/BayesXdll.cpp)
add_definitions(-D__BUILDING_GNU -DTEMPL_INCL_DEF -D_MSC_VER2 -DNO_TEMPLATE_FRIENDS -DINCLUDE_REML -DINCLUDE_MCMC)

if(NOT WIN32)
  add_definitions(-D__BUILDING_LINUX)
  find_package(Readline REQUIRED)
  target_link_libraries(BayesXdll ${READLINE_LIBRARY})
else(WIN32)
  add_definitions(-D_JNI_IMPLEMENTATION_)
  set_target_properties(BayesXdll PROPERTIES PREFIX "") 
  if(MINGW AND 32BIT)
    set_target_properties(BayesXdll PROPERTIES LINK_FLAGS "-Wl,--kill-at") 
  endif()
endif()

find_package(JNI REQUIRED)
include_directories(${JNI_INCLUDE_DIRS})
include_directories(../bib ../alex ../adaptiv ../andrea ../dag ../graph ../mcmc ../psplines ../samson ../leyre ../structadd ../java)

set(JavaSources 
BayesX.java
Dense.java
FileFormatException.java
File_Writer.java
HelpFilter.java
Interpolate.java
MapFilter.java
MapPanel.java
Matrix3D.java
Model3D.java
OutputFilter.java
OutputFilter2.java
PSGr2.java
PSGrBase.java
PageableText.java
ParseData.java
ParseThread.java
Plot3D.java
PrintableComponent.java
Sort.java
Sort_ind.java
Surfplot.java
SwingWorker.java
ThreeD.java
finalmap.java
)
find_package(Java)
include(UseJava)

# set(CMAKE_JAVA_EXECUTABLE BayesX)
# set(CMAKE_JAVA_MANIFEST manifest.mf)
#set(CMAKE_JAVA_JAR_ENTRY_POINT BayesX)
add_jar(BayesX ${JavaSources}) 


