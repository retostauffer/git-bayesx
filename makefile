# Project: BayesX

CPP  = g++
CC   = gcc
RES  = 
OBJ  = gnuobj/clstring.o gnuobj/adminparse_gnu.o gnuobj/Random.o gnuobj/dataobj.o gnuobj/model.o gnuobj/statobj.o gnuobj/option.o gnuobj/use.o gnuobj/bandmat.o gnuobj/bandmat_penalty.o gnuobj/bayesreg2.o gnuobj/bayesreg3.o gnuobj/bayesreg.o gnuobj/data.o gnuobj/envmatrix.o gnuobj/envmatrix_penalty.o gnuobj/graph.o gnuobj/map.o gnuobj/tvariance2dim.o gnuobj/tvariance.o gnuobj/variance_nonp.o gnuobj/variance_nonp_vector.o gnuobj/baseline_reml.o gnuobj/distribution.o gnuobj/fullcond.o gnuobj/fullcond_merror.o gnuobj/fullcond_mult.o gnuobj/fullcond_nonp_gaussian.o gnuobj/fullcond_nonp_gaussian_stepwise.o gnuobj/gaussian_heteroskedastic.o gnuobj/kriging2.o gnuobj/kriging.o gnuobj/mcmc.o gnuobj/mcmc_const.o gnuobj/mcmc_const_stepwise.o gnuobj/mcmc_nonp.o gnuobj/mcmc_nonpbasis.o gnuobj/mcmcsimul2.o gnuobj/mcmcsimul2_multi.o gnuobj/mcmcsimul.o gnuobj/randomeffect.o gnuobj/randomeffect_stepwise.o gnuobj/remlest.o gnuobj/remlest_multi2.o gnuobj/remlest_multi3.o gnuobj/remlest_multi.o gnuobj/fullcond_adaptiv.o gnuobj/mixture.o gnuobj/baseline.o gnuobj/cox.o gnuobj/IWLS_baseline.o gnuobj/multibaseline.o gnuobj/multistate.o gnuobj/nbinomial.o gnuobj/zip.o gnuobj/multgaussian.o gnuobj/bsplinemat.o gnuobj/fullcond_pspline_gaussian.o gnuobj/fullcond_pspline_stepwise.o gnuobj/fullcond_pspline_surf_gaussian.o gnuobj/fullcond_pspline_surf_stepwise.o gnuobj/IWLS_pspline.o gnuobj/mcmc_pspline.o gnuobj/mcmc_pspline_surf.o gnuobj/spline_basis.o gnuobj/spline_basis_surf.o gnuobj/command.o gnuobj/sparsemat.o gnuobj/statmat.o gnuobj/statmat_penalty.o gnuobj/realobs.o gnuobj/realvar.o gnuobj/remlreg.o gnuobj/model_remlreg.o gnuobj/stepwisereg.o gnuobj/model_stepwise.o gnuobj/adjacency.o gnuobj/dagobject.o gnuobj/fullcond_dag.o gnuobj/fullcond_dag_d.o gnuobj/fullcond_dag_ia.o gnuobj/fullcond_dag_ia_mixed.o gnuobj/fullcond_rj.o gnuobj/fullcond_rj_int.o gnuobj/fullcond_rj_mix.o gnuobj/func_dag.o gnuobj/ia.o gnuobj/ia_mixed.o gnuobj/mapobject.o gnuobj/hrandom.o gnuobj/variance_nonp_vector_nigmix.o gnuobj/FC.o gnuobj/FC_hrandom.o gnuobj/FC_hrandom_variance.o gnuobj/FC_linear.o gnuobj/FC_mult.o gnuobj/FC_nonp.o gnuobj/FC_nonp_variance.o gnuobj/FC_predict.o gnuobj/GENERAL_OPTIONS.o gnuobj/MASTER_obj.o gnuobj/design.o gnuobj/design_hrandom.o gnuobj/design_mrf.o gnuobj/design_pspline.o gnuobj/distr.o gnuobj/distr_categorical.o gnuobj/mcmcsim.o gnuobj/model_parameters.o gnuobj/superbayesreg.o gnuobj/design_kriging.o gnuobj/FC_predictive_check.o gnuobj/main.o $(RES)
LINKOBJ  = gnuobj/main.o gnuobj/clstring.o gnuobj/adminparse_gnu.o gnuobj/Random.o gnuobj/dataobj.o gnuobj/model.o gnuobj/statobj.o gnuobj/option.o gnuobj/use.o gnuobj/bandmat.o gnuobj/bandmat_penalty.o gnuobj/bayesreg2.o gnuobj/bayesreg3.o gnuobj/bayesreg.o gnuobj/data.o gnuobj/envmatrix.o gnuobj/envmatrix_penalty.o gnuobj/graph.o gnuobj/map.o gnuobj/tvariance2dim.o gnuobj/tvariance.o gnuobj/variance_nonp.o gnuobj/variance_nonp_vector.o gnuobj/baseline_reml.o gnuobj/distribution.o gnuobj/fullcond.o gnuobj/fullcond_merror.o gnuobj/fullcond_mult.o gnuobj/fullcond_nonp_gaussian.o gnuobj/fullcond_nonp_gaussian_stepwise.o gnuobj/gaussian_heteroskedastic.o gnuobj/kriging2.o gnuobj/kriging.o gnuobj/mcmc.o gnuobj/mcmc_const.o gnuobj/mcmc_const_stepwise.o gnuobj/mcmc_nonp.o gnuobj/mcmc_nonpbasis.o gnuobj/mcmcsimul2.o gnuobj/mcmcsimul2_multi.o gnuobj/mcmcsimul.o gnuobj/randomeffect.o gnuobj/randomeffect_stepwise.o gnuobj/remlest.o gnuobj/remlest_multi2.o gnuobj/remlest_multi3.o gnuobj/remlest_multi.o gnuobj/fullcond_adaptiv.o gnuobj/mixture.o gnuobj/baseline.o gnuobj/cox.o gnuobj/IWLS_baseline.o gnuobj/multibaseline.o gnuobj/multistate.o gnuobj/nbinomial.o gnuobj/zip.o gnuobj/multgaussian.o gnuobj/bsplinemat.o gnuobj/fullcond_pspline_gaussian.o gnuobj/fullcond_pspline_stepwise.o gnuobj/fullcond_pspline_surf_gaussian.o gnuobj/fullcond_pspline_surf_stepwise.o gnuobj/IWLS_pspline.o gnuobj/mcmc_pspline.o gnuobj/mcmc_pspline_surf.o gnuobj/spline_basis.o gnuobj/spline_basis_surf.o gnuobj/command.o gnuobj/sparsemat.o gnuobj/statmat.o gnuobj/statmat_penalty.o gnuobj/realobs.o gnuobj/realvar.o gnuobj/remlreg.o gnuobj/model_remlreg.o gnuobj/stepwisereg.o gnuobj/model_stepwise.o gnuobj/adjacency.o gnuobj/dagobject.o gnuobj/fullcond_dag.o gnuobj/fullcond_dag_d.o gnuobj/fullcond_dag_ia.o gnuobj/fullcond_dag_ia_mixed.o gnuobj/fullcond_rj.o gnuobj/fullcond_rj_int.o gnuobj/fullcond_rj_mix.o gnuobj/func_dag.o gnuobj/ia.o gnuobj/ia_mixed.o gnuobj/mapobject.o gnuobj/hrandom.o gnuobj/variance_nonp_vector_nigmix.o gnuobj/FC.o gnuobj/FC_hrandom.o gnuobj/FC_hrandom_variance.o gnuobj/FC_linear.o gnuobj/FC_mult.o gnuobj/FC_nonp.o gnuobj/FC_nonp_variance.o gnuobj/FC_predict.o gnuobj/GENERAL_OPTIONS.o gnuobj/MASTER_obj.o gnuobj/design.o gnuobj/design_hrandom.o gnuobj/design_mrf.o gnuobj/design_pspline.o gnuobj/distr.o gnuobj/distr_categorical.o gnuobj/mcmcsim.o gnuobj/model_parameters.o gnuobj/superbayesreg.o gnuobj/design_kriging.o gnuobj/FC_predictive_check.o  /usr/lib/libreadline.so $(RES)
LIBS =  -L"/usr/lib" -L"usr/share/readline" -g3
INCS =  -I"bib"  -I"alex"  -I"adaptiv"  -I"andrea"  -I"dag"  -I"graph"  -I"mcmc"  -I"psplines"  -I"samson"  -I"leyre"  -I"structadd" -I"/usr/include/c++/3.4/backward/" -I"/usr/local/include/readline/"
CXXINCS =  -I"bib"  -I"alex"  -I"adaptiv"  -I"andrea"  -I"dag"  -I"graph"  -I"mcmc"  -I"psplines"  -I"samson"  -I"leyre"  -I"structadd" -I"/usr/include/c++/3.4/backward/" -I"/usr/local/include/readline/"

BIN  = BayesX
CXXFLAGS = $(CXXINCS)   -w -D__BUILDING_GNU -DTEMPL_INCL_DEF -D_MSC_VER2 -DNO_TEMPLATE_FRIENDS -DINCLUDE_REML -DINCLUDE_MCMC -D__BUILDING_LINUX -g3
CFLAGS = $(INCS)   -w -D__BUILDING_GNU -DTEMPL_INCL_DEF -D_MSC_VER2 -DNO_TEMPLATE_FRIENDS -DINCLUDE_REML -DINCLUDE_MCMC -D__BUILDING_LINUX -g3
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before BayesX all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o "BayesX" $(LIBS) 

gnuobj/fullcond_adaptiv.o: adaptiv/fullcond_adaptiv.cpp
	$(CPP) -c adaptiv/fullcond_adaptiv.cpp -o gnuobj/fullcond_adaptiv.o $(CXXFLAGS)


gnuobj/mixture.o: alex/mixture.cpp
	$(CPP) -c alex/mixture.cpp -o gnuobj/mixture.o $(CXXFLAGS)


gnuobj/baseline.o: andrea/baseline.cpp
	$(CPP) -c andrea/baseline.cpp -o gnuobj/baseline.o $(CXXFLAGS)

gnuobj/cox.o: andrea/cox.cpp
	$(CPP) -c andrea/cox.cpp -o gnuobj/cox.o $(CXXFLAGS)

gnuobj/IWLS_baseline.o: andrea/IWLS_baseline.cpp
	$(CPP) -c andrea/IWLS_baseline.cpp -o gnuobj/IWLS_baseline.o $(CXXFLAGS)

gnuobj/multibaseline.o: andrea/multibaseline.cpp
	$(CPP) -c andrea/multibaseline.cpp -o gnuobj/multibaseline.o $(CXXFLAGS)

gnuobj/multistate.o: andrea/multistate.cpp
	$(CPP) -c andrea/multistate.cpp -o gnuobj/multistate.o $(CXXFLAGS)


gnuobj/Random.o: bib/Random.cpp
	$(CPP) -c bib/Random.cpp -o gnuobj/Random.o $(CXXFLAGS)

gnuobj/adminparse_gnu.o: bib/adminparse_gnu.cpp
	$(CPP) -c bib/adminparse_gnu.cpp -o gnuobj/adminparse_gnu.o $(CXXFLAGS)

gnuobj/bandmat.o: bib/bandmat.cpp
	$(CPP) -c bib/bandmat.cpp -o gnuobj/bandmat.o $(CXXFLAGS)

gnuobj/bandmat_penalty.o: bib/bandmat_penalty.cpp
	$(CPP) -c bib/bandmat_penalty.cpp -o gnuobj/bandmat_penalty.o $(CXXFLAGS)

gnuobj/bayesreg.o: bib/bayesreg.cpp
	$(CPP) -c bib/bayesreg.cpp -o gnuobj/bayesreg.o $(CXXFLAGS)

gnuobj/bayesreg2.o: bib/bayesreg2.cpp
	$(CPP) -c bib/bayesreg2.cpp -o gnuobj/bayesreg2.o $(CXXFLAGS)

gnuobj/bayesreg3.o: bib/bayesreg3.cpp
	$(CPP) -c bib/bayesreg3.cpp -o gnuobj/bayesreg3.o $(CXXFLAGS)

gnuobj/clstring.o: bib/clstring.cpp
	$(CPP) -c bib/clstring.cpp -o gnuobj/clstring.o $(CXXFLAGS)

gnuobj/command.o: bib/command.cpp
	$(CPP) -c bib/command.cpp -o gnuobj/command.o $(CXXFLAGS)

gnuobj/data.o: bib/data.cpp
	$(CPP) -c bib/data.cpp -o gnuobj/data.o $(CXXFLAGS)

gnuobj/dataobj.o: bib/dataobj.cpp
	$(CPP) -c bib/dataobj.cpp -o gnuobj/dataobj.o $(CXXFLAGS)

gnuobj/envmatrix.o: bib/envmatrix.cpp
	$(CPP) -c bib/envmatrix.cpp -o gnuobj/envmatrix.o $(CXXFLAGS)

gnuobj/envmatrix_penalty.o: bib/envmatrix_penalty.cpp
	$(CPP) -c bib/envmatrix_penalty.cpp -o gnuobj/envmatrix_penalty.o $(CXXFLAGS)

gnuobj/graph.o: bib/graph.cpp
	$(CPP) -c bib/graph.cpp -o gnuobj/graph.o $(CXXFLAGS)

gnuobj/map.o: bib/map.cpp
	$(CPP) -c bib/map.cpp -o gnuobj/map.o $(CXXFLAGS)

gnuobj/mapobject.o: bib/mapobject.cpp
	$(CPP) -c bib/mapobject.cpp -o gnuobj/mapobject.o $(CXXFLAGS)

gnuobj/model.o: bib/model.cpp
	$(CPP) -c bib/model.cpp -o gnuobj/model.o $(CXXFLAGS)

gnuobj/model_remlreg.o: bib/model_remlreg.cpp
	$(CPP) -c bib/model_remlreg.cpp -o gnuobj/model_remlreg.o $(CXXFLAGS)

gnuobj/model_stepwise.o: bib/model_stepwise.cpp
	$(CPP) -c bib/model_stepwise.cpp -o gnuobj/model_stepwise.o $(CXXFLAGS)

gnuobj/option.o: bib/option.cpp
	$(CPP) -c bib/option.cpp -o gnuobj/option.o $(CXXFLAGS)

gnuobj/realobs.o: bib/realobs.cpp
	$(CPP) -c bib/realobs.cpp -o gnuobj/realobs.o $(CXXFLAGS)

gnuobj/realvar.o: bib/realvar.cpp
	$(CPP) -c bib/realvar.cpp -o gnuobj/realvar.o $(CXXFLAGS)

gnuobj/remlreg.o: bib/remlreg.cpp
	$(CPP) -c bib/remlreg.cpp -o gnuobj/remlreg.o $(CXXFLAGS)

gnuobj/sparsemat.o: bib/sparsemat.cpp
	$(CPP) -c bib/sparsemat.cpp -o gnuobj/sparsemat.o $(CXXFLAGS)

gnuobj/statmat.o: bib/statmat.cpp
	$(CPP) -c bib/statmat.cpp -o gnuobj/statmat.o $(CXXFLAGS)

gnuobj/statmat_penalty.o: bib/statmat_penalty.cpp
	$(CPP) -c bib/statmat_penalty.cpp -o gnuobj/statmat_penalty.o $(CXXFLAGS)

gnuobj/statobj.o: bib/statobj.cpp
	$(CPP) -c bib/statobj.cpp -o gnuobj/statobj.o $(CXXFLAGS)

gnuobj/stepwisereg.o: bib/stepwisereg.cpp
	$(CPP) -c bib/stepwisereg.cpp -o gnuobj/stepwisereg.o $(CXXFLAGS)

gnuobj/use.o: bib/use.cpp
	$(CPP) -c bib/use.cpp -o gnuobj/use.o $(CXXFLAGS)

gnuobj/vectorn.o: bib/vectorn.cpp
	$(CPP) -c bib/vectorn.cpp -o gnuobj/vectorn.o $(CXXFLAGS)


gnuobj/adjacency.o: dag/adjacency.cpp
	$(CPP) -c dag/adjacency.cpp -o gnuobj/adjacency.o $(CXXFLAGS)

gnuobj/dagobject.o: dag/dagobject.cpp
	$(CPP) -c dag/dagobject.cpp -o gnuobj/dagobject.o $(CXXFLAGS)

gnuobj/fullcond_dag.o: dag/fullcond_dag.cpp
	$(CPP) -c dag/fullcond_dag.cpp -o gnuobj/fullcond_dag.o $(CXXFLAGS)

gnuobj/fullcond_dag_d.o: dag/fullcond_dag_d.cpp
	$(CPP) -c dag/fullcond_dag_d.cpp -o gnuobj/fullcond_dag_d.o $(CXXFLAGS)

gnuobj/fullcond_dag_ia.o: dag/fullcond_dag_ia.cpp
	$(CPP) -c dag/fullcond_dag_ia.cpp -o gnuobj/fullcond_dag_ia.o $(CXXFLAGS)

gnuobj/fullcond_dag_ia_mixed.o: dag/fullcond_dag_ia_mixed.cpp
	$(CPP) -c dag/fullcond_dag_ia_mixed.cpp -o gnuobj/fullcond_dag_ia_mixed.o $(CXXFLAGS)

gnuobj/fullcond_rj.o: dag/fullcond_rj.cpp
	$(CPP) -c dag/fullcond_rj.cpp -o gnuobj/fullcond_rj.o $(CXXFLAGS)

gnuobj/fullcond_rj_int.o: dag/fullcond_rj_int.cpp
	$(CPP) -c dag/fullcond_rj_int.cpp -o gnuobj/fullcond_rj_int.o $(CXXFLAGS)

gnuobj/fullcond_rj_mix.o: dag/fullcond_rj_mix.cpp
	$(CPP) -c dag/fullcond_rj_mix.cpp -o gnuobj/fullcond_rj_mix.o $(CXXFLAGS)

gnuobj/func_dag.o: dag/func_dag.cpp
	$(CPP) -c dag/func_dag.cpp -o gnuobj/func_dag.o $(CXXFLAGS)

gnuobj/ia.o: dag/ia.cpp
	$(CPP) -c dag/ia.cpp -o gnuobj/ia.o $(CXXFLAGS)

gnuobj/ia_mixed.o: dag/ia_mixed.cpp
	$(CPP) -c dag/ia_mixed.cpp -o gnuobj/ia_mixed.o $(CXXFLAGS)



gnuobj/nbinomial.o: leyre/nbinomial.cpp
	$(CPP) -c leyre/nbinomial.cpp -o gnuobj/nbinomial.o $(CXXFLAGS)

gnuobj/zip.o: leyre/zip.cpp
	$(CPP) -c leyre/zip.cpp -o gnuobj/zip.o $(CXXFLAGS)



gnuobj/baseline_reml.o: mcmc/baseline_reml.cpp
	$(CPP) -c mcmc/baseline_reml.cpp -o gnuobj/baseline_reml.o $(CXXFLAGS)

gnuobj/distribution.o: mcmc/distribution.cpp
	$(CPP) -c mcmc/distribution.cpp -o gnuobj/distribution.o $(CXXFLAGS)

gnuobj/fullcond.o: mcmc/fullcond.cpp
	$(CPP) -c mcmc/fullcond.cpp -o gnuobj/fullcond.o $(CXXFLAGS)

gnuobj/fullcond_merror.o: mcmc/fullcond_merror.cpp
	$(CPP) -c mcmc/fullcond_merror.cpp -o gnuobj/fullcond_merror.o $(CXXFLAGS)

gnuobj/fullcond_mult.o: mcmc/fullcond_mult.cpp
	$(CPP) -c mcmc/fullcond_mult.cpp -o gnuobj/fullcond_mult.o $(CXXFLAGS)

gnuobj/fullcond_nonp_gaussian.o: mcmc/fullcond_nonp_gaussian.cpp
	$(CPP) -c mcmc/fullcond_nonp_gaussian.cpp -o gnuobj/fullcond_nonp_gaussian.o $(CXXFLAGS)

gnuobj/fullcond_nonp_gaussian_stepwise.o: mcmc/fullcond_nonp_gaussian_stepwise.cpp
	$(CPP) -c mcmc/fullcond_nonp_gaussian_stepwise.cpp -o gnuobj/fullcond_nonp_gaussian_stepwise.o $(CXXFLAGS)

gnuobj/gaussian_heteroskedastic.o: mcmc/gaussian_heteroskedastic.cpp
	$(CPP) -c mcmc/gaussian_heteroskedastic.cpp -o gnuobj/gaussian_heteroskedastic.o $(CXXFLAGS)

gnuobj/hrandom.o: mcmc/hrandom.cpp
	$(CPP) -c mcmc/hrandom.cpp -o gnuobj/hrandom.o $(CXXFLAGS)

gnuobj/kriging2.o: mcmc/kriging2.cpp
	$(CPP) -c mcmc/kriging2.cpp -o gnuobj/kriging2.o $(CXXFLAGS)

gnuobj/kriging.o: mcmc/kriging.cpp
	$(CPP) -c mcmc/kriging.cpp -o gnuobj/kriging.o $(CXXFLAGS)

gnuobj/mcmc.o: mcmc/mcmc.cpp
	$(CPP) -c mcmc/mcmc.cpp -o gnuobj/mcmc.o $(CXXFLAGS)

gnuobj/mcmc_const.o: mcmc/mcmc_const.cpp
	$(CPP) -c mcmc/mcmc_const.cpp -o gnuobj/mcmc_const.o $(CXXFLAGS)

gnuobj/mcmc_const_stepwise.o: mcmc/mcmc_const_stepwise.cpp
	$(CPP) -c mcmc/mcmc_const_stepwise.cpp -o gnuobj/mcmc_const_stepwise.o $(CXXFLAGS)

gnuobj/mcmc_nonp.o: mcmc/mcmc_nonp.cpp
	$(CPP) -c mcmc/mcmc_nonp.cpp -o gnuobj/mcmc_nonp.o $(CXXFLAGS)

gnuobj/mcmc_nonpbasis.o: mcmc/mcmc_nonpbasis.cpp
	$(CPP) -c mcmc/mcmc_nonpbasis.cpp -o gnuobj/mcmc_nonpbasis.o $(CXXFLAGS)

gnuobj/mcmcsimul.o: mcmc/mcmcsimul.cpp
	$(CPP) -c mcmc/mcmcsimul.cpp -o gnuobj/mcmcsimul.o $(CXXFLAGS)

gnuobj/mcmcsimul2.o: mcmc/mcmcsimul2.cpp
	$(CPP) -c mcmc/mcmcsimul2.cpp -o gnuobj/mcmcsimul2.o $(CXXFLAGS)

gnuobj/mcmcsimul2_multi.o: mcmc/mcmcsimul2_multi.cpp
	$(CPP) -c mcmc/mcmcsimul2_multi.cpp -o gnuobj/mcmcsimul2_multi.o $(CXXFLAGS)

gnuobj/randomeffect.o: mcmc/randomeffect.cpp
	$(CPP) -c mcmc/randomeffect.cpp -o gnuobj/randomeffect.o $(CXXFLAGS)

gnuobj/randomeffect_stepwise.o: mcmc/randomeffect_stepwise.cpp
	$(CPP) -c mcmc/randomeffect_stepwise.cpp -o gnuobj/randomeffect_stepwise.o $(CXXFLAGS)

gnuobj/remlest.o: mcmc/remlest.cpp
	$(CPP) -c mcmc/remlest.cpp -o gnuobj/remlest.o $(CXXFLAGS)

gnuobj/remlest_multi.o: mcmc/remlest_multi.cpp
	$(CPP) -c mcmc/remlest_multi.cpp -o gnuobj/remlest_multi.o $(CXXFLAGS)

gnuobj/remlest_multi2.o: mcmc/remlest_multi2.cpp
	$(CPP) -c mcmc/remlest_multi2.cpp -o gnuobj/remlest_multi2.o $(CXXFLAGS)

gnuobj/remlest_multi3.o: mcmc/remlest_multi3.cpp
	$(CPP) -c mcmc/remlest_multi3.cpp -o gnuobj/remlest_multi3.o $(CXXFLAGS)

gnuobj/tvariance.o: mcmc/tvariance.cpp
	$(CPP) -c mcmc/tvariance.cpp -o gnuobj/tvariance.o $(CXXFLAGS)

gnuobj/tvariance2dim.o: mcmc/tvariance2dim.cpp
	$(CPP) -c mcmc/tvariance2dim.cpp -o gnuobj/tvariance2dim.o $(CXXFLAGS)

gnuobj/variance_nonp.o: mcmc/variance_nonp.cpp
	$(CPP) -c mcmc/variance_nonp.cpp -o gnuobj/variance_nonp.o $(CXXFLAGS)

gnuobj/variance_nonp_vector.o: mcmc/variance_nonp_vector.cpp
	$(CPP) -c mcmc/variance_nonp_vector.cpp -o gnuobj/variance_nonp_vector.o $(CXXFLAGS)

gnuobj/variance_nonp_vector_nigmix.o: mcmc/variance_nonp_vector_nigmix.cpp
	$(CPP) -c mcmc/variance_nonp_vector_nigmix.cpp -o gnuobj/variance_nonp_vector_nigmix.o $(CXXFLAGS)



gnuobj/IWLS_pspline.o: psplines/IWLS_pspline.cpp
	$(CPP) -c psplines/IWLS_pspline.cpp -o gnuobj/IWLS_pspline.o $(CXXFLAGS)

gnuobj/bsplinemat.o: psplines/bsplinemat.cpp
	$(CPP) -c psplines/bsplinemat.cpp -o gnuobj/bsplinemat.o $(CXXFLAGS)

gnuobj/fullcond_pspline_gaussian.o: psplines/fullcond_pspline_gaussian.cpp
	$(CPP) -c psplines/fullcond_pspline_gaussian.cpp -o gnuobj/fullcond_pspline_gaussian.o $(CXXFLAGS)

gnuobj/fullcond_pspline_stepwise.o: psplines/fullcond_pspline_stepwise.cpp
	$(CPP) -c psplines/fullcond_pspline_stepwise.cpp -o gnuobj/fullcond_pspline_stepwise.o $(CXXFLAGS)

gnuobj/fullcond_pspline_surf_gaussian.o: psplines/fullcond_pspline_surf_gaussian.cpp
	$(CPP) -c psplines/fullcond_pspline_surf_gaussian.cpp -o gnuobj/fullcond_pspline_surf_gaussian.o $(CXXFLAGS)

gnuobj/fullcond_pspline_surf_stepwise.o: psplines/fullcond_pspline_surf_stepwise.cpp
	$(CPP) -c psplines/fullcond_pspline_surf_stepwise.cpp -o gnuobj/fullcond_pspline_surf_stepwise.o $(CXXFLAGS)

gnuobj/mcmc_pspline.o: psplines/mcmc_pspline.cpp
	$(CPP) -c psplines/mcmc_pspline.cpp -o gnuobj/mcmc_pspline.o $(CXXFLAGS)

gnuobj/mcmc_pspline_surf.o: psplines/mcmc_pspline_surf.cpp
	$(CPP) -c psplines/mcmc_pspline_surf.cpp -o gnuobj/mcmc_pspline_surf.o $(CXXFLAGS)

gnuobj/spline_basis.o: psplines/spline_basis.cpp
	$(CPP) -c psplines/spline_basis.cpp -o gnuobj/spline_basis.o $(CXXFLAGS)

gnuobj/spline_basis_surf.o: psplines/spline_basis_surf.cpp
	$(CPP) -c psplines/spline_basis_surf.cpp -o gnuobj/spline_basis_surf.o $(CXXFLAGS)



gnuobj/multgaussian.o: samson/multgaussian.cpp
	$(CPP) -c samson/multgaussian.cpp -o gnuobj/multgaussian.o $(CXXFLAGS)



gnuobj/FC.o: structadd/FC.cpp
	$(CPP) -c structadd/FC.cpp -o gnuobj/FC.o $(CXXFLAGS)

gnuobj/FC_hrandom.o: structadd/FC_hrandom.cpp
	$(CPP) -c structadd/FC_hrandom.cpp -o gnuobj/FC_hrandom.o $(CXXFLAGS)

gnuobj/FC_hrandom_variance.o: structadd/FC_hrandom_variance.cpp
	$(CPP) -c structadd/FC_hrandom_variance.cpp -o gnuobj/FC_hrandom_variance.o $(CXXFLAGS)

gnuobj/FC_linear.o: structadd/FC_linear.cpp
	$(CPP) -c structadd/FC_linear.cpp -o gnuobj/FC_linear.o $(CXXFLAGS)

gnuobj/FC_mult.o: structadd/FC_mult.cpp
	$(CPP) -c structadd/FC_mult.cpp -o gnuobj/FC_mult.o $(CXXFLAGS)

gnuobj/FC_nonp.o: structadd/FC_nonp.cpp
	$(CPP) -c structadd/FC_nonp.cpp -o gnuobj/FC_nonp.o $(CXXFLAGS)

gnuobj/FC_nonp_variance.o: structadd/FC_nonp_variance.cpp
	$(CPP) -c structadd/FC_nonp_variance.cpp -o gnuobj/FC_nonp_variance.o $(CXXFLAGS)

gnuobj/FC_predict.o: structadd/FC_predict.cpp
	$(CPP) -c structadd/FC_predict.cpp -o gnuobj/FC_predict.o $(CXXFLAGS)

gnuobj/GENERAL_OPTIONS.o: structadd/GENERAL_OPTIONS.cpp
	$(CPP) -c structadd/GENERAL_OPTIONS.cpp -o gnuobj/GENERAL_OPTIONS.o $(CXXFLAGS)

gnuobj/MASTER_obj.o: structadd/MASTER_obj.cpp
	$(CPP) -c structadd/MASTER_obj.cpp -o gnuobj/MASTER_obj.o $(CXXFLAGS)

gnuobj/design.o: structadd/design.cpp
	$(CPP) -c structadd/design.cpp -o gnuobj/design.o $(CXXFLAGS)

gnuobj/design_hrandom.o: structadd/design_hrandom.cpp
	$(CPP) -c structadd/design_hrandom.cpp -o gnuobj/design_hrandom.o $(CXXFLAGS)

gnuobj/design_mrf.o: structadd/design_mrf.cpp
	$(CPP) -c structadd/design_mrf.cpp -o gnuobj/design_mrf.o $(CXXFLAGS)

gnuobj/design_pspline.o: structadd/design_pspline.cpp
	$(CPP) -c structadd/design_pspline.cpp -o gnuobj/design_pspline.o $(CXXFLAGS)

gnuobj/distr.o: structadd/distr.cpp
	$(CPP) -c structadd/distr.cpp -o gnuobj/distr.o $(CXXFLAGS)

gnuobj/distr_categorical.o: structadd/distr_categorical.cpp
	$(CPP) -c structadd/distr_categorical.cpp -o gnuobj/distr_categorical.o $(CXXFLAGS)

gnuobj/mcmcsim.o: structadd/mcmcsim.cpp
	$(CPP) -c structadd/mcmcsim.cpp -o gnuobj/mcmcsim.o $(CXXFLAGS)

gnuobj/model_parameters.o: structadd/model_parameters.cpp
	$(CPP) -c structadd/model_parameters.cpp -o gnuobj/model_parameters.o $(CXXFLAGS)

gnuobj/superbayesreg.o: structadd/superbayesreg.cpp
	$(CPP) -c structadd/superbayesreg.cpp -o gnuobj/superbayesreg.o $(CXXFLAGS)

gnuobj/design_kriging.o: structadd/design_kriging.cpp
	$(CPP) -c structadd/design_kriging.cpp -o gnuobj/design_kriging.o $(CXXFLAGS)

gnuobj/FC_predictive_check.o: structadd/FC_predictive_check.cpp
	$(CPP) -c structadd/FC_predictive_check.cpp -o gnuobj/FC_predictive_check.o $(CXXFLAGS)

gnuobj/main.o: main.cpp
	$(CPP) -c main.cpp -o gnuobj/main.o $(CXXFLAGS)


