
g++ -Wno-deprecated -c -o gnuobj/clstring.o bib/clstring.cpp

g++ -Wno-deprecated -c -o gnuobj/tlinklst.o bib/tlinklst.cpp
g++ -Wno-deprecated -c -o gnuobj/tarray.o bib/tarray.cpp
g++ -Wno-deprecated -c -o gnuobj/tarray2d.o bib/tarray2d.cpp
g++ -Wno-deprecated -c -o gnuobj/tpremat.o bib/tpremat.cpp
g++ -Wno-deprecated -c -o gnuobj/tmatrix.o bib/tmatrix.cpp
g++ -Wno-deprecated -c -o gnuobj/random.o bib/random.cpp

g++ -w -Wno-deprecated -c -o gnuobj/graph.o bib/graph.cpp
g++ -Wno-deprecated -c -o gnuobj/adminparse_basic.o bib/adminparse_basic.cpp
g++ -w -Wno-deprecated -c -o gnuobj/map.o bib/map.cpp

g++ -Wno-deprecated -c -o gnuobj/sparsemat.o bib/sparsemat.cpp
g++ -Wno-deprecated -c -o gnuobj/statmat.o bib/statmat.cpp
g++ -Wno-deprecated -c -o gnuobj/statmat_penalty.o bib/statmat_penalty.cpp

g++ -Wno-deprecated -c -o gnuobj/realobs.o bib/realobs.cpp
g++ -Wno-deprecated -c -o gnuobj/vectorn.o bib/vectorn.cpp
g++ -Wno-deprecated -c -o gnuobj/realvar.o bib/realvar.cpp
g++ -w -Wno-deprecated -c -o gnuobj/data.o bib/data.cpp
g++ -Wno-deprecated -c -o gnuobj/option.o bib/option.cpp
g++ -Wno-deprecated -c -o gnuobj/model.o bib/model.cpp
g++ -Wno-deprecated -c -o gnuobj/use.o bib/use.cpp
g++ -Wno-deprecated -c -o gnuobj/command.o bib/command.cpp

g++ -Wno-deprecated -c -o gnuobj/adminparse_pointers.o bib/adminparse_pointers.cpp
g++ -Wno-deprecated -c -o gnuobj/statobj.o bib/statobj.cpp

g++ -w -Wno-deprecated -c -o gnuobj/dataobj.o bib/dataobj.cpp

g++ -Wno-deprecated -c -o gnuobj/mapobject.o bib/mapobject.cpp

g++ -Wno-deprecated -c -Ibib -o gnuobj/graphobj.o graph/graphobj.cpp

g++ -Wno-deprecated -c -Ibib -o gnuobj/mcmc.o mcmc/mcmc.cpp
g++ -w -Wno-deprecated -c -Ibib -o gnuobj/fullcond.o mcmc/fullcond.cpp
g++ -w -Wno-deprecated -c -Ibib -o gnuobj/distribution.o mcmc/distribution.cpp

g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/nbinomial.o leyre/nbinomial.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/mcmc_const.o mcmc/mcmc_const.cpp
g++ -w -Wno-deprecated -c -o gnuobj/bandmat.o bib/bandmat.cpp
g++ -w -Wno-deprecated -c -o gnuobj/envmatrix.o bib/envmatrix.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/mcmc_nonpbasis.o mcmc/mcmc_nonpbasis.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/mcmc_nonp.o mcmc/mcmc_nonp.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/randomeffect.o mcmc/randomeffect.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/fullcond_nonp_gaussian.o mcmc/fullcond_nonp_gaussian.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/tvariance.o mcmc/tvariance.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/variance_nonp.o mcmc/variance_nonp.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/fullcond_nonp_gaussian_stepwise.o mcmc/fullcond_nonp_gaussian_stepwise.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/randomeffect_stepwise.o mcmc/randomeffect_stepwise.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/bsplinemat.o psplines/bsplinemat.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/spline_basis.o psplines/spline_basis.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/spline_basis_surf.o psplines/spline_basis_surf.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/fullcond_pspline_gaussian.o psplines/fullcond_pspline_gaussian.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/fullcond_pspline_surf_gaussian.o psplines/fullcond_pspline_surf_gaussian.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/IWLS_pspline.o psplines/IWLS_pspline.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/mcmc_pspline.o psplines/mcmc_pspline.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/mcmc_pspline_surf.o psplines/mcmc_pspline_surf.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -Isamson -o gnuobj/multgaussian.o samson/multgaussian.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -o gnuobj/fullcond_adaptiv.o adaptiv/fullcond_adaptiv.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/tvariance2dim.o mcmc/tvariance2dim.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/cox.o andrea/cox.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/baseline.o andrea/baseline.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/multibaseline.o andrea/multibaseline.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -Iandrea -o gnuobj/mcmcsimul.o mcmc/mcmcsimul.cpp

g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -Iandrea -Isamson -Ialex -Iadaptiv -o gnuobj/bayesreg.o bib/bayesreg.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -Iandrea -Isamson -Ialex -Iadaptiv -o gnuobj/bayesreg2.o bib/bayesreg2.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -Iandrea -Isamson -Ialex -Iadaptiv -o gnuobj/bayesreg3.o bib/bayesreg3.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/zip.o leyre/zip.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/kriging2.o mcmc/kriging2.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/mixture.o alex/mixture.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/IWLS_baseline.o andrea/IWLS_baseline.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/multistate.o andrea/multistate.cpp

g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -o gnuobj/model_remlreg.o bib/model_remlreg.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -o gnuobj/kriging.o mcmc/kriging.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -o gnuobj/baseline_reml.o mcmc/baseline_reml.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -o gnuobj/remlest.o mcmc/remlest.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -o gnuobj/remlest_multi.o mcmc/remlest_multi.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -o gnuobj/remlest_multi2.o mcmc/remlest_multi2.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -o gnuobj/remlest_multi3.o mcmc/remlest_multi3.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -o gnuobj/remlest_cox.o mcmc/remlest_cox.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -o gnuobj/remlest_multistate.o mcmc/remlest_multistate.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -o gnuobj/remlreg.o bib/remlreg.cpp

g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -o gnuobj/model_stepwise.o bib/model_stepwise.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -o gnuobj/mcmc_const_stepwise.o mcmc/mcmc_const_stepwise.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -Iandrea -o gnuobj/mcmcsimul2.o mcmc/mcmcsimul2.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -Iandrea -o gnuobj/mcmcsimul2_multi.o mcmc/mcmcsimul2_multi.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/fullcond_pspline_stepwise.o psplines/fullcond_pspline_stepwise.cpp
g++ -w -Wno-deprecated -c -Ibib -Ileyre -Imcmc -Ipsplines -o gnuobj/fullcond_pspline_surf_stepwise.o psplines/fullcond_pspline_surf_stepwise.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -Ipsplines -Ileyre -Iandrea -o gnuobj/stepwisereg.o bib/stepwisereg.cpp

g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/adjacency.o dag/adjacency.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/fullcond_dag.o dag/fullcond_dag.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/fullcond_dag_d.o dag/fullcond_dag_d.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/func_dag.o dag/func_dag.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/ia.o dag/ia.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/fullcond_dag_ia.o dag/fullcond_dag_ia.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/fullcond_dag_ia_mixed.o dag/fullcond_dag_ia_mixed.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/ia_ok.o dag/ia_ok.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/ia_mixed.o dag/ia_mixed.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/functions.o dag/functions.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/fullcond_rj_mix.o dag/fullcond_rj_mix.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/fullcond_rj.o dag/fullcond_rj.cpp
g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/fullcond_rj_int.o dag/fullcond_rj_int.cpp
//g++ -w -Wno-deprecated -c -Ibib -Imcmc -o gnuobj/fullcond_rj_ia.o dag/fullcond_rj_ia.cpp

g++ -w -Wno-deprecated -c -Ibib -Imcmc -Iandrea -Ipsplines -Ileyre -o gnuobj/dagobject.o dag/dagobject.cpp

g++ -w -Wno-deprecated -c -Ibib -Imcmc -Igraph -Iandrea -Ipsplines -Ialex -Isamson -Iadaptiv -Idag -Ileyre -o gnuobj/adminparse.o bib/adminparse.cpp

