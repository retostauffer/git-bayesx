# Project: BayesX

ANDREA_SRC = \
	andrea/baseline.cpp\
	andrea/cox.cpp\
	andrea/IWLS_baseline.cpp\
	andrea/multibaseline.cpp\
	andrea/multistate.cpp
BIB_SRC = \
	bib/Random.cpp\
	bib/adminparse_gnu.cpp\
	bib/bandmat.cpp\
	bib/bandmat_penalty.cpp\
	bib/bayesreg.cpp\
	bib/bayesreg2.cpp\
	bib/bayesreg3.cpp\
	bib/clstring.cpp\
	bib/command.cpp\
	bib/data.cpp\
	bib/dataobj.cpp\
	bib/envmatrix.cpp\
	bib/envmatrix_penalty.cpp\
	bib/graph.cpp\
	bib/map.cpp\
	bib/mapobject.cpp\
	bib/model.cpp\
	bib/model_remlreg.cpp\
	bib/model_stepwise.cpp\
	bib/option.cpp\
	bib/realobs.cpp\
	bib/realvar.cpp\
	bib/remlreg.cpp\
	bib/sparsemat.cpp\
	bib/statmat.cpp\
	bib/statmat_penalty.cpp\
	bib/statobj.cpp\
	bib/stepwisereg.cpp\
	bib/use.cpp\
	bib/vectorn.cpp
DAG_SRC = \
	dag/adjacency.cpp\
	dag/dagobject.cpp\
	dag/fullcond_dag.cpp\
	dag/fullcond_dag_d.cpp\
	dag/fullcond_dag_ia.cpp\
	dag/fullcond_dag_ia_mixed.cpp\
	dag/fullcond_rj.cpp\
	dag/fullcond_rj_int.cpp\
	dag/fullcond_rj_mix.cpp\
	dag/func_dag.cpp\
	dag/ia.cpp\
	dag/ia_mixed.cpp
LEYRE_SRC = \
	leyre/nbinomial.cpp\
	leyre/zip.cpp
MCMC_SRC = \
	mcmc/baseline_reml.cpp\
	mcmc/distribution.cpp\
	mcmc/fullcond.cpp\
	mcmc/fullcond_merror.cpp\
	mcmc/fullcond_mult.cpp\
	mcmc/fullcond_nonp_gaussian.cpp\
	mcmc/fullcond_nonp_gaussian_stepwise.cpp\
	mcmc/gaussian_heteroskedastic.cpp\
	mcmc/hrandom.cpp\
	mcmc/kriging2.cpp\
	mcmc/kriging.cpp\
	mcmc/mcmc.cpp\
	mcmc/mcmc_const.cpp\
	mcmc/mcmc_const_stepwise.cpp\
	mcmc/mcmc_nonp.cpp\
	mcmc/mcmc_nonpbasis.cpp\
	mcmc/mcmcsimul.cpp\
	mcmc/mcmcsimul2.cpp\
	mcmc/mcmcsimul2_multi.cpp\
	mcmc/randomeffect.cpp\
	mcmc/randomeffect_stepwise.cpp\
	mcmc/remlest.cpp\
	mcmc/remlest_multi.cpp\
	mcmc/remlest_multi2.cpp\
	mcmc/remlest_multi3.cpp\
	mcmc/tvariance.cpp\
	mcmc/tvariance2dim.cpp\
	mcmc/variance_nonp.cpp\
	mcmc/variance_nonp_vector.cpp\
	mcmc/variance_nonp_vector_nigmix.cpp
PSPLINES_SRC = \
	psplines/IWLS_pspline.cpp\
	psplines/bsplinemat.cpp\
	psplines/fullcond_pspline_gaussian.cpp\
	psplines/fullcond_pspline_stepwise.cpp\
	psplines/fullcond_pspline_surf_gaussian.cpp\
	psplines/fullcond_pspline_surf_stepwise.cpp\
	psplines/mcmc_pspline.cpp\
	psplines/mcmc_pspline_surf.cpp\
	psplines/spline_basis.cpp\
	psplines/spline_basis_surf.cpp
STRUCTADD_SRC = \
	structadd/FC.cpp\
	structadd/FC_hrandom.cpp\
	structadd/FC_hrandom_variance.cpp\
	structadd/FC_hrandom_variance_vec.cpp\
	structadd/FC_linear.cpp\
	structadd/FC_mult.cpp\
	structadd/FC_nonp.cpp\
	structadd/FC_nonp_variance.cpp\
	structadd/FC_predict.cpp\
	structadd/FC_cv.cpp\
	structadd/FC_variance_pen_vector.cpp\
	structadd/GENERAL_OPTIONS.cpp\
	structadd/MASTER_obj.cpp\
	structadd/design.cpp\
	structadd/design_hrandom.cpp\
	structadd/design_mrf.cpp\
	structadd/design_pspline.cpp\
	structadd/distr.cpp\
	structadd/distr_categorical.cpp\
	structadd/distr_categorical_mult.cpp\
	structadd/distr_mixture.cpp\
	structadd/mcmcsim.cpp\
	structadd/model_parameters.cpp\
	structadd/superbayesreg.cpp\
	structadd/design_kriging.cpp\
	structadd/FC_predictive_check.cpp\
	structadd/FC_predict_predictor.cpp\
	structadd/FC_nonp_variance_vec.cpp\
	structadd/distr_gamlss.cpp\
    structadd/distr_gamlss_nadja.cpp
SRC = \
	$(ANDREA_SRC)\
	$(BIB_SRC)\
	$(DAG_SRC)\
	$(LEYRE_SRC)\
	$(MCMC_SRC)\
	$(PSPLINES_SRC)\
	$(STRUCTADD_SRC)\
	main.cpp\
       	samson/multgaussian.cpp\
	adaptiv/fullcond_adaptiv.cpp\
	alex/mixture.cpp

OBJ = $(patsubst %.cpp,%.o,$(SRC))

INCS =  -I. -I"bib"  -I"alex"  -I"adaptiv"  -I"andrea"  -I"dag"  -I"graph"  -I"mcmc"  -I"psplines"  -I"samson"  -I"leyre"  -I"structadd"
CXXINCS = $(INCS)

WARNINGS = -Wcomment -Wdeprecated -Woverflow -Wpointer-arith\
           -Wuninitialized -Wunused-function\
           -Wwrite-strings -Wunknown-pragmas
           
ifeq ($(MAKECMDGOALS),intel)
    CPP  = icpc 
    CC   = icc 

    # icc warning flags yet to be added
    # -Wcast-qual
    # -Wconversion
    # -Wformat
    # -Wmissing-declarations
    # -Wmissing-prototypes
    # -Woverloaded-virtual
    # -Wp64
    # -Wreturn-type
    # -Wshadow
    # -Wshorten-64-to-32
    # -Wstrict-prototypes
    WARNINGS := $(WARNINGS)
else
    CPP  = g++ 
    CC   = gcc 

    # gcc warning flags still yet to be added
    # -Wcast-qual
    # -Wconversion
    # -Wformat
    # -Wmissing-declarations
    # -Wmissing-prototypes
    # -Woverloaded-virtual
    # -Wparentheses
    # -Wpointer-to-int-cast
    # -Wredundant-decls
    # -Wreturn-type
    # -Wshadow
    # -Wsign-compare
    # -Wsign-promo
    # -Wtype-limits
    # -Wunreachable-code
    # -Wunused-value
    # -Wall  -> and then reduce the number of redundant warnings
    WARNINGS := $(WARNINGS) -Waddress -Wendif-labels\
       	-Wmissing-braces -Wmissing-declarations -Wpragmas -Wundef -Wunused-macros\
       	-Wunused-variable 
    # -Wempty-body 
    # -Wtype-limits
endif


BIN  = BayesX
SYSTEM := $(shell uname -s)
DEFINES = -D__BUILDING_GNU -DTEMPL_INCL_DEF -D_MSC_VER2 -DNO_TEMPLATE_FRIENDS -DINCLUDE_REML -DINCLUDE_MCMC
ifeq ($(SYSTEM),Linux)
    DEFINES := $(DEFINES) -D__BUILDING_LINUX
    LIBS := ${LIBS} -lreadline
endif
ifeq ($(SYSTEM),CYGWIN_NT-5.1)
    DEFINES := $(DEFINES) -D__BUILDING_LINUX
    LIBS := ${LIBS} -lreadline
endif
ifeq ($(SYSTEM),Darwin)
    DEFINES := $(DEFINES) -D__BUILDING_LINUX
    LIBS := ${LIBS} -lreadline
    # Building Universal Binaries..
    #TARGET_ARCH = -arch i386 -arch x86_64
endif
ifeq ($(SYSTEM),FreeBSD)
    DEFINES := $(DEFINES) -D__BUILDING_LINUX
    LIBS := ${LIBS} -lreadline
endif
ifeq ($(SYSTEM),OpenBSD)
    DEFINES := $(DEFINES) -D__BUILDING_LINUX
    LIBS := ${LIBS} -lreadline -lncurses
endif
CFLAGS = $(WARNINGS) $(INCS) $(DEFINES) -O3 -ansi
CXXFLAGS = $(CFLAGS)
RM = rm -f

.PHONY: all clean clean-custom

all: $(BIN)

clean: clean-custom
	${RM} $(OBJ) $(BIN) $(subst .cpp,.d,$(SRC))

$(BIN): $(OBJ)
	$(CPP) $(OBJ) -o "BayesX" $(LIBS) $(TARGET_ARCH) 

# automatic dependency generation
ifneq ($(MAKECMDGOALS),clean)
    -include $(subst .cpp,.d,$(SRC))
endif

%.d: %.cpp
	$(CPP) -M $(CXXFLAGS) $< > $@.$$$$; \
	    base=`basename $*`; \
	    sed "s,$$base\.o[ :]*,$*.o : ,g" < $@.$$$$ > $@; \
	    rm -f $@.$$$$
