//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
USEFORM("statwinwindows\StatwinFrame.cpp", Frame);
USEFORM("statwinwindows\StatReview.cpp", Review);
USEFORM("statwinwindows\StatResults.cpp", Results);
USEUNIT("mcmc\randomeffect.cpp");
USEUNIT("mcmc\fullcond.cpp");
USEUNIT("mcmc\mcmc.cpp");
USEUNIT("mcmc\mcmc_const.cpp");
USEUNIT("mcmc\mcmc_nonp.cpp");
USEUNIT("mcmc\mcmcsimul.cpp");
USEUNIT("mcmc\distribution.cpp");
USEFORM("bib\describe_map.cpp", mapform);
USEUNIT("bib\data.cpp");
USEUNIT("bib\dataobj.cpp");
USEFORM("bib\describe_dataset.cpp", datasetform);
USEUNIT("bib\clstring.cpp");
USEUNIT("bib\errorm.cpp");
USEUNIT("bib\map.cpp");
USEUNIT("bib\mapobject.cpp");
USEUNIT("bib\realobs.cpp");
USEUNIT("bib\realvar.cpp");
USEUNIT("bib\remlreg.cpp");
USEUNIT("bib\sparsemat.cpp");
USEUNIT("bib\statobj.cpp");
USEUNIT("dag\dagobject.cpp");
USEUNIT("dag\fullcond_dag.cpp");
USEUNIT("dag\fullcond_rj.cpp");
USEUNIT("mcmc\fullcond_nonp_gaussian.cpp");
USEUNIT("mcmc\mcmc_nonpbasis.cpp");
USEUNIT("mcmc\variance_nonp.cpp");
USEUNIT("..\cprog\adaptiv\fullcond_adaptiv.cpp");
USEUNIT("bib\graph.cpp");
USEUNIT("dag\fullcond_rj_int.cpp");
USEUNIT("dag\fullcond_dag_ia.cpp");
USEUNIT("dag\fullcond_dag_d.cpp");
USEUNIT("dag\func_dag.cpp");
USEUNIT("dag\ia.cpp");
USEUNIT("mcmc\tvariance2dim.cpp");
USEUNIT("mcmc\tvariance.cpp");
USEUNIT("psplines\spline_basis_surf.cpp");
USEUNIT("psplines\IWLS_pspline.cpp");
USEUNIT("dag\ia_ok.cpp");
USEUNIT("dag\fullcond_rj_mix.cpp");
USEUNIT("dag\fullcond_dag_ia_mixed.cpp");
USEUNIT("dag\ia_mixed.cpp");
USEUNIT("andrea\cox.cpp");
USEUNIT("andrea\baseline.cpp");
USEUNIT("psplines\bsplinemat.cpp");
USEUNIT("leyre\zip.cpp");
USEUNIT("leyre\nbinomial.cpp");
USEUNIT("bib\bayesreg3.cpp");
USEUNIT("bib\bayesreg2.cpp");
USEUNIT("bib\bayesreg.cpp");
USEUNIT("mcmc\mcmcsimul2.cpp");
USEUNIT("samson\multgaussian.cpp");
USEUNIT("bib\model_remlreg.cpp");
USEUNIT("bib\model.cpp");
USEUNIT("bib\command.cpp");
USEUNIT("bib\model_stepwise.cpp");
USEUNIT("bib\option.cpp");
USEUNIT("bib\Random.cpp");
USEUNIT("bib\stepwisereg.cpp");
USEUNIT("bib\use.cpp");
USEUNIT("mcmc\remlest_multi2.cpp");
USEUNIT("mcmc\kriging.cpp");
USEUNIT("mcmc\kriging2.cpp");
USEUNIT("mcmc\remlest_multi.cpp");
USEUNIT("mcmc\baseline_reml.cpp");
USEFORM("statwinwindows\StatObjects.cpp", Objectbrowser);
USEUNIT("psplines\fullcond_pspline_gaussian.cpp");
USEUNIT("psplines\fullcond_pspline_surf_gaussian.cpp");
USEUNIT("psplines\mcmc_pspline.cpp");
USEUNIT("psplines\mcmc_pspline_surf.cpp");
USEUNIT("psplines\spline_basis.cpp");
USEUNIT("dag\adjacency.cpp");
USEFORM("statwinwindows\statwin_haupt.cpp", hauptformular);
USEUNIT("bib\statmat_penalty.cpp");
USEUNIT("mcmc\remlest.cpp");
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
    try
    {
        Application->Initialize();
        Application->Title = "BayesX";
        Application->CreateForm(__classid(TFrame), &Frame);
        Application->CreateForm(__classid(TReview), &Review);
        Application->CreateForm(__classid(TResults), &Results);
        Application->CreateForm(__classid(Tmapform), &mapform);
        Application->CreateForm(__classid(Tdatasetform), &datasetform);
        Application->CreateForm(__classid(TObjectbrowser), &Objectbrowser);
        Application->CreateForm(__classid(Thauptformular), &hauptformular);
        Application->Run();
    }
    catch (Exception &exception)
    {
        Application->ShowException(&exception);
    }
    return 0;
}
//---------------------------------------------------------------------------
