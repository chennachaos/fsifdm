#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineBase.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"


extern DomainTree domain;


int macro230(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "fsi";
    macro.type = "chen";
    macro.what = " FSI schemes ";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*StagForcePred","StagDispPred","MonoFixedPointFP","MonoFixedPointDP");

    macro.db.frameButtonBox();

    macro.db.addTextField("Steps = ",2);

    macro.db.addTextField("Tol = ",0.0001,6);

    macro.db.frameRadioBox();

    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

  int  domType=27, id=0, Steps, soltype;
  double  tol=0.01;

  //domType  = roundToInt(macro.p[0]);
  //id       = roundToInt(macro.p[1]) - 1;
  //soltype  = roundToInt(macro.p[2]);
  Steps    = roundToInt(macro.p[2]);

  //cout << " domType = " << domType << endl;
  //cout << " id = " << id << endl;
  //cout << " soltype = " << soltype << endl;
  //cout << " Steps = " << Steps << endl;
  //cout << macro.p[0] << '\t' << macro.p[1] << '\t' << macro.p[2] << '\t' << macro.p[3] << '\t' << macro.p[4] << endl;

  if(domType == 27)
  {
    hbsplineFEM(domain(domType,id)).solve_fsi(Steps, tol);
    /*
    if(soltype == 111)
      hbsplineFEM(domain(domType,id)).fsi_monolithic(Steps, tol);
    else if(soltype == 1)
      hbsplineFEM(domain(domType,id)).fsi_staggered_force_predictor(Steps, tol);
    else if(soltype == 2)
      hbsplineFEM(domain(domType,id)).fsi_staggered_displacement_predictor(Steps, tol);
    else if(soltype == 3)
      hbsplineFEM(domain(domType,id)).fsi_monolithic_fixedpoint_forcePred(Steps, tol);
    else if(soltype == 4)
      hbsplineFEM(domain(domType,id)).fsi_monolithic_fixedpoint_dispPred(Steps, tol);
    else
    {
      cout << " macro230 ... FSI schemes ... undefined for this domain " << endl;
    }
    */
  }
  else if(domType == 28)
  {
    if(soltype == 1)
      hbsplineCutFEM(domain(domType,id)).fsi_staggered_force_predictor(Steps, tol);
    else if(soltype == 2)
      hbsplineCutFEM(domain(domType,id)).fsi_staggered_displacement_predictor(Steps, tol);
    else if(soltype == 3)
      hbsplineCutFEM(domain(domType,id)).fsi_monolithic_fixedpoint_forcePred(Steps, tol);
    else if(soltype == 4)
      hbsplineCutFEM(domain(domType,id)).fsi_monolithic_fixedpoint_dispPred(Steps, tol);
    else
    {
      cout << " macro230 ... FSI schemes ... undefined for this domain " << endl;
    }
  }

  //--------------------------------------------------------------------------------------------------

  return 0;  
}

