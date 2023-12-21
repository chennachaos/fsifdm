

#include "HBSplineFEM.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "myDataIntegrateCutFEM.h"
#include "ContactElementPointToPoint2D.h"
#include "QuadratureUtil.h"
#include "headersEigen.h"


#include <omp.h>


extern ComputerTime       computerTime;
extern MpapTime mpapTime;


void HBSplineFEM::setInitialConditions()
{
    double* tmp;

    solverEigen->zeroMtx();
/*
    for(int ii=0;ii<activeElements.size();ii++)
    {
       node *nd = elems[activeElements[ii]];

       nd->setInitialProfile();

       //cout << " AAAAAAAAAAAAAAAAA " << endl;

       //cout << " BBBBBBBBBBBBBBBBB " << endl;
       //printMatrix(Klocal);
       //printf("\n\n");
       //printVector(Flocal);

      //nd->assembleElementMatrix2(0, solverEigen->mtx);
      //nd->assembleElementVector(0, 0, &(solverEigen->rhsVec(0)));
    }

    solverEigen->currentStatus = ASSEMBLY_OK;

    factoriseSolveAndUpdate();
*/
    //
    int resln1[3]; resln1[0] = resln1[1] = 5;

    //postProcess2D(1, 1, 10, 1, 0.0, 1.0, resln1);
    //postProcess1D(1, 1, 10, 1, 0.0, 1.0, resln1);
    //

    //printVector(solver->soln);
    
    //solnInit = soln;
    //SolnData.var1 = solnInit;
    //SolnData.var1Prev  = SolnData.var1;
    //SolnData.var1Prev2 = SolnData.var1;
    //SolnData.var1Prev3 = SolnData.var1;
    //SolnData.var1Prev4 = SolnData.var1;

    //printVector(SolnData.var1);

   return;
}




int HBSplineFEM::calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes)
{
    //cout << "     HBSplineFEM: generating coefficient Matrices ...\n\n";

    char fct[] = "HBSplineFEM::calcStiffnessAndResidual";
    computerTime.go(fct);

    int bb, ee, ii, jj, kk, ll, nr, nc, aa, dd, numPoints, size2, r, c, gp;

    time_t tstart, tend;
    
    IterNum   = (iterCount == 1);

    if(firstIter)
    {
      rNorm = -1.0;
    }

    solverEigen->zeroMtx();

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the background fluid grid
    //
    ////////////////////////////////////////////////////////

    tstart = time(0);
    for(ii=0;ii<activeElements.size();ii++)
    {
      node *nd = elems[activeElements[ii]];
      //cout << " nd->getID() " <<  nd->getID() << '\t' <<  nd->getLevel() << endl;

      //printVector(nd->forAssyVec);
      //printVector(nd->forAssyVec2);

      int nr = nd->forAssyVec.size();

      MatrixXd  Klocal;
      VectorXd  Flocal;

      Klocal = MatrixXd::Zero(nr, nr);
      Flocal = VectorXd::Zero(nr);

      nd->calcStiffnessAndResidualGFEM(Klocal, Flocal);
      nd->applyDirichletBCsGFEM(Klocal, Flocal);
      //nd->applyNeumannBCsGFEM(Klocal, Flocal);

      solverEigen->assembleMatrixAndVector(0, 0, nd->forAssyVec, Klocal, Flocal);
    }

    //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec(0)), totalDOF);

    //applyBoundaryConditions();
    //addExternalForces();

    tend = time(0); 
    //cout << "HBSplineFEM::calcStiffnessAndResidual()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed boundary points
    // Lagrange multipliers
    //  or
    // Penalty method
    //
    ////////////////////////////////////////////////////////

    if(DIM == 2)
      applyInterfaceTerms2D();

    if(DIM == 3)
      applyInterfaceTerms3D();

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed solid body
    // only for the MONOLITHIC_SCHEME
    //
    ////////////////////////////////////////////////////////

    if(!STAGGERED)
    {
      //cout << " Terms related to monolithic scheme " << endl;

      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
          kk = fluidDOF + IBDOF;
          ImmersedBodyObjects[bb]->assembleGlobalMatrixAndVector(fluidDOF, kk, solverEigen->mtx, &(solverEigen->rhsVec(0)));
      }

      kk = fluidDOF + IBDOF;
      for(bb=0;bb<contactElementObjects.size();bb++)
      {
        contactElementObjects[bb]->calcStiffnessAndResidual();
        contactElementObjects[bb]->assembleMatrixAndVector(kk, solverEigen->mtx, &(solverEigen->rhsVec(0)));
      }
    }

    //printf("\n rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    //cout << " rhsVec " << endl;
    //for(int ii=fluidDOF-10;ii<totalDOF;ii++)
      //printf("%5d \t %12.8f \n", ii, rhsVec(ii));

    rNormPrev = rNorm;
    rNorm     = solverEigen->rhsVec.norm();

    solverEigen->currentStatus = ASSEMBLY_OK;

    ctimCalcStiffRes += computerTime.stop(fct);
    //computerTime.stopAndPrint(fct);

    //cout << "     HBSplineFEM: generating coefficient Matrices ...DONE  \n\n";

    return 0;
}





void  HBSplineFEM::applyBoundaryConditions()
{
    //cout << "     HBSplineFEM: applyBoundaryConditions ... STARTED  \n\n";
    //firstIter = true;

    //cout << "     HBSplineFEM: applyBoundaryConditions ...DONE  \n\n";
    return;
}



void HBSplineFEM::addExternalForces()
{
//  cout << "HBSplineFEM::addExternalForces() " << endl;
  // if(firstIter)     calcAndAssyLoadVector(1.0, 0.0);

  return;
}



int HBSplineFEM::factoriseSolveAndUpdate()
{
   cout << "     HBSplineFEM::factoriseSolveAndUpdate() ...  \n\n";
   //char fct[] = "HBSplineFEM::factoriseSolveAndUpdate";
   //computerTime.go(fct);

   time_t tstart, tend;

   int ii, kk, ee;

  //cout << " rhsVec " << endl;        printVector(&(solverEigen->rhsVec(0)), totalDOF);
  //cout << " rhsVec " << endl;
  //for(int ii=fluidDOF-20;ii<totalDOF;ii++)
    //printf("%5d \t %12.8f \n", ii, solverEigen->rhsVec(ii));

  tstart = time(0);

  solverEigen->factoriseAndSolve();

  tend = time(0);
  printf("HBSplineFEM::factoriseSolveAndUpdate() took %8.4f second(s) \n ", difftime(tend, tstart) );

  //cout << " result " << endl;
  //for(int ii=fluidDOF;ii<totalDOF;ii++)
    //printf("%5d \t %12.8f \n", ii, solverEigen->soln(ii));

   
   double  *sln = &(solverEigen->soln[0]);

   for(ii=0;ii<velDOF;ii++)
     SolnData.var1[ii] += sln[ii];

   kk = velDOF;
   for(ii=0;ii<presDOF;ii++)
     SolnData.var2[ii] += sln[kk+ii];

   kk = velDOF + presDOF;
   for(ii=0;ii<IBDOF;ii++)
     SolnData.var3[ii] += sln[kk+ii];

   if(!STAGGERED)
   {
     kk = velDOF + presDOF + IBDOF;
     for(ii=0;ii<solidDOF;ii++)
       SolnData.var4[ii] += sln[kk+ii];
   }

   //computerTime.stopAndPrint(fct);

   //ctimFactSolvUpdt += computerTime.stop(fct);

   //cout << " Lagrange multipliers " << endl;        printVector(SolnData.var3);

   //for(ii=0;ii<IBDOF;ii++)
     //printf("%5d \t %20.16f \n", ii, SolnData.var3(kk));

   //cout << "     HBSplineFEM: solving the matrix system ...DONE  \n\n";

   //computerTime.stopAndPrint(fct);

   return 0;
}



void HBSplineFEM::computeElementErrors(int index)
{
    int  ii, ee, count=0;
    node  *nd;

    totalError = 0.0;
    
    if(index < 4) // L2 or H1 norm based errors
    {
       for(ii=0;ii<activeElements.size();ii++)
       {
          ee = activeElements[ii];

          //cout << ee << '\t' << elems[ee]->getError() << endl;
          //totalError += ( elems[ee]->getError() * elems[ee]->getError() );
          totalError +=  elems[ee]->calcError(index);
       }
       totalError = sqrt(totalError);
    }
    else // gradient based error
    {
       for(ii=0;ii<activeElements.size();ii++)
       {
          nd = elems[activeElements[ii]];

          totalError += nd->calcError(index);

          count++;
       }
       totalError /= count;
    }
    
    //printf(" \n\n \t totalError = %12.6E \n\n " , totalError);
    
    if(index < 3)
      printf(" \n\n \t L2 Error = %12.6E \n\n " , totalError);
    else
      printf(" \n\n \t H1 Error = %12.6E \n\n " , totalError);

    char        tmp[50];
    MyString    tmpStr, fname;

    sprintf(tmp,"\n %5d \t %5d \t %12.6E ", CURRENT_LEVEL, totalDOF, totalError);
    tmpStr.append(tmp);

    prgWriteToTFile(tmpStr);

    return;
}



void HBSplineFEM::computeConditionNumber()
{
  solverEigen->computeConditionNumber();

  return;
}





