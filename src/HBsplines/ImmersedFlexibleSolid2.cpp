
#include "ImmersedFlexibleSolid.h"
#include "SolverMA41Eigen.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "ImmersedIntegrationElement.h"
#include "SolverPetsc.h"
#include "QuadratureUtil.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"


extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



void ImmersedFlexibleSolid::initialise()
{
  SolnData.STAGGERED = STAGGERED;

  setSolver(1);
  setTimeParam();
  computeInitialAcceleration();

  return;
}




void ImmersedFlexibleSolid::setSolver(int slv, int *parm, bool cIO)
{
    //if(solver != NULL)
      //delete solver;
    //solver = NULL;

    int numProc=0;

    switch(slv)
    {
        case  1: // MA41 ..........................

            solver = (SolverEigen*) new SolverMA41Eigen;

            prepareMatrixPattern();

            //solver->printInfo();

            if(solver->initialise(0,0,totalDOF) != 0)
              return;

        break;

        case  4: // SolverEigen ..........................

            solver = (SolverEigen*) new SolverEigen;

            solver->setAlgorithmType(1);

            prepareMatrixPattern();

            if(solver->initialise(0,0,totalDOF) != 0)
              return;

            solver->printInfo();

        break;

        default: // invalid slv ...................

             cout << " this solver has not been implemented yet!\n\n";

        break;
    }

    solverOK = true;

    if(solver != NULL)
      solver->checkIO = cIO;

    //if(tis > 0)
      //setInitialConditions();

    for(int ii=0;ii<nElem;ii++)
    {
      elems[ii]->prepareElemData();
    }

    return;
}



void ImmersedFlexibleSolid::prepareMatrixPattern()
{
    printf("\n     ImmersedFlexibleSolid::prepareMatrixPattern()  .... STARTED ...\n");

    int  r=0, c=0, r1=0, c1=0, count=0, count1=0, count2=0;
    int  ii=0, jj=0, iii=0, e=0, ee=0, ind=0, val1=0, val2=0;
    int  *tt, *tt1, *tt2, n1=0, n2=0, kk=0, e1=0, a=0, b=0, ll=0, pp=0, nnz=0, nn;
    int  nRow=0, nCol=0, ind1=0, ind2=0, ndof_temp1=0, npElem, dof, dim, row, col;

    ndof = elems[0]->getNdofPerNode();
    ndof_temp1 = ndof;

    //totalDOF = nNode*ndof;     totalDOF += nElem_Constraint;

    node_map_new_to_old.resize(nNode, 0);
    node_map_old_to_new.resize(nNode, 0);

    dof_map_new_to_old.resize(totalDOF, 0);
    dof_map_old_to_new.resize(totalDOF, 0);

    kk=0;
    for(ii=0; ii<nNode; ii++)
    {
      node_map_new_to_old[ii] = ii;
      node_map_old_to_new[ii] = ii;

      for(jj=0; jj<ndof; jj++)
      {
        dof_map_new_to_old[kk] = kk;
        dof_map_old_to_new[kk] = kk;
        kk++;
      }
    }

    SolnData.node_map_new_to_old = node_map_new_to_old;
    SolnData.node_map_old_to_new = node_map_old_to_new;

    GeomData.node_map_new_to_old = node_map_new_to_old;
    GeomData.node_map_old_to_new = node_map_old_to_new;

    //cout << " nelem and  ndof " << nElem << '\t' << ndof << '\t' << npElem << endl;
    cout << " ImmersedFlexibleSolid.... nNode and  ndof " << nNode << '\t' << ndof << endl;

    totalDOF = 0;
    for(ii=0; ii<nNode; ii++)
    {
      for(jj=0; jj<ndof; jj++)
      {
        //cout << ii << '\t' << jj << '\t' << NodeType[ii][jj] << endl;
        if(NodeType[ii][jj] == (int) -7777)
        {
          ID[ii][jj] = totalDOF++;

          assy4r.push_back(ii*ndof + jj);
        }
      }
    }
    cout << " totalDOF = " << totalDOF << endl;

    // adjust assy4r array to account for the constraints

    ind = nNode*ndof_temp1;
    for(ii=0; ii<nElem_Constraint; ii++)
    {
      assy4r.push_back(ind+ii);
    }

    GeomData.assy4r = assy4r;

    // adjust ID array to account for the constraints

    for(ee=0; ee<nElem; ee++)
    {
      ii = elems[ee]->nodeNums[0];

      pp = elems[ee]->getElmTypeNameNum();

      //cout << " ee = " << ee << '\t' << pp << endl;

      if( (pp == 33) || (pp == 35) ) // contact element along X axis
      {
        jj = ndof_temp1 + 0;
        ID[ii][jj] = totalDOF++;
      }

      if( (pp == 34) || (pp == 36) ) // contact element along Y axis
      {
        jj = ndof_temp1 + 1;
        ID[ii][jj] = totalDOF++;
      }

      if( pp == 37 )  // contact element along Z axis
      {
        jj = ndof_temp1 + 2;
        ID[ii][jj] = totalDOF++;
      }
    }

    cout << " totalDOF " << totalDOF << endl;


    for(ee=0; ee<nElem; ee++)
    {
      npElem     = elems[ee]->getNodesPerElement();
      ndof_temp1 = elems[ee]->getNdofPerNode();
      ind        = ndof_temp1*npElem;

      //printVector(IEN[ee]);

      pp = elems[ee]->getElmTypeNameNum();

      if( pp == 33 ) // contact element along X-axis, 2D
      {
        LM[ee].resize(2);  // to account for the Lagrange multiplier

        kk = IEN[ee][0];

        LM[ee][0] = ID[kk][0];
        LM[ee][1] = ID[kk][2];
      }
      else if( pp == 34 ) // contact element along Y-axis, 2D
      {
        LM[ee].resize(2);  // to account for the Lagrange multiplier

        kk = IEN[ee][0];

        LM[ee][0] = ID[kk][1];
        LM[ee][1] = ID[kk][3];
      }
      else if( pp == 35 ) // contact element along X-axis, 3D
      {
        LM[ee].resize(2);  // to account for the Lagrange multiplier

        kk = IEN[ee][0];

        LM[ee][0] = ID[kk][0];
        LM[ee][1] = ID[kk][3];
      }
      else if( pp == 36 ) // contact element along Y-axis, 3D
      {
        LM[ee].resize(2);  // to account for the Lagrange multiplier

        kk = IEN[ee][0];

        LM[ee][0] = ID[kk][1];
        LM[ee][1] = ID[kk][4];
      }
      else if( pp == 37 ) // contact element along Z-axis, 3D
      {
        LM[ee].resize(2);  // to account for the Lagrange multiplier

        kk = IEN[ee][0];

        LM[ee][0] = ID[kk][2];
        LM[ee][1] = ID[kk][5];
      }
      else // standard elements
      {
        LM[ee].resize(ind);
        for(ii=0;ii<npElem;ii++)
        {
          ind = ndof_temp1*ii;

          kk = IEN[ee][ii];

          for(jj=0;jj<ndof_temp1;jj++)
          {
            LM[ee][ind+jj] = ID[kk][jj];
          }
        }
      }
    }

    //cout << " totalDOF " << totalDOF << endl;
    //cout << " nElem " << nElem << endl;
    for(ee=0; ee<nElem; ee++)
    {
      elems[ee]->forAssyVec = LM[ee];
    }

    //pp=false;
    pp=true;
    if(pp)
    {
       printf("   IEN array \n\n");
       for(ee=0; ee<nElem; ee++)
       {
          npElem = elems[ee]->getNodesPerElement();
          for(ii=0; ii<npElem; ii++)
            cout << '\t' << IEN[ee][ii];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   ID array \n\n");
       for(ii=0; ii<nNode; ii++)
       {
          for(jj=0; jj<ID[ii].size(); jj++)
            cout << '\t' << ID[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   LM array \n\n");
       for(ee=0; ee<nElem; ee++)
       {
          for(jj=0; jj<LM[ee].size(); jj++)
            cout << '\t' << LM[ee][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("  assy4r array \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << ii << '\t' << assy4r[ii] << endl;
       }
       printf("\n\n\n");
    }

    //printf("\n element DOF values initialised \n\n");
    //printf("\n Preparing matrix pattern \n\n");

    vector<int>::const_iterator location;
    set<int>::iterator it;

    forAssyMat.resize(totalDOF);

    for(ee=0;ee<nElem;ee++)
    {
       tt = &(LM[ee][0]);
       nsize = LM[ee].size();

       for(ii=0;ii<nsize;ii++)
       {
          count1 = tt[ii];

          if(tt[ii] != -1)
          {
            for(jj=0;jj<nsize;jj++)
            {
              if(tt[jj] != -1)
                forAssyMat[count1].push_back(tt[jj]);
            }
          }
       }
    }

    bool pp1=false;
    //pp1=true;
    if(pp1)
    {
       printf("   Number of non-zeros = %5d \n\n", nnz);
       printf("   dof to dof connectivity ...:  \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << " dof # " << ii << " : ";
          for(jj=0;jj<forAssyMat[ii].size();jj++)
            cout << '\t' << forAssyMat[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");
    }

    //printf("\n Preparing matrix pattern DONE \n\n");

    VectorXd  nnzVec(totalDOF);

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      forAssyMat[ii].push_back(ii);
      findUnique(forAssyMat[ii]);
      nnzVec[ii] = forAssyMat[ii].size();
      nnz += nnzVec[ii];
    }
    //cout << " nnz " << nnz << endl;

    nRow = nCol = totalDOF;

    //cout << " AAAAAAAAAA " << endl;
    solver->mtx.setZero();

    solver->mtx.resize(nRow, nCol);
    solver->mtx.reserve(nnz);
    solver->mtx.reserve(nnzVec);

    for(ii=0;ii<totalDOF;ii++)
    {
      for(jj=0;jj<forAssyMat[ii].size();jj++)
      {
        //cout << ii << '\t' << forAssyMat[ii][jj] << endl;
        solver->mtx.coeffRef(ii, forAssyMat[ii][jj]) = 0.0;
      }
    }

    solver->mtx.makeCompressed();

    solver->currentStatus = PATTERN_OK;

    cout << " STAGGERED = " << STAGGERED << endl;


    if(!STAGGERED)
    {
      forAssyCoupledHorz.resize(totalDOF);
      forAssyCoupledVert.resize(nNode*DIM);
      
      vector<int>  nodeNums, pointNums;
      vector<int>  forAssyVec, posIndices;

      for(ee=0; ee<nElem; ee++)
      {
        nodeNums  = elems[ee]->nodeNums;
        //pointNums = ImmIntgElems[ee]->pointNums;

        forAssyVec  = elems[ee]->forAssyVec;
        //posIndices = ImmIntgElems[ee]->posIndices;

        printVector(nodeNums);
        printVector(forAssyVec);
        //printVector(posIndices);
        
        int size1 = forAssyVec.size();
        int nlbS  = nodeNums.size();

        for(ii=0; ii<size1; ii++)
        {
          row = forAssyVec[ii];

          if(row > -1)
          {
              for(jj=0; jj<nlbS; jj++)
              {
                ind2 = nodeNums[jj]*DIM;
                //cout << ii << '\t' << row << '\t' << jj << '\t' << ind2 << endl;

                forAssyCoupledHorz[row].push_back(ind2);
                forAssyCoupledHorz[row].push_back(ind2+1);

                forAssyCoupledVert[ind2].push_back(row);
                forAssyCoupledVert[ind2+1].push_back(row);
              }
          }
        }
      }
      //cout << " forAssyCoupledHorz ...... " << endl;      cout << endl;
      for(ii=0;ii<totalDOF;ii++)
      {
        findUnique(forAssyCoupledHorz[ii]);
        printVector(forAssyCoupledHorz[ii]);
      }
    }

/*
    if(!STAGGERED)
    {
      forAssyCoupledHorz.resize(totalDOF);
      forAssyCoupledVert.resize(nNode*DIM);
      
      vector<int>  nodeNums, pointNums;
      vector<int>  forAssyVec, posIndices;

      for(ee=0; ee<nElem; ee++)
      {
        nodeNums  = elems[ee]->nodeNums;
        pointNums = ImmIntgElems[ee]->pointNums;

        forAssyVec  = elems[ee]->forAssyVec;
        posIndices = ImmIntgElems[ee]->posIndices;
        
        printVector(forAssyVec);
        printVector(posIndices);
        
        int nlbS = forAssyVec.size();
        int nlbL = posIndices.size();

        for(ii=0; ii<nlbS; ii++)
        {
          row = forAssyVec[ii];

          if(row > -1)
          {
              for(jj=0; jj<nlbL; jj++)
              {
                //cout << ii << '\t' << row << '\t' << jj << '\t' << ind2 << endl;

                forAssyCoupledHorz[row].push_back(posIndices[jj]);

                forAssyCoupledVert[posIndices[jj]].push_back(row);
              }
          }
        }
      }
      //cout << " forAssyCoupledHorz ...... " << endl;      cout << endl;
      for(ii=0;ii<totalDOF;ii++)
      {
        findUnique(forAssyCoupledHorz[ii]);
        printVector(forAssyCoupledHorz[ii]);
      }
    }
*/

    printf("\n     ImmersedFlexibleSolid::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return;
}




int ImmersedFlexibleSolid::solveTimeStep()
{
    printf("\n Solving Immersed Flexible Solid \n");
    
    bool CONVERGED = false;

    for(int iter=1; iter<=10;iter++)
    {
      updateIterStep();

      calcStiffnessAndResidual(1, 1, 1);

      printf("\t ImmersedFlexibleSolid:   %5d \t %12.6E \n", iter, rNorm);

      if(converged())
      {
        CONVERGED = true;
        break;
      }

      factoriseSolveAndUpdate();
    }

    printf("\n Solving Immersed Flexible Solid ..... DONE  \n\n");

    if(CONVERGED)
      return 0;
    else
      return 1;
}






int ImmersedFlexibleSolid::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
    cout << "     ImmersedFlexibleSolid::calcStiffnessAndResidual ..." << endl;

    if(solver == NULL)
    {
      COUT << "You need to select a solver first!\n\n";
      return 1;
    }

    solver->zeroMtx();

    if(firstIter)
      rNorm = -1.0;

    MatrixXd  Klocal;
    VectorXd  Flocal;

    for(int ee=0;ee<nElem;ee++)  // loop over all the elements
    {
      //cout << "       elem... : " << (ee+1) << endl;

      elems[ee]->calcStiffnessAndResidual(Klocal, Flocal, firstIter);

      solver->assembleMatrixAndVector(0, 0, elems[ee]->forAssyVec, Klocal, Flocal);
    }

    //cout << " solver->rhsVec " << endl;        printVector(solver->rhsVec);

    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    //applyBoundaryConditions(0, 0, solver->mtx, &(solver->rhsVec(0)));

    applyExternalForces();

    if(mpapTime.cur >= 100000.0)  addControlTerms();


    //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), totalDOF);

    firstIter = false;
    SolnData.firstIter = firstIter;
    rNormPrev = rNorm;
    rNorm     = solver->rhsVec.norm();

    solver->currentStatus = ASSEMBLY_OK;

    return 0;
}




void  ImmersedFlexibleSolid::addControlTerms()
{
    int  nboundnodes = 11, nn, n1, n2, ii;
    int boundnodes[] = {201, 402, 603, 804, 1005, 1206, 1407, 1608, 1809, 2010, 2211};
    //double  Kd = 0.0/2.0, Kv = 0.5/2.0;
    double  Kd = SolnData.MatlProp[0].data[4]/2.0, Kv = SolnData.MatlProp[0].data[5]/2.0;
    //cout <<  Kd << '\t' << Kv << endl;


    double  stiffnessFact = SolnData.td[2];
    double  dampingFact   = SolnData.td[6];

    for(ii=0; ii<nboundnodes; ii++)
    {
      nn = boundnodes[ii] - 1;

      n1 = nn*2;
      n2 = ID[nn][0];

      if(n2 != -1)
      {
        solver->mtx.coeffRef(n2,  n2)   += (stiffnessFact*Kd);
        solver->rhsVec[n2]              -= (Kd*SolnData.var1Cur[n1]);

        solver->mtx.coeffRef(n2,  n2)   += (dampingFact*Kv);
        solver->rhsVec[n2]              -= (Kv*SolnData.var1DotCur[n1]);
      }

      n1 = nn*2+1;
      n2 = ID[nn][1];

      if(n2 != -1)
      {
        solver->mtx.coeffRef(n2,  n2)   += (stiffnessFact*Kd);
        solver->rhsVec[n2]              -= (Kd*SolnData.var1Cur[n1]);

        solver->mtx.coeffRef(n2,  n2)   += (dampingFact*Kv);
        solver->rhsVec[n2]              -= (Kv*SolnData.var1DotCur[n1]);
      }
    }

    return;
}



int ImmersedFlexibleSolid::applyBoundaryConditions(int start1, int start2, SparseMatrixXd& globalK, double* rhs)
{
    int  ii=0, jj=0, nn=0, dof=0, aa=0, ind=0;
    double  specVal=0.0, PENALTY=1.0e4;
    double  af = SolnData.td(2);

    for(aa=0;aa<DirichletBCs.size();aa++)
    {
      nn  = (int) (DirichletBCs[aa][0] - 1);
      dof = (int) (DirichletBCs[aa][1] - 1);
      specVal = DirichletBCs[aa][2];

      ind = nn*ndof+dof;
      specVal  -=  SolnData.var1DotCur[ind];
      //specVal  -=  SolnData.var1Cur[ind];

      //cout << start1 << '\t' << nn << '\t' << ind << '\t' << specVal << endl;

      ind += start2;
      globalK.coeffRef(ind, ind) += af*PENALTY;
      rhs[ind]   += (PENALTY*specVal);
    }

    return 0;
}



int  ImmersedFlexibleSolid::applyExternalForces()
{
    //printVector(SolnData.forceCur);

    for(int ii=0;ii<totalDOF;ii++)
    {
      solver->rhsVec[ii] += SolnData.forceCur[assy4r[ii]];
      //solver->rhsVec[ii] += fluidAcceCur[assy4r[ii]];
    }

    return 0;
}



void ImmersedFlexibleSolid::calcForceVector(double* rhs)
{
    return;
}




void ImmersedFlexibleSolid::calcForceVectorMonolithic(double* rhs)
{
  int  ee, gp, ii, jj, TI, TIp1, TIp2, TJ, TJp1;

  int  nlbS = ImmIntgElems[0]->pointNums.size();
  int  nlbL = elems[0]->nodeNums.size();

  vector<double>  NL(nlbL), dN(nlbL), dN_dx(nlbL), xNode(nlbL), yNode(nlbL), NS(nlbS);
  double  dvol, detJ, dvol1, dvol2, fact1, fact2, xc;

  vector<int>  nodeNums, pointNums;

  VectorXd  rhsTemp(nNode*ndof);
  rhsTemp.setZero();

  bool axsy = false;//(GeomDataHBS->FluidProps[2] == 1);

  int nGP=2;
  vector<double>  gausspoints(nGP), gaussweights(nGP);
  getGaussPoints1D(nGP, gausspoints, gaussweights);


  for(ee=0; ee<nElem; ee++)
  {
      nodeNums = elems[ee]->nodeNums;
      pointNums = ImmIntgElems[ee]->pointNums;
      
      for(ii=0;ii<nlbL;ii++)
      {
        xNode[ii] = GeomData.NodePosCur[nodeNums[ii]][0];
        yNode[ii] = GeomData.NodePosCur[nodeNums[ii]][1];
      }

      for(gp=0;gp<gausspoints.size();gp++)
      {
        computeLagrangeBFsLine2D(nlbL-1, gausspoints[gp], &xNode[0], &yNode[0], &NL[0], &dN[0], detJ);

        //cout << " detJ " << detJ << '\t' << gaussweights[gp] << endl;

        dvol  = gaussweights[gp] * detJ;

        xc = 0.0;
        for(ii=0;ii<nlbL;ii++)
          xc += NL[ii] * xNode[ii];

        if(axsy)
          dvol1 = dvol * 2.0*PI*xc;


        if(nlbS == nlbL)
        {
          for(ii=0; ii<nlbL; ii++)
            NS[ii] = NL[ii];
        }
        else
        {
          NS[0] = 1.0;
        }

        for(ii=0; ii<nlbS; ii++)
        {
          fact1 = NS[ii]*dvol1;
          fact2 = NS[ii]*dvol2;

          TI   = ii*3;
          TIp1 = TI+1;
          TIp2 = TI+2;
        }

        for(ii=0; ii<nlbL; ii++)
        {
          fact1 = NS[ii]*dvol1;
          fact2 = NS[ii]*dvol2;

          TI   = ii*3;
          TIp1 = TI+1;
          TIp2 = TI+2;

        }
    }
  }

  return;
}




int ImmersedFlexibleSolid::factoriseSolveAndUpdate()
{
    cout << " ImmersedFlexibleSolid::factoriseSolveAndUpdate ... " << endl;
    //cout << " residue_new " << endl;        printVector(&(solver->rhsVec[0]), totalDOF);

    solver->factoriseAndSolve();

    soln.setZero();
    // update solution vector
    for(int kk=0;kk<totalDOF;kk++)
      soln[assy4r[kk]] = solver->soln[kk];

    SolnData.var1 += soln;

    //printVector(SolnData.var1);

    //printf("\n\n\n");

    //cout << " result " << endl;        printVector(&(SolnData.var1[0]), totalDOF);

    return 0;
}


/*
void ImmersedFlexibleSolid::calcCouplingMatrices()
{
    //////////////////////////////////////////
    // off-diagonal matrices
    //////////////////////////////////////////

    int  aa=0, nlb=0, ii=0, jj=0, kk=0, nlbS=0, nlbL=0;
    int  TI=0, TIp1=0, TIp2=0, TJ=0, TJp1=0;
    MatrixXd  Ktemp, Ktemp2;

    int  ind1 = nNode*ndof;
    int  ind2 = nNode*DIM;

    Khorz.resize(ind1, ind2);
    Khorz.setZero();

    Kvert.resize(ind2, ind1);
    Kvert.setZero();

    for(aa=0;aa<ImmIntgElems.size();aa++)
    {
      nlbL = ImmIntgElems[aa]->pointNums.size();
      nlbS = elems[aa]->nodeNums.size();

      //cout << " aa = " << aa << endl;
      ImmIntgElems[aa]->computeKhorzKvertFlexible(0, 0, Ktemp, Ktemp2);
      //printMatrix(Ktemp);
      //printf("\n\n");
      //printMatrix(Ktemp2);
      //printf("\n\n");

      for(ii=0; ii<nlbS; ii++)
      {
        TI   = 3*elems[aa]->nodeNums[ii];
        TIp1 = TI+1;
        TIp2 = TI+2;
        
        ind1 = 3*ii;

        for(jj=0; jj<nlbL; jj++)
        {
          TJ   = ImmIntgElems[aa]->pointNums[jj] * DIM;
          TJp1 = TJ+1;

          ind2 = jj*DIM;

          Khorz(TI,   TJ)    +=  Ktemp(ind1,   ind2);
          Khorz(TI,   TJp1)  +=  Ktemp(ind1,   ind2+1);

          Khorz(TIp1, TJ)    +=  Ktemp(ind1+1, ind2);
          Khorz(TIp1, TJp1)  +=  Ktemp(ind1+1, ind2+1);

          Khorz(TIp2, TJ)    +=  Ktemp(ind1+2, ind2);
          Khorz(TIp2, TJp1)  +=  Ktemp(ind1+2, ind2+1);

          Kvert(TJ,   TI)    +=  Ktemp2(ind2,   ind1);
          Kvert(TJp1, TI)    +=  Ktemp2(ind2+1, ind1);

          Kvert(TJ,   TIp1)  +=  Ktemp2(ind2,   ind1+1);
          Kvert(TJp1, TIp1)  +=  Ktemp2(ind2+1, ind1+1);

          Kvert(TJ,   TIp2)  +=  Ktemp2(ind2,   ind1+2);
          Kvert(TJp1, TIp2)  +=  Ktemp2(ind2+1, ind1+2);
        }
      }
    }

    Khorz *= -SolnData.td[2];
    Kvert *= -SolnData.td[2];

    //printMatrix(Khorz);
    //printMatrix(Kvert);
}
*/





void ImmersedFlexibleSolid::calcCouplingMatrices()
{
    //////////////////////////////////////////
    // off-diagonal matrices
    //////////////////////////////////////////

    int  ee=0, ii=0, jj=0, kk=0, nlbS=0, nlbL=0;
    int  TI=0, TIp1=0, TIp2=0, TJ=0, TJp1=0;
    MatrixXd  Ktemp(6,4), Ktemp2(4,6);
    Ktemp.setZero();
    Ktemp2.setZero();

    int  ind1 = nNode*ndof;
    int  ind2 = nNode*DIM;

    Khorz.resize(ind1, ind2);
    Khorz.setZero();

    Kvert.resize(ind2, ind1);
    Kvert.setZero();
    
    vector<int>  nodeNums;
    
    double  xNode[2], yNode[2], he, fact1, fact2;

    nlbL = ImmIntgElems[0]->pointNums.size();

    for(ee=0; ee<nElem; ee++)
    {
      nodeNums = elems[ee]->nodeNums;

      nlbS = nodeNums.size();

      xNode[0] = GeomData.NodePosCur[nodeNums[0]][0];
      yNode[0] = GeomData.NodePosCur[nodeNums[0]][1];
      xNode[1] = GeomData.NodePosCur[nodeNums[1]][0];
      yNode[1] = GeomData.NodePosCur[nodeNums[1]][1];

      he = sqrt( (xNode[1]-xNode[0])*(xNode[1]-xNode[0]) + (yNode[1]-yNode[0])*(yNode[1]-yNode[0]));

      ImmIntgElems[ee]->computeKhorzKvertFlexible(0, 0, Ktemp, Ktemp2);

      //printMatrix(Ktemp);
      //printf("\n\n");
      
      fact1 = he/3.0;
      fact2 = he/6.0;

      //Ktemp(0,0) = fact1;      Ktemp(0,2) = fact2;
      //Ktemp(1,1) = fact1;      Ktemp(1,3) = fact2;

      //Ktemp(3,0) = fact2;      Ktemp(3,2) = fact1;
      //Ktemp(4,1) = fact2;      Ktemp(4,3) = fact1;

      //printMatrix(Ktemp);
      //printf("\n\n");
      //printMatrix(Ktemp2);
      //printf("\n\n");

      for(ii=0; ii<nlbS; ii++)
      {
        TI   = 3*nodeNums[ii];
        TIp1 = TI+1;
        TIp2 = TI+2;
        
        ind1 = 3*ii;

        for(jj=0; jj<nlbL; jj++)
        {
          //TJ   = ImmIntgElems[ee]->pointNums[jj] * DIM;
          TJ   = nodeNums[jj] * DIM;
          TJp1 = TJ+1;

          ind2 = jj*DIM;

          Khorz(TI,   TJ)    +=  Ktemp(ind1,   ind2);
          Khorz(TI,   TJp1)  +=  Ktemp(ind1,   ind2+1);

          Khorz(TIp1, TJ)    +=  Ktemp(ind1+1, ind2);
          Khorz(TIp1, TJp1)  +=  Ktemp(ind1+1, ind2+1);

          Khorz(TIp2, TJ)    +=  Ktemp(ind1+2, ind2);
          Khorz(TIp2, TJp1)  +=  Ktemp(ind1+2, ind2+1);

          //Kvert(TJ,   TI)    +=  Ktemp2(ind2,   ind1);
          //Kvert(TJp1, TI)    +=  Ktemp2(ind2+1, ind1);

          //Kvert(TJ,   TIp1)  +=  Ktemp2(ind2,   ind1+1);
          //Kvert(TJp1, TIp1)  +=  Ktemp2(ind2+1, ind1+1);

          //Kvert(TJ,   TIp2)  +=  Ktemp2(ind2,   ind1+2);
          //Kvert(TJp1, TIp2)  +=  Ktemp2(ind2+1, ind1+2);
        }
      }
    }

    Khorz *= -SolnData.td[2];

    Kvert = Khorz.transpose();

    //printMatrix(Khorz);
    //printMatrix(Kvert);
}





int ImmersedFlexibleSolid::assembleGlobalMatrixAndVector(int start1, int start2, SparseMatrixXd& mtx, double* rhs)
{
    if(totalDOF <= 0)
      return 0;

    int ee=0, ii=0, jj=0, k1=0, k2=0, r=0, c=0, aa, bb;

    vector<int> forAssyVec;
    MatrixXd  Klocal;
    VectorXd  Flocal;

    for(ee=0;ee<nElem;ee++)  // loop over all the elements
    {
      //cout << "       elem... : " << (ee+1) << endl;

      elems[ee]->calcStiffnessAndResidual(Klocal, Flocal, false);
      //printMatrix(Klocal);
      //printVector(Flocal);

      //cout << " MMMMMMMMMMM " << endl;
      //elems[ee]->assembleElementMatrixAndVector(start2, mtx, rhs);
      
      forAssyVec = elems[ee]->forAssyVec;
      //printVector(forAssyVec);

      for(ii=0;ii<nsize;ii++)
      {
        aa = forAssyVec[ii];
        if( aa != -1 )
        {
          r = start2 + aa;
          rhs[r] += Flocal(ii);
          for(jj=0;jj<nsize;jj++)
          {
            bb = forAssyVec[jj];
            if( bb != -1 )
              mtx.coeffRef(r, start2+bb) += Klocal(ii,jj);
          }
        }
      }
    }

    //cout << " aaaaaaaaaaaaaaa " << endl;
    //printVector(SolnData.forceCur);

    updateForce();

    SolnData.force = SolnData.forceTemp;
    SolnData.forceCur = SolnData.td[2]*SolnData.force + (1.0-SolnData.td[2])*SolnData.forcePrev;

    //printVector(SolnData.forceTemp);

    for(ii=0;ii<totalDOF;ii++)
    {
      //rhs[start2+ii] += SolnData.forceTemp[assy4r[ii]];
      rhs[start2+ii] += SolnData.forceCur[assy4r[ii]];
    }

    //cout << " applyBoundaryConditions " << endl;
    //applyBoundaryConditions(start1, start2, mtx, rhs);


    //cout << " calcCouplingMatrices " << endl;
    calcCouplingMatrices();
    //cout << " calcCouplingMatrices " << endl;
    //printMatrix(Khorz);


    for(ii=0;ii<totalDOF;ii++)
    {
      r  = start2 + ii;
      k1 = assy4r[ii];

      for(jj=0;jj<forAssyCoupledHorz[ii].size();jj++)
      {
        k2 = forAssyCoupledHorz[ii][jj];
        c = start1 + k2;

        mtx.coeffRef(r, c) += Khorz(k1, k2);
        mtx.coeffRef(c, r) += Khorz(k1, k2);
      }
    }

    return 0;
}




int ImmersedFlexibleSolid::assembleGlobalMatrixAndVectorCutFEM(int start1, int start2, SolverPetsc* solverTemp)
{
    int ee=0, ii=0, jj=0, ind=0;

    MatrixXd  Klocal;
    VectorXd  Flocal;

    for(ee=0; ee<nElem; ee++)  // loop over all the elements
    {
      //cout << "       elem... : " << (ee+1) << endl;

      elems[ee]->calcStiffnessAndResidual(Klocal, Flocal);

      solverTemp->assembleMatrixAndVector(start1, start2, elems[ee]->forAssyVec, elems[ee]->forAssyVec, Klocal, Flocal);
    }

    //cout << " solver->rhsVec " << endl;        printVector(solver->rhsVec);

    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    applyBoundaryConditions(start1, start2, solverTemp->mtx, solverTemp->rhsVec);

    //applyExternalForces();

    //cout << " aaaaaaaaaaaaaaa " << endl;
    for(ii=0;ii<totalDOF;ii++)
    {
      ind = start2+ii;
      VecSetValue(solverTemp->rhsVec, ind, SolnData.forceCur[assy4r[ii]], ADD_VALUES);
    }

    return 0;
}




int ImmersedFlexibleSolid::applyBoundaryConditions(int start1, int start2, Mat mtxTemp, Vec rhsTemp)
{
    int  ii=0, jj=0, nn=0, dof=0, aa=0, ind=0;
    double  specVal=0.0, PENALTY=1.0e8, stiff=0.0, res=0.0;

    double  af = SolnData.td(2);
    double  velFact = SolnData.td(10);

    for(aa=0;aa<DirichletBCs.size();aa++)
    {
      nn  = (int) (DirichletBCs[aa][0] - 1);
      dof = (int) (DirichletBCs[aa][1] - 1);
      specVal = DirichletBCs[aa][2];

      ind = nn*ndof+dof;
      //specVal  -=  SolnData.var1DotCur[ind];
      specVal  -=  SolnData.var1Cur[ind];

      //cout << start << '\t' << nn << '\t' << ind << '\t' << specVal << endl;

      ind += start1;

      stiff = af*PENALTY/velFact;
      res   = specVal*PENALTY;

      MatSetValue(mtxTemp, ind, ind, stiff, ADD_VALUES);
      VecSetValue(rhsTemp, ind, res, ADD_VALUES);
    }

    return 0;
}



