
#include "HBSplineFEM.h"
#include "MpapTime.h"
#include "ImmersedIntegrationElement.h"
#include "myDataIntegrateCutFEM.h"
#include "QuadratureUtil.h"
#include "BasisFunctionsLagrange.h"

extern MpapTime mpapTime;


void  HBSplineFEM::computeTotalForce(int index)
{
  return;
}



void  HBSplineFEM::computeTotalBodyForce(int index)
{
  double val, forceTotal[3];
  int dir, ee;
   
  for(dir=0;dir<DIM;dir++)
  {
    val = 0.0;
    for(ee=0;ee<activeElements.size();ee++)
      val += elems[activeElements[ee]]->computeTotalBodyForce(index, dir);

    forceTotal[dir] = val;
  }

  char        tmp[200];
  MyString    tmpStr;
    
  sprintf(tmp," \t %12.6E \t %12.6E ", forceTotal[0], forceTotal[1]);
  tmpStr.append(tmp);

  prgWriteToTFile(tmpStr);

  return;
}




void  HBSplineFEM::solveSolidProblem()
{
  IB_MOVED = false;
  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
  {
    //cout << " Lagrange multipliers " << endl;        printVector(&(FluidSolnData.var3(0)), IBDOF);
    //ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));

    //cout << " kkkkkkkkkkk " << endl;
    ImmersedBodyObjects[bb]->solveTimeStep();

    IB_MOVED = (IB_MOVED || ImmersedBodyObjects[bb]->updatePointPositions() );
  }

  return;
}




void  HBSplineFEM::writeNodalData()
{
  writeFluidOutput();
  writeImmersedSolidOutput();

  return;
}


void  HBSplineFEM::writeFluidOutput()
{
  return;
}



void  HBSplineFEM::writeImmersedSolidOutput()
{
  for(int bb=0; bb<nImmSolids; bb++)
    ImmersedBodyObjects[bb]->writeOutput();

  return;
}



void  HBSplineFEM::printResultAtPoint(int index, double u1, double v1, double w1)
{
  return;
}



void  HBSplineFEM::applyInterfaceTerms2D()
{
    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed boundary points
    // Lagrange multipliers
    //  or
    // Penalty method
    //
    ////////////////////////////////////////////////////////

    int ii, jj, aa, bb, c, r, gp, nr;

    MatrixXd  Klocal, Klocal2;
    VectorXd  Flocal;
    node *nd;

    double PENALTY;
    ImmersedIntegrationElement  *lme;
    
    int nlf=(degree[0]+1)*(degree[0]+1);

    int nlb, ind1, ind2, nGauss;

      VectorXd  NN(nlf), dNN_dx(nlf), dNN_dy(nlf), dN_dx, dN_dy, Nf;
      VectorXd  Flocal2, vel(DIM), vel2(DIM), lagmults(DIM), Nb, dNb, xx, yy,  specValx, specValy, res(DIM);
      MatrixXd  Khorz;
      myPoint  knotIncr, knotBegin, knotEnd;

      double  detJ, af, dvol, fact, fact1, fact2;

      af = SolnData.td(2);
      
      bool axsy = (GeomData.FluidProps[2] == 1);
      double  rho = GeomData.FluidProps[3];
      double  mu = GeomData.FluidProps[4];


      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        nlb = ImmersedBodyObjects[0]->ImmIntgElems[0]->pointNums.size();

        Nb.resize(nlb);
        dNb.resize(nlb);
        xx.resize(nlb);
        yy.resize(nlb);
        specValx.resize(nlb);
        specValy.resize(nlb);

        if(ImmersedBodyObjects[bb]->isBoundaryConditionTypeLagrange())
        {
          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][1];

              specValx[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][0];
              specValy[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][1];
            }

            //cout << specValx[0] << '\t' << specValy[0] << endl;
            //cout << specValx[1] << '\t' << specValy[1] << endl;

            //cout << " lme->gausspoints.size() = " << lme->gausspoints.size() << endl;

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine2D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << lme->gaussweights[gp] << '\t' << Nb[0] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              lagmults.setZero();
              geom.setZero();
              vel2.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0]     += Nb[ii] * xx[ii];
                geom[1]     += Nb[ii] * yy[ii];

                vel2[0]     += Nb[ii] * specValx[ii];
                vel2[1]     += Nb[ii] * specValy[ii];

                lagmults[0] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM);
                lagmults[1] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM+1);
              }

              if(axsy)
              {
                dvol *= 2.0*PI*geom[0];
              }

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel2[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d \n", Nb[0], nlb);
              //printf("vel2[0] = %12.6f, vel2[1] = %12.6f \n", vel2[0], vel2[1]);
              //printf("lagmults[0] = %12.6f, lagmults[1] = %12.6f \n", lagmults[0], lagmults[1]);
              //printf("detJ = %12.6f, gw = %12.6f \n", detJ, lme->gaussweights[gp]);

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              knotBegin = nd->getKnotBegin();
              knotEnd   = nd->getKnotEnd();
              knotIncr  = nd->getKnotIncrement();

              GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

              if(nd->getParent() == NULL )
              {
                Nf = NN;
                dN_dx = dNN_dx;
                dN_dy = dNN_dy;
              }
              else
              {
                Nf    = nd->SubDivMat * NN;
                dN_dx = nd->SubDivMat * dNN_dx;
                dN_dy = nd->SubDivMat * dNN_dy;
              }

              //printVector(Nf);
              //printf("\n\n");

              ind1 = lme->pointNums.size();
              ind2 = nd->GlobalBasisFuncs.size();

              Khorz.resize(ind1*DIM, ind2*ndof);    Khorz.setZero();

              Flocal.resize(ind2*ndof);   Flocal.setZero();
              Flocal2.resize(ind1*DIM);   Flocal2.setZero();

              vel(0) = vel2(0) - nd->computeValueCur(0, Nf);
              vel(1) = vel2(1) - nd->computeValueCur(1, Nf);

              for(ii=0;ii<ind1;ii++)
              {
                r = DIM*ii;
                fact1 = Nb[ii]*dvol;

                Flocal2(r)   += fact1*vel(0);
                Flocal2(r+1) += fact1*vel(1);

                fact1 = af*fact1;
                //cout << ii << '\t' << r << '\t' << fact1 << endl;

                for(jj=0;jj<ind2;jj++)
                {
                  c = ndof*jj;

                  fact2 = fact1*Nf[jj];
                  //cout << jj << '\t' << c << '\t' << fact2 << endl;
                  Khorz(r,   c)   += fact2;
                  Khorz(r+1, c+1) += fact2;
                }
              }

              for(ii=0;ii<ind2;ii++)
              {
                r = ndof*ii;

                fact = Nf[ii]*dvol;

                Flocal(r)   -= fact*lagmults[0];
                Flocal(r+1) -= fact*lagmults[1];
              }

              //printMatrix(Khorz);
              //printf("\n\n");
              //printVector(Flocal);
              //printf("\n\n");
              //printVector(Flocal3);
              //printf("\n\n");

              ind1 = lme->posIndices.size();
              ind2 = nd->forAssyVec.size();

              //printVector(lme->posIndices);
              //printf("\n\n");
              //printVector(nd->forAssyVec);
              //printf("\n\n");

              for(ii=0;ii<ind2;ii++)
              {
                c = nd->forAssyVec[ii];

                solverEigen->rhsVec[c] += Flocal(ii) ;
              }

              for(ii=0;ii<ind1;ii++)
              {
                r = fluidDOF + lme->posIndices[ii];

                solverEigen->rhsVec[r] += Flocal2(ii);

                for(jj=0;jj<ind2;jj++)
                {
                  c = nd->forAssyVec[jj];

                  fact = Khorz(ii, jj);

                  solverEigen->mtx.coeffRef(r, c) += fact;
                  solverEigen->mtx.coeffRef(c, r) += fact;
                }
              }

            }//for(gp=0...
          }//for(aa=0...
        }//if(
        else
        {
          myDataIntegrateCutFEM  myData;

          PENALTY = ImmersedBodyObjects[bb]->getPenaltyParameter();

          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][1];
            }

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine2D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << gaussweights[gp] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              vel.setZero();

              geom.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0] += Nb[ii] * xx[ii];
                geom[1] += Nb[ii] * yy[ii];
              }

              if(axsy)
              {
                dvol *= 2.0*PI*geom[0];
              }

              //cout << " uuuuuuuuuuu " << endl;

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              //cout << " uuuuuuuuuuu " << endl;
              //cout << " dvol " << dvol << endl;

              nr = nd->forAssyVec.size();

              myData.K1 = MatrixXd::Zero(nr, nr);
              myData.F1 = VectorXd::Zero(nr);

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //if(geom[1] == 1.0 )
                //vel[0] = 1.0;

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d, ee = %5d \n", Nb[0], nlb, lme->elemNums[gp]);

              //cout << " uuuuuuuuuuu " << endl;

              myData.geom  = geom;
              myData.param = param;
              myData.dvol  = PENALTY * dvol;

              for(jj=0;jj<DIM;jj++)
              {
                myData.dir = jj;
                myData.specVal[0] = vel[jj];

                //nd->applyBoundaryConditionsAtApoint(jj, param, vel[jj], fact, myData.K1, myData.F1);
                nd->applyBoundaryConditionsAtApoint(myData);
              }

              //solverEigen->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, nd->forAssyVec2, Klocal, Flocal);
              solverEigen->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, myData.K1, myData.F1);
              //cout << " uuuuuuuuuuu " << endl;
            }//for(gp=0...
          }//for(aa=0...
        }//else
      }//for(bb=0;...
      //

  return;
}




void  HBSplineFEM::applyInterfaceTerms3D()
{
    int ii, jj, aa, bb, c, r, gp, nr;

    MatrixXd  Klocal, Klocal2;
    VectorXd  Flocal;
    node *nd;

    double PENALTY;
    ImmersedIntegrationElement  *lme;
    
    int nlf=(degree[0]+1)*(degree[0]+1)*(degree[0]+1);
    int nlb, ind1, ind2, nGauss;

      VectorXd  NN(nlf), dNN_dx(nlf), dNN_dy(nlf), dNN_dz(nlf), dN_dx, dN_dy, dN_dz, Nf;
      VectorXd  Flocal2, vel(DIM), vel2(DIM), lagmults(DIM), Nb, dNb, res(DIM);
      VectorXd  xx, yy, zz, specValx, specValy, specValz;
      MatrixXd  Khorz;
      myPoint  knotIncr, knotBegin;

      double  detJ, af, dvol, fact, fact1, fact2;

      af = SolnData.td(2);
      double  rho = GeomData.FluidProps[3];
      double  mu = GeomData.FluidProps[4];

      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        nlb = ImmersedBodyObjects[0]->ImmIntgElems[0]->pointNums.size();

        Nb.resize(nlb);
        dNb.resize(nlb);
        xx.resize(nlb);
        yy.resize(nlb);
        zz.resize(nlb);
        specValx.resize(nlb);
        specValy.resize(nlb);
        specValz.resize(nlb);

        if(ImmersedBodyObjects[bb]->isBoundaryConditionTypeLagrange())
        {
          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][1];
              zz[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][2];
              
              specValx[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][0];
              specValy[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][1];
              specValz[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][2];
            }

            //cout << " uuuuuuuuuuu " << endl;
            //cout << specValx[0] << '\t' << specValy[0] << endl;
            //cout << specValy[0] << '\t' << specValy[1] << endl;

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine3D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &zz(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << lme->gaussweights[gp] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              lagmults.setZero();
              geom.setZero();
              vel2.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0]     += Nb[ii] * xx[ii];
                geom[1]     += Nb[ii] * yy[ii];
                geom[2]     += Nb[ii] * zz[ii];
                
                vel2[0]     += Nb[ii] * specValx[ii];
                vel2[1]     += Nb[ii] * specValy[ii];
                vel2[2]     += Nb[ii] * specValz[ii];

                lagmults[0] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM);
                lagmults[1] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM+1);
                lagmults[2] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM+2);
              }

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel2[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d \n", Nb[0], nlb);
              //printf("lagmults[0] = %12.6f, lagmults[1] = %12.6f \n", lagmults[0], lagmults[1]);
              //printf("detJ = %12.6f, gw = %5d \n", detJ, lme->gaussweights[gp]);

              //lme->computePointAtGP(gp, geom);

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              knotBegin = nd->getKnotBegin();
              knotIncr  = nd->getKnotIncrement();

              GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

              if(nd->getParent() == NULL )
              {
                Nf = NN;
                dN_dx = dNN_dx;
                dN_dy = dNN_dy;
                dN_dz = dNN_dz;
              }
              else
              {
                Nf    = nd->SubDivMat * NN;
                dN_dx = nd->SubDivMat * dNN_dx;
                dN_dy = nd->SubDivMat * dNN_dy;
                dN_dz = nd->SubDivMat * dNN_dz;
              }

              //printVector(Nf);
              //printf("\n\n");

              ind1 = lme->pointNums.size();
              ind2 = nd->GlobalBasisFuncs.size();

              //cout << " ind1 = " << ind1 << endl;
              //cout << " ind2 = " << ind2 << endl;

              Khorz.resize(ind1*DIM, ind2*ndof);    Khorz.setZero();

              Flocal.resize(ind2*ndof);   Flocal.setZero();
              Flocal2.resize(ind1*DIM);   Flocal2.setZero();

              vel(0) = vel2(0) - nd->computeValueCur(0, Nf);
              vel(1) = vel2(1) - nd->computeValueCur(1, Nf);
              vel(2) = vel2(2) - nd->computeValueCur(2, Nf);

              for(ii=0;ii<ind1;ii++)
              {
                r = DIM*ii;
                fact1 = Nb[ii]*dvol;

                Flocal2(r)   += fact1*vel(0);
                Flocal2(r+1) += fact1*vel(1);
                Flocal2(r+2) += fact1*vel(2);

                fact1 = af*fact1;
                //cout << ii << '\t' << r << '\t' << fact1 << endl;

                for(jj=0;jj<ind2;jj++)
                {
                  c = ndof*jj;

                  fact2 = fact1*Nf[jj];
                  //cout << jj << '\t' << c << '\t' << fact2 << endl;
                  Khorz(r,   c)   += fact2;
                  Khorz(r+1, c+1) += fact2;
                  Khorz(r+2, c+2) += fact2;
                }
              }

              for(ii=0;ii<ind2;ii++)
              {
                r = ndof*ii;

                fact = Nf[ii]*dvol;

                Flocal(r)   -= fact*lagmults[0];
                Flocal(r+1) -= fact*lagmults[1];
                Flocal(r+2) -= fact*lagmults[2];
              }

              //cout << " assembling " << endl;
              //printMatrix(Khorz);
              //printf("\n\n");
              //printVector(Flocal);
              //printf("\n\n");
              //printVector(Flocal3);
              //printf("\n\n");

              ind1 = lme->posIndices.size();
              ind2 = nd->forAssyVec.size();

              //printVector(lme->posIndices);
              //printf("\n\n");
              //printVector(nd->forAssyVec);
              //printf("\n\n");

              //cout << " ind1 = " << ind1 << endl;
              //cout << " ind2 = " << ind2 << endl;

              for(ii=0;ii<ind2;ii++)
              {
                c = nd->forAssyVec[ii];
                
                solverEigen->rhsVec[c] += Flocal(ii) ;
              }

              for(ii=0;ii<ind1;ii++)
              {
                r = fluidDOF + lme->posIndices[ii];
                
                solverEigen->rhsVec[r] += Flocal2(ii);

                for(jj=0;jj<ind2;jj++)
                {
                  c = nd->forAssyVec[jj];
                  
                  fact = Khorz(ii, jj);
                  
                  //cout << ii << '\t' << r << '\t' << jj << '\t' << c << endl;

                  solverEigen->mtx.coeffRef(r, c) += fact;
                  solverEigen->mtx.coeffRef(c, r) += fact;
                }
              }

            }//for(gp=0...
          }//for(aa=0...
        }//if(
        else
        {
          myDataIntegrateCutFEM  myData;
          
          PENALTY = ImmersedBodyObjects[bb]->getPenaltyParameter();

          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][1];
              zz[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][2];
            }

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine3D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &zz(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << gaussweights[gp] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              vel.setZero();

              geom.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0] += Nb[ii] * xx[ii];
                geom[1] += Nb[ii] * yy[ii];
                geom[2] += Nb[ii] * zz[ii];
              }

              //cout << " uuuuuuuuuuu " << endl;

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              //cout << " uuuuuuuuuuu " << endl;
              //cout << " dvol " << dvol << endl;

              nr = nd->forAssyVec.size();

              myData.K1 = MatrixXd::Zero(nr, nr);
              myData.F1 = VectorXd::Zero(nr);

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //if(geom[1] == 1.0 )
                //vel[0] = 1.0;

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d, ee = %5d \n", Nb[0], nlb, lme->elemNums[gp]);

              //cout << " uuuuuuuuuuu " << endl;

              myData.geom  = geom;
              myData.param = param;
              myData.dvol  = PENALTY * dvol;

              for(jj=0;jj<DIM;jj++)
              {
                myData.dir = jj;
                myData.specVal[0] = vel[jj];

                //nd->applyBoundaryConditionsAtApoint(jj, param, vel[jj], fact, myData.K1, myData.F1);
                nd->applyBoundaryConditionsAtApoint(myData);
              }

              //solverEigen->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, nd->forAssyVec2, Klocal, Flocal);
              solverEigen->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, myData.K1, myData.F1);
              //cout << " uuuuuuuuuuu " << endl;
            }//for(gp=0...
          }//for(aa=0...
        }//else
      }//for(bb=0;...


  return;
}















