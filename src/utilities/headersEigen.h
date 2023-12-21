/**
 * @file
 * @brief Header which adds support for the Eigen vectors/matrices to our code.
 *
 */

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_SUPERLU_SUPPORT

#ifndef incl_headersEigen_h
#define incl_headersEigen_h

#undef Sucess

#include <Eigen/Core>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeLinearSolvers>

using namespace Eigen;


//
// Need the sparse matrix to be row-major to be compatible with PARDISO solver.
// 
typedef  SparseMatrix<double, RowMajor>  SparseMatrixXd;
typedef  SparseMatrix<double>  SparseMatrixXdCM;
//typedef  SparseMatrix<double>  SparseMatrixXd;
typedef  SparseMatrix<float>   SparseMatrixXf;
typedef  SparseMatrix<int>     SparseMatrixXi;

typedef  DynamicSparseMatrix<double>  DynamicSparseMatrixXd;
typedef  DynamicSparseMatrix<float>   DynamicSparseMatrixXf;
typedef  DynamicSparseMatrix<int>     DynamicSparseMatrixXi;

typedef  Matrix<double, Dynamic, Dynamic, RowMajor>  MatrixXdRM;

typedef  Eigen::Vector3d  myPoint;


#endif
