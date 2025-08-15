#ifndef ATLASQTL_TYPES_H_
#define ATLASQTL_TYPES_H_

#include <RcppEigen.h>

// These typedefs have to be in a separate header file to be properly copied to RcppExports

typedef Eigen::ArrayXd Arr1D;
typedef Eigen::ArrayXXd Arr2D;
typedef Eigen::Map<Eigen::ArrayXd> MapArr1D;
typedef Eigen::Map<Eigen::ArrayXXd> MapArr2D;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;


#endif // ATLASQTL_TYPES_H_
