#include <igl/boundary_loop.h>
#include <igl/cat.h>
#include <igl/colon.h>
#include <igl/slice.h>
#include <igl/upsample.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/harmonic.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>

using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using namespace igl;

typedef VectorXi Integer_Vector;
typedef VectorXd Double_Vector;
typedef MatrixXd Double_Matrix;
typedef MatrixXi Integer_Matrix;