#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "DTLZ1_function.h"
#include "DTLZ1_Para.h"
using namespace std;
using namespace Eigen;

int main(){

    Matrix<float,1,7> pos;
//    pos=MatrixXf::Random(1, 7).array() * 0.5 + 0.5;
pos.fill(0);
    cout<<pos<<endl;
    Matrix<float,1,3> cost_;
    cost_=DTLZ_1(pos);
    cout<<cost_;
    return 0;
}
