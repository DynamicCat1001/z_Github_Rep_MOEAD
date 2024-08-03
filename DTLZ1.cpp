/**DTZL1 DTLZ1 multi-objective function
%   This function represents a hyper-plane.
%   Using k = 5, the number of dimensions must be n = (M - 1) + k.
%   The Pareto front of this function happens at xm := x(n-k+1:end) = 0.5,
%   that is, when the last k variables are all equal to 0.5. The first
%   M-1 variables are varied to map the whole Pareto front. Of course, this
%   is not the easiest thing to do for more than M = 3 objectives. However,
%   the hyper-surface is a hyper-plane satisfying the equation
%      f1 + f2 + f3 + ... + fm = 0.5
%
%   Syntax:
%      f = dtzl1(x, M)
%
%   Input arguments:
%      x: a n x mu matrix with mu points and n dimensions
%      M: a scalar with the number of objectives
%
%   Output argument:
%      f: a m x mu matrix with mu points and their m objectives computed at
%         the input
*/
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
# define M_PI           3.14159265358979323846

MatrixXf DTLZ_1(MatrixXf x);
void disp(MatrixXf(*f)(MatrixXf));


int main()//(MatrixXf x)
{
//    MatrixXf x=MatrixXf::Constant(1,10,0.5);
    disp(DTLZ_1);
//    DTLZ_1;
    return 0;
}

void disp(MatrixXf (*f)(MatrixXf)){
    MatrixXf x=MatrixXf::Constant(1,10,0.5);
    (*f)(x);
    cout<<"display!"<<endl;
}
MatrixXf DTLZ_1(MatrixXf x){
//    MatrixXf x=MatrixXf::Random(1,10);
//    MatrixXf x=MatrixXf::Constant(1,10,0.5);


    float pi_=M_PI;
    int M = 3,
        k = 5; //as suggested by Deb
    int n = M-1 + k; // Error check: the number of dimensions must be M-1+k

    MatrixXf Mat_f(1,3), Mat_g;
    MatrixXf xm;
    float const_g;

    xm = x.rightCols(k);//xm contains the last k variables
//    cout<<xm<<endl;

    Mat_g=pow((xm.array() - 0.5),2)-(20*pi_*(xm.array() - 0.5)).cos();//3*x.size
    const_g=100*(Mat_g.sum()+k);//3*1

    /* Now, computes the functions. The first and the last will be written separately to facilitate things*/

    Mat_f.col(0)=0.5*(1+const_g)*(x.col(0)*x.col(1));

    for(int i=2;i<M;++i){//M=5,M-2=3,M-3=2,M-4=1
        Mat_f(i-1)=0.5*(1+const_g)*(x.leftCols(M-i).prod())*(1-x(M-i));
        //x.leftCols(M-i)==x(1:1)
    }
    Mat_f(M-1)=0.5*(1-x(0))*(1+const_g);
    cout<<Mat_f<<endl;
//    cout<<Mat_f.sum()<<endl;

    return Mat_f;
    /* best sol will be 0.5 sum hyperplane
    f(1)=0.5* prod(x(1:2))*(1+g)
    // f(2)=0.5* prod(x(1:1))*(1-x(3-2+1))*(1+g)
    // f(3)=0.5*(1-x(1))*(1+g)
    */
}










