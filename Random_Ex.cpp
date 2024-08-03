#include <iostream>
//#include <numeric>

#include <random>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;


void best_random();
void uniform_real_distribution_EX();

int main()
{
    best_random();

}
//func.
MatrixXf unifrnd(float VarMin, float VarMax, int VarSize)
{


    MatrixXf random_num(1,VarSize);

    std::random_device rd;
    std::mt19937 generator( rd() );
    /* normal distribution*/
    std::uniform_real_distribution<float> distribution(VarMin, VarMax);

    /* generate normal distribution random nums */
    for (int i=0;i<VarSize;++i){
       random_num(i) = distribution(generator);
    }

    std::cout << "position = " << random_num << endl;
    return random_num;

}


void best_random()
{
    float VarMin=0;
    float VarMax=10;
    int VarSize=10;
    MatrixXf position(1,VarSize);

    std::random_device rd;
    std::mt19937 generator( rd() );
    /* 標準常態分布 */
    std::uniform_real_distribution<float> distribution(VarMin, VarMax);

    /* 產生標準常態分布的亂數 */
    for (int i=0;i<VarSize;++i){
       position(i) = distribution(generator);
    }

    std::cout << "position = " << position << endl;

}
void uniform_real_distribution_EX()
{
    std::uniform_real_distribution<float> unif(0.0, 6.0);



//    Eigen::Matrix<double, 3, 4>::Random().array() * 0.5 + 0.5
    MatrixXf X=MatrixXf::Random(3,4).array() * 0.5 + 0.5;

//    MatrixXf::
    cout<<X;
}
