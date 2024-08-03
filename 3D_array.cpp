#include <iostream>
#include <numeric>
#include <vector>
#include <random>
#include <cstdlib> //rand()
#include <ctime> //time()
#include <algorithm>//random suffle
#include <Eigen/Dense>
#include "DTLZ1_Para.h"
#include "DTLZ1_function.h"
#include "unifrnd.h"

#include <cmath>

using namespace std;
using namespace Eigen;

struct Array_w_idx
{
    MatrixXf Val_;
    MatrixXf idx_;
};
struct empty_Subproblem
{

    MatrixXf lambda;
    MatrixXf Neighbors;
};

void min_array();
void sum_test2();
void min_coeff();
void AssertionEX();
DTLZ1_Para_F header_ex();
void header_ex_func();
void norm_ex();
empty_Subproblem CreateSubProblems();
Eigen::MatrixXf pdist2(const Eigen::MatrixXf& A, const Eigen::MatrixXf& B);
struct Array_w_idx sort_idx(MatrixXf);
void random_NoRepeat();
void AssertionEX2();
void random_in_Range();
int main()
{

//Matrix<int, 3, 3> A;               // Fixed rows and cols. Same as Matrix3d.
//    Matrix<double, 3, Dynamic> B;         // Fixed rows, dynamic cols.
//    Matrix<double, Dynamic, Dynamic> C;   // Full dynamic. Same as MatrixXd.
    Matrix<double, 3, 3, RowMajor> E;     // Row major; default is column-major.
    MatrixXf  Q, R;                     // 3x3 float matrix.
    Vector3f x, y, z;                     // 3x1 float matrix.
    RowVector3f a, b, c;                  // 1x3 float matrix.



//DTLZ1_Para_F AA;
//AA=header_ex();
//AA.nPop=1000;
//cout<<AA.nPop;

DTLZ1_Para_F header_ex();
    random_in_Range();

    return 0;
}

//function

void min_array()
{
    Vector3f A= {4.3,  2.2, 7.11};
    Vector3f B= {5.2, 1.7,  6.3};
    A=A.eval().cwiseMin(B);
    cout<<A;
}

Array_w_idx sort_idx(MatrixXf unsort_arr)
{
    //create array to storage Val & Indx
    Array_w_idx Arr_sorted;
    Arr_sorted.Val_.resize(unsort_arr.rows(),unsort_arr.cols());/**Assertion Error if no initialization**/
    Arr_sorted.idx_.resize(unsort_arr.rows(),unsort_arr.cols());

    //Sort unsort_arr & input to Arr_sorted.Val_ and Arr_sorted.idx_

    for(int j=0; j<unsort_arr.rows(); ++j)
    {
        Eigen::Array<pair<float, int>,1, Eigen::Dynamic >index(1,unsort_arr.cols());//1x6
        for(int i=0; i<unsort_arr.cols(); ++i) //6
        {
            index[i]=make_pair(unsort_arr(j,i),i);
//            cout<<i<< " ,UnSorted, "<<index[i].first<<","<<index[i].second<<endl;
        }
        std::sort(index.begin(),index.end());

        for(int i=0; i<unsort_arr.cols(); ++i)
        {
            Arr_sorted.Val_(j,i)=index[i].first;
            Arr_sorted.idx_(j,i)=index[i].second+1;
        }
    }

//    cout<< "Val="<<Arr_sorted.Val_<<endl;
//    cout<<"idx="<<Arr_sorted.idx_<<endl;
    return Arr_sorted;

};

empty_Subproblem CreateSubProblems()//wrong type. it's array, not "empty_Subproblem"
{

    const int nObj=3, nPop=100, T=20;
//Create Subproblem


    Eigen::Array<empty_Subproblem, 1, Eigen::Dynamic> sp(1,nPop);


    MatrixXf lambda_;
    MatrixXf LAMBDA(nPop,nObj);
    MatrixXf Distance(nPop,nPop);

    for(int i=0; i<nPop; ++i)
    {
        lambda_=MatrixXf::Random(1, nObj).array() * 0.5 + 0.5;/**random type TBC*/
//cout<<"lambda  "<<lambda_<<endl;
//cout<<"lambda_.norm  "<<lambda_.norm()<<endl;
        lambda_=lambda_.array()/lambda_.norm();//Norm­n

        sp[i].lambda=lambda_;
        LAMBDA.row(i)=lambda_;
//        cout<<"sp["<<i<<"].lambda " <<sp[i].lambda<<endl;

    }
//    cout<<"LAMBDA"<<LAMBDA<<endl;
    Distance=pdist2(LAMBDA,LAMBDA);
//    cout<<"dist"<<Distance<<endl;


//    MatrixXf Sorted_Idx(1,nPop);
    for(int i=0; i<nPop; ++i)
    {
        sp[i].Neighbors.resize(1,T);
        Array_w_idx Sorted_Array=sort_idx(Distance.row(i));
//        cout<<i<<" : "<<Sorted_Array.idx_<<endl;
//        cout<<i<<" : "<<Sorted_Array.Val_<<endl;
        sp[i].Neighbors=Sorted_Array.idx_.leftCols(T);

        cout<<i<<" : "<<sp[i].Neighbors<<endl;
        //    Q=MatrixXf::Random(2, 4).array() * 0.5 + 0.5;
//    Array_w_idx AA=sort_idx(Q);


//        [~, Sorted_Idx]=sort_index(D(i,:));//1*100
//        sp.row(i).Neighbors=Sorted_Idx(1:T);//1~20

    }
//    return sp;

}


Eigen::MatrixXf pdist2(const Eigen::MatrixXf& A, const Eigen::MatrixXf& B)
{
    if(B.cols()!=A.cols())cout<<"error: Matrix column are different in pdist2!";


    int A_r= A.rows();
    int B_r = B.rows();
//    int A_c = A.cols();

    MatrixXf V_dist;

    Eigen::MatrixXf C(A_r, B_r);

    for (int i = 0; i < A_r; i++)
    {
        for (int j = 0; j < B_r; j++)
        {
//           dist = 0.0;
//            for (int k = 0; k < d; k++) {
//                float diff = A(i, k) - B(j, k);
//                dist += diff * diff;

            V_dist=(A.row(i)-B.row(j)).array().pow(2);
//            cout<<i<< j<<":"<<V_dist<<endl;
            C(i, j)=(V_dist.sum());
            //.sqrt()
            //(ref_pt_matrix-T_dist).array().pow(2)
//            }
//            C(i, j) = std::sqrt(dist);
        }
    }
//    cout<<C;
    return C;
}
void norm_ex()
{
    Eigen::Vector2d vec(0.3,0.4);
//    float A=vec.norm();
    vec=vec.array()/vec.norm();
    cout<<vec;

}
void header_ex_func()
{
//    MatrixXf xx=MatrixXf::Constant(1,7,0.5);
    MatrixXf xx=unifrnd(0,1,10);
    MatrixXf(*CostFunction)(MatrixXf)=DTLZ_1;

    MatrixXf  DTLZ_F=CostFunction(xx);

    MatrixXf z_min_pt=MatrixXf::Zero(1,3);

    z_min_pt=z_min_pt.eval().cwiseMin(DTLZ_F);

}
DTLZ1_Para_F header_ex()
{
    DTLZ1_Para_F dtlz1_instance;
    int objfun_dim = dtlz1_instance.objfun_dim;
    int nobjfun = dtlz1_instance.nobjfun;
    MatrixXf searchspaceUB = dtlz1_instance.searchspaceUB;
    MatrixXf searchspaceLB = dtlz1_instance.searchspaceLB;
    int MaxIt = dtlz1_instance.MaxIt;
    int nPop = dtlz1_instance.nPop;
    int nbox = dtlz1_instance.nbox;
    string dtlz_name = dtlz1_instance.dtlz_name;
    return dtlz1_instance;

}

void multi_calculation()
{

    MatrixXf temp_rand(1,3);
    temp_rand<<10,20,30;
    MatrixXf UB=MatrixXf::Constant(1,3,1);
    MatrixXf LB=MatrixXf::Constant(1,3,0.5);

    MatrixXf Position_(2,3);//=A.array()-B.array();
    for(int i=0; i<2; ++i)
    {

//        temp_rand=MatrixXf::Random(1,3);
        (Position_.row(i)=(UB.array()-LB.array())*temp_rand.array()).matrix();
        Position_.row(i)=Position_.row(i).array()+LB.array();
    }
//      (P.array() = (A.array() + B.array()) * C.array()).matrix();//calculate plus & * at same code line


    cout<<Position_<<endl;
}

void min_coeff()
{
    Eigen::MatrixXf mat(2,4);

    mat<<2,4,6,8,
        7,7,11,7;

    int r,c;//row idx, col idx
    float s;//min val
    s = mat.minCoeff(&r, &c);

    cout<<s<<", ["<<r<<","<<c<<"]"<<endl;

}
void sum_test2() //from chatgpt
{
    Eigen::MatrixXd mat(3, 3);

    mat << 1, 2, 3,
        4, 5, 6,
        7, 8, 9;
    mat=mat.array().pow(2);

    Eigen::MatrixXd mat2= mat.rowwise().sum(); //only this one 'll be good.
    mat.conservativeResize(3, 1);
    mat=mat2;
    cout<<mat;


//    MatrixXf NPOP_mat=MatrixXf::Constant(99,1,1);
//    MatrixXf pop(1,3);
//    pop<<8.8, 9.9, 4.4;
//    NPOP_mat=NPOP_mat.eval()*pop;

//    std::cout << NPOP_mat << std::endl;

}
void sum_test1()
{
    MatrixXf A;
    A.setOnes(3,4);

    Eigen::MatrixXf W(4,2);
    W.fill(3);

    MatrixXf AW;
    MatrixXf AX,AZ;
    AW=A*W;
    AW=AW.array().pow(2);

    AX=AW.colwise().sum();
    AZ=AW.rowwise().sum();

    cout<<"AW"<<AW<<endl;
    cout<<"AZ"<<AZ<<endl;
    cout<<"AX"<<AX;
}
void case_1()
{

    Eigen::Matrix<float, 3, 3> n_a,n_b,n_c,n_d;

    Eigen::MatrixXf A(6,6);
    n_a.fill(1);
    n_b.fill(5);
    n_c.fill(9);
    n_d.fill(27);
    A<<n_a,n_b,n_c,n_d;


}

void random_NoRepeat()
{
    Matrix<int, 1, Eigen::Dynamic> vec(20);
    vec<<1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        11,12,13,14,15,16,17,18,19,20;


    cout << "list  : ";
    for (int i = 0; i < vec.size(); i++)
    {
        cout << vec[i] << " ";
    }
    cout << endl;

    //std::srand(unsigned(std::time(0)));
//    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//    std::srand(seed);

    // 方法1 沒使用std::srand的話每次亂數結果都一樣
//    std::random_shuffle(vec.begin(), vec.end()); // using built-in random generator

    std::random_device rd;
    std::shuffle(vec.begin(), vec.end(), std::default_random_engine(rd()));


    cout << "result: ";
    for (int i = 0; i < vec.size(); i++)
    {
        cout << vec[i] << " ";
    }
    cout << endl;
}

void random_in_Range()
{
    srand(time(NULL));
    int random_num;
    int VarMin=0,
        VarMax=20;

    for(int i=0; i<5; ++i)
    {
        random_num = rand()%(VarMax-VarMin+1)+VarMin;
        cout<<i<<" : "<<random_num<<endl;
    }

}
void AssertionEX()
{
//    Eigen::MatrixXf P;
//    P.setRandom(2,4);
    MatrixXf A;
    A=MatrixXf::Constant(1,3,1.5);

    MatrixXf B=MatrixXf::Constant(1,3,1);

    MatrixXf C=MatrixXf::Constant(1,3,5);
    MatrixXf P;//=A.array()-B.array();

    (P.array() = (A.array() + B.array()) * C.array()).matrix();//calculate plus & * at same code line



//    Eigen::MatrixXf C;
//    C=P.middleCols<2>(1);
//
//    cout<<C<<endl;
}

void AssertionEX2()
{
    int nPop=10;
    MatrixXf Cost_mat_origin;
    Cost_mat_origin=MatrixXf::Random(nPop, 3).array() * 0.5 + 0.5;
    MatrixXf Cost_mat(nPop-1,3);

    int j=2;
    Cost_mat.topRows(j-1)=Cost_mat_origin.topRows(j-1);
    Cost_mat.bottomRows(nPop-j)=Cost_mat_origin.bottomRows(nPop-j);

    for(int i=0; i<nPop-1; ++i)
    {

        cout<<i<<" || original:"<<Cost_mat_origin.row(i)<<", rematch:"<<Cost_mat.row(i)<<endl;
//if(i%2==0)cout<<endl;
    }
    cout<<nPop<<" || original:"<<Cost_mat_origin.row(nPop-1);
}
void transposeEX()
{

    Eigen::MatrixXf Ref_Mat(3,136);
    Eigen::MatrixXf W;
    int count=0;
    for(int e=15; e>=0; --e)
    {
        W.resize(e+1,3);
        for(int i=0; i<=e; ++i)
        {
            W.row(i)<<15-e,i,e-i;
        }
        W.transposeInPlace() ;

        cout<<e<<endl;
        cout<<W<<endl;

        Ref_Mat.middleCols(count,e+1)=W;

        cout<<Ref_Mat<<endl;
        count+=(e+1);

        //A.transposeInPlace() ;
        //A/=2;


    }
    cout<<Ref_Mat;
}

