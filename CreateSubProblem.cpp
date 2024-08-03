#include <iostream>
#include <numeric>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include "DTLZ1_Para.h"
#include "DTLZ1_function.h"

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
empty_Subproblem CreateSubProblems();
Eigen::MatrixXf pdist2(const Eigen::MatrixXf& A, const Eigen::MatrixXf& B);
struct Array_w_idx sort_idx(MatrixXf);

int main()
{

//    Q=MatrixXf::Random(2, 4).array() * 0.5 + 0.5;
//    Array_w_idx AA=sort_idx(Q);
    empty_Subproblem sp=CreateSubProblems();
}

//function

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

empty_Subproblem CreateSubProblems()
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



void header_ex_func()
{
    MatrixXf xx=MatrixXf::Constant(1,10,0.5);
    MatrixXf(*CostFunction)(MatrixXf)=DTLZ_1;

    MatrixXf  DTLZ_F=CostFunction(xx);
    cout<<DTLZ_F;
}

/*
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

}*/


