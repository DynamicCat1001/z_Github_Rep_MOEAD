#include <iostream>
#include <vector>
#include <Eigen/Dense>



using namespace std;
using namespace Eigen;

Eigen::MatrixXf GenerateReferencePoints(int, int);
Eigen::MatrixXf GetFixedRowSumIntegerMatrix(int num_obj, int num_obj_dim_depart);


int main()//GenerateReferencePoints
{
    int num_obj=3;
    int num_obj_dim_depart=15;
    Eigen::MatrixXf matrices;

    matrices=GenerateReferencePoints(num_obj,num_obj_dim_depart);

}

Eigen::MatrixXf  GenerateReferencePoints(int num_obj, int num_obj_dim_depart)
{

    Eigen::MatrixXf Ref_Mat;


    Ref_Mat=GetFixedRowSumIntegerMatrix(num_obj,num_obj_dim_depart);

    Ref_Mat=Ref_Mat/num_obj_dim_depart;

    cout<<Ref_Mat;
    return Ref_Mat;

}

Eigen::MatrixXf  GetFixedRowSumIntegerMatrix(int num_obj, int RowSum)
{

    int total=(2+RowSum)*(RowSum+1)/2;
    Eigen::MatrixXf Ref_Mat(num_obj, total);
    int count=0;

    for(int e=RowSum; e>=0; --e)
    {

        Eigen::MatrixXf W(e+1,num_obj);

        for(int i=0; i<=e; ++i)
        {
            W.row(i)<<RowSum-e,i,e-i;
        }
        W.transposeInPlace();


        Ref_Mat.middleCols(count,e+1)=W;
        count+=(e+1);


    }
//    cout<<Ref_Mat;
    return Ref_Mat;

}

/*GetFixedRowSum(Matrix, 3, 15)

    for(0:15)
    i=0
        B=GetFixedRowSum(2,15)
                for(0:15)
                i=0
                    B=GetFixedRowSum(1,15)
                        A=15; return;
                    B=15
                    A={ , {0, 15}} => A ={0,15}
                i=1
                    B=GetFixedRowSum(1,14)
                        A=14; return;
                    B=14
                    A={{0,15},{1,14}}
                i=2:15....
                    A={{0,15}, {1,14}, {2,13}, {3,12},{4,11},{},{}{}{},{15,0}};
                    return A;
        B[16][2]={
        {0,15}, {1,14}, {2,13}, {3,12}, {4,11},
        {5,10}, {6,9}, {7,8}, {8,7}, {9,6},
        {10,5}, {11,4}, {12,3}, {13,2}, {14,1}, {15,0}
        }
        float temp[size(B)]=temp[16]={0}
        A={ ; 0*temp[16][1] B[16][2]} => A[16][3] //resize()??
        A0={
        0 0 15
        0 1 14
        0 2 13
        0 3,12
        0 4,11
        0 5,10
        0 6,9
        0 7,8
        0 8,7
        ...
        0,15,0
        }
        return A;
    i=1
       B=GetFixedRowSum(2,14)
            i=0
                B=GetFixedRowSum(1,14)
                    A=14
                B=14
                A={ , {0,14}}
            ...
            i=14
                A[15][2]={{0,14},{1,13}, ... ,{13,1},{14,0}}
                return A;
        B[15][2]={}
        float temp[size(B)]=temp[15]={1}
        A={ A[16][3]; 1*temp[15][1] B[15][2]} => A={ A0[16][3]; A1[15][3]}
        A1={
        1 0 14
        1 1 13
        1 2 12
        1 3 13
        1 4 14
        ...
        1 14 0
        }

    ...
    i=14

    A14[2][3]={
    14 0 1
    14 1 0
    }


    i=15
        B=GetFixedRowSum(2,15-15)
            for(0:0)
                i=0
                    B=GetFixedRowSum(1,0)
                        A=0; return;
                    B=0
                    A={ , {0, 0}} => A ={0,0}
                    return A;
        B={0,0}
                float temp[size(B)]=temp[1]={1}
        A={ ; 15*1 {0,0}} => A[3]={15,0,0} //resize()??
        A15={15,0,0}
        return A;


        A[][][] = {A[16][3]; temp[15][1] B[15][2]}

        =>

        A={ A0[16][3], A1[15][3], ..., A15[1][3]


*/






