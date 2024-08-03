#ifndef TEST_HEAD_H_INCLUDED
#define TEST_HEAD_H_INCLUDED

#include <iostream>
#include <fstream>
#include <random>
#include <algorithm> //for std::shuffle as randsample in matlab
#include <cstdlib> //rand()
#include <ctime> //time()
//#include <conio.h>
#include <Eigen/Dense>
#include <type_traits>//for isfield>>has_cost
#include <memory>//std::unique_ptr

#include "DTLZ1_Para.h"
#include "DTLZ1_function.h"
#include "unifrnd.h"


using namespace std;
using namespace Eigen;

Eigen::MatrixXf  GetFixedRowSumIntegerMatrix(int, int);



#endif // TEST_HEAD_H_INCLUDED

Eigen::MatrixXf Read_file1(int col_, int row_)//col=2500, row=3
{

    ifstream file2("DTLZ1ReferencePoints.dat");
    if (!file2)
    {
        cout << "file can't open" <<endl;
    }
    else
    {
        Eigen::MatrixXf p(col_,row_);
        Eigen::RowVector3f row_element;
        float A;
        for(int i=0; i<col_; ++i)//2500
        {
            for(int j=0; j<row_; ++j)//3 dim.
            {
                A=0;
                file2>>A;
                row_element[j]=A;
            }
            p.row(i)=row_element;
        }


        cout<<"data read finished"<<endl;

        file2.close();
        return p;
    }
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
