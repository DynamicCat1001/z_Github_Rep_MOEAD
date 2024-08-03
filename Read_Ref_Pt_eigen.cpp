#include <iostream>
#include <fstream>
#include <conio.h>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

Eigen::MatrixXf Read_file(int col_, int row_);//col=2500, row=3
Eigen::MatrixXf Slope(Eigen::MatrixXf Zr_matrix, Eigen::MatrixXf ref_pt_matrix);


int main()
{

    Eigen::MatrixXf Zr_matrix = GenerateReferencePoints(3, 15);//3x136

    Eigen::MatrixXf ref_pt_matrix=Read_file(2500,3);
    Eigen::MatrixXf s=Slope(Zr_matrix,ref_pt_matrix,15);

//    cout<<ref_pt_matrix.rows()<<","<<ref_pt_matrix.cols();
//    cout<<" ,matrix"<<ref_pt_matrix<<endl;

    return 0;
//    getchar();

}

Eigen::MatrixXf Slope(Eigen::MatrixXf Zr_matrix, Eigen::MatrixXf ref_pt_matrix, int psdiv)
{
    Eigen::MatrixXf temp_2,slope_0,slope_1,slope_2,s_dist,min_matrix;
    int min_idxR,min_idxC;
    float min_val;

    int ps=0;
    float div_=1.00f/psdiv;
    Eigen::MatrixXf rps;
    Eigen::MatrixXf T_dist(ref_pt_matrix.rows(),ref_pt_matrix.cols());

    for(int i=0; i<Zr_matrix.cols(); ++i)
    {
        temp_2=(Zr_matrix.col(i)).transpose();//temp_2.transposeInPlace();
        slope_0=MatrixXf::Constant(ref_pt_matrix.rows(), 1,1) * (Zr_matrix.col(i)).transpose();//2500x3

        slope_1=slope_0.array()*ref_pt_matrix.array();
        slope_1 = slope_1.rowwise().sum().eval();


        slope_2=slope_0.array().square();
        slope_2 =slope_2.rowwise().sum().eval();


        s_dist=slope_1.array()/slope_2.array();

        T_dist<<s_dist,s_dist,s_dist;
        T_dist=T_dist.array()*slope_0.array();

        min_matrix= (ref_pt_matrix-T_dist).array().pow(2);
        min_matrix=(min_matrix.rowwise().sum().eval()).array().sqrt();//sqrt(sum((reference_ps-t).^2,2));

        min_val = min_matrix.minCoeff(&min_idxR,&min_idxC);
//        cout<<"1_over_15: "<<div_<<", min val:"<<min_val<<endl;

        if(min_val<div_)
        {
            rps.conservativeResize(ps+1,3);
            rps.row(ps)=ref_pt_matrix.row(min_idxR);
            ps++;
        }
    }
//    cout<<ps<<endl;
//    cout<<rps;

    ref_pt_matrix.conservativeResize(rps.rows(), rps.cols());
    ref_pt_matrix=rps;

    return ref_pt_matrix;
}

Eigen::MatrixXf Read_file(int col_, int row_)//col=2500, row=3
{

    ifstream file2("DTLZ1ReferencePoints.dat");
    if (!file2)
    {
        cout << "file can't open" <<endl;
    }
    else
    {
        Eigen::MatrixXf ref_mat(col_,row_);
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
//            cout<<endl;
//            cout<<"vector"<<row_element<<endl;
            ref_mat.row(i)=row_element;
        }


        cout<<"data read finished"<<endl;

        file2.close();
        return ref_mat;
    }
}
