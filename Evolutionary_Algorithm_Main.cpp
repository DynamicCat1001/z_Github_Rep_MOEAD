#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "DTLZ1_Para.h"
#include "DTLZ1_function.h"
#include "MOEAD_function.h"
#include"Generate_Ref_Pts.h"
//#include"Algorithm_Plot.h"
#include"Slope.h"
#include"Read_Ref_Pt.h"

#include<numeric>
using namespace std;

// Test Problem 4: DTLZ1



//Eigen::MatrixXf Read_file(int col_, int row_);//col=2500, row=3



int main()
{

    DTLZ1_Para_F dtlz1_para;
    int objfun_dim = dtlz1_para.objfun_dim;//7
    int nobjfun = dtlz1_para.nobjfun;//3
    MatrixXf searchspaceUB = dtlz1_para.searchspaceUB;
    MatrixXf searchspaceLB = dtlz1_para.searchspaceLB;
    string testfun_name = dtlz1_para.dtlz_name;


    int psdiv=15;//ps=1
    Eigen::MatrixXf  Zr_matrix, reference_ps_, MOEAD_IGD, MOEAD_hv;
    Eigen::MatrixXf  ref_pt_matrix(2500,3);//memory Dynamic may be need

    //initial data

    Zr_matrix = GenerateReferencePoints(dtlz1_para.nobjfun, psdiv);//3,15 Zr(136,3)
//    int col_Zr=sizeof(Zr_matrix)/sizeof(Zr_matrix[0]);//24/4

    //Read Reference Point
    ref_pt_matrix=Read_file(2500,3);

    //cross pt of shortest distance
    reference_ps_=Slope(Zr_matrix,ref_pt_matrix,psdiv);

    /**Algorithm Parameters*/
    dtlz1_para.MaxIt=300;
    dtlz1_para.nPop=100;//Population Size
    dtlz1_para.nbox=25;
    dtlz1_para.nRep=round(1.5*(pow(dtlz1_para.nbox, dtlz1_para.nobjfun-1))/(factorial(dtlz1_para.nobjfun-1)));// Repository Size factorial:5!=120



//    MOEAD_IGD.setZero(1,data_iter);
//    MOEAD_hv=MOMFAPref_IGD;

    dtlz1_para.Position=MatrixXf::Constant(dtlz1_para.nPop, dtlz1_para.objfun_dim, 0);

    for (int i=0; i<20; ++i)//dtlz1_para.MaxIt
    {


        for(int j=0; j<dtlz1_para.nPop; ++j)//dtlz1_para.nPop
        {
            dtlz1_para.Position.row(i)=dtlz1_para.searchspaceLB.array() + (dtlz1_para.searchspaceUB.array() - dtlz1_para.searchspaceLB.array()) *unifrnd(0,1,dtlz1_para.objfun_dim).array();

//            cout<<j<<"_"<<dtlz1_para.Position.row(i)<<endl;
        }


        /**Algorithm execution*/



            std::vector<empty_individual_class> rep;
            rep=MOEAD_function(dtlz1_para);

            cout<< rep.size()<<endl;

            Eigen::Matrix<float,Eigen::Dynamic,nobjfun> obtained_cost_MOEAD;
            Eigen::Matrix<float,Eigen::Dynamic,objfun_dim> obtained_position_MOEAD;
            for(int i_rep=0;i_rep<rep.size();++i_rep){
                obtained_cost_MOEAD.row(i)=rep[i].Cost;
                obtained_position_MOEAD.row(i)=rep[i].Position;
            }
//
//            cout<<"cost: "<<obtained_cost_MOEAD.rows()<<","<<obtained_cost_MOEAD.cols()<<endl;
//            cout<<"pos: "<<obtained_position_MOEAD.rows()<<","<<obtained_position_MOEAD.cols()<<endl;


        //        MOEAD_IGD(di)=IGD_calculation(obtained_cost_MOEAD,reference_ps);
        //        NadirPtMOEA=max(obtained_cost_MOEAD,[],1);
        //        MOEAD_hv(di)=hypervolume(obtained_cost_MOEAD,NadirPtMOEA,1000);

    }

    return 0;

}


//Eigen::MatrixXf Read_file(int col_, int row_)//col=2500, row=3
//{
//
//    ifstream file2("DTLZ1ReferencePoints.dat");
//    if (!file2)
//    {
//        cout << "file can't open" <<endl;
//    }
//    else
//    {
//        Eigen::MatrixXf p(col_,row_);
//        Eigen::RowVector3f row_element;
//        float A;
//        for(int i=0; i<col_; ++i)//2500
//        {
//            for(int j=0; j<row_; ++j)//3 dim.
//            {
//                A=0;
//                file2>>A;
//                row_element[j]=A;
//            }
//            p.row(i)=row_element;
//        }
//
//
//        cout<<"data read finished"<<endl;
//
//        file2.close();
//        return p;
//    }
//}


// IGD

//    indicator=[MOEAD_IGD MOEAD_hv];

// saves files






