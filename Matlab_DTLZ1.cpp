#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "DTLZ1_Para.h"
#include "DTLZ1_function.h"
//#include "test_head.h"
#include "MOEAD_function.h"
#include"Generate_Ref_Pts.h"


#include<numeric>
using namespace std;

// Test Problem 4: DTLZ1



Eigen::MatrixXf Read_file(int col_, int row_);//col=2500, row=3
//
Eigen::MatrixXf Slope(Eigen::MatrixXf Zr_matrix, Eigen::MatrixXf ref_pt_matrix, int psdiv);


int main()
{

    DTLZ1_Para_F dtlz1_para;
    int objfun_dim = dtlz1_para.objfun_dim;//7
    int nobjfun = dtlz1_para.nobjfun;//3
    MatrixXf searchspaceUB = dtlz1_para.searchspaceUB;
    MatrixXf searchspaceLB = dtlz1_para.searchspaceLB;
    string testfun_name = dtlz1_para.dtlz_name;


    int ps=1, psdiv=15;
    Eigen::MatrixXf  Zr_matrix, reference_ps_, MOEAD_IGD, MOEAD_hv;
    Eigen::MatrixXf  ref_pt_matrix(2500,3);//memory Dynamic may be need

    //initial data

    Zr_matrix = GenerateReferencePoints(dtlz1_para.nobjfun, psdiv);//3,15 Zr(136,3)
    int col_Zr=sizeof(Zr_matrix)/sizeof(Zr_matrix[0]);//24/4

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

//    MatrixXf temp_rand;
    dtlz1_para.Position=MatrixXf::Constant(dtlz1_para.nPop, dtlz1_para.objfun_dim, 0);

    for (int i=0; i<10; ++i)//dtlz1_para.MaxIt
    {

    cout<<"para.MaxItinEX"<<i<<endl;
        for(int j=0; j<10; ++j)//dtlz1_para.nPop
        {

//            cout<<"searchspaceLB: "<<dtlz1_para.searchspaceLB<<endl;
//            cout<<"searchspaceLB: "<<dtlz1_para.searchspaceUB<<endl;
//            cout<<"Position.rows"<<dtlz1_para.Position.rows()<<endl;
//            cout<<"Position.cols"<<dtlz1_para.Position.cols()<<endl;
            dtlz1_para.Position.row(i)=dtlz1_para.searchspaceLB.array() + (dtlz1_para.searchspaceUB.array() - dtlz1_para.searchspaceLB.array()) *unifrnd(0,1,dtlz1_para.objfun_dim).array();

            cout<<dtlz1_para.Position.row(i)<<endl;
        }


        /**Algorithm execution*/



            std::vector<empty_individual_class> rep;
            rep=MOEAD_function(dtlz1_para);

            cout<< rep.size()<<endl;
            Eigen::MatrixXf obtained_ps_MOEAD(rep.size(),nobjfun);
            Eigen::MatrixXf MOEAD_position(rep.size(), objfun_dim);
            for(int i_rep=0;i_rep<rep.size();++i_rep){
                obtained_ps_MOEAD.row(i_rep)=rep[i_rep].Cost;
                MOEAD_position.row(i_rep)=rep[i_rep].Position;
            }
//
            cout<<"cost: "<<obtained_ps_MOEAD.rows()<<","<<obtained_ps_MOEAD.cols()<<endl;
            cout<<"pos: "<<MOEAD_position.rows()<<","<<MOEAD_position.cols()<<endl;


        //        MOEAD_IGD(di)=IGD_calculation(obtained_ps_MOEAD,reference_ps);
        //        NadirPtMOEA=max(obtained_ps_MOEAD,[],1);
        //        MOEAD_hv(di)=hypervolume(obtained_ps_MOEAD,NadirPtMOEA,1000);

    }

    return 0;

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



// IGD

//    indicator=[MOEAD_IGD MOEAD_hv];

// saves files






