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
/**pre-define variable*/

DTLZ1_Para_F MOP;
int nVar=MOP.objfun_dim;     // Number of Decision Variables
//    Eigen::MatrixXf Lb(1,MOP.nobjfun), Ub(1,MOP.nobjfun);
//    auto Lb=MOP.searchspaceLB;
//    auto Ub=MOP.searchspaceUB;
float VarMin = 0;//MOP.searchspaceLB(0);         // Decision Variables Lower Bound
float VarMax = 1;//MOP.searchspaceUB(0);         // Decision Variables Upper Bound
int nObj=MOP.nobjfun;
//    float sigma=0.1*(VarMax-VarMin);//0.1 ??
int CrsOver_T=20; // from NSGAIII (Deb) =neighbor size
/**pre-define function*/
class crossover_params_class
{
public:
    int number=20;
    float gamma=0.5; // from NSGAIII (Deb)
    float VarMin_cross;
    float VarMax_cross;
};
class mutation_params_class
{
public:
    float number=20;
    float possibility;
};

crossover_params_class Crossover_params;
mutation_params_class Mutation_params;

struct empty_Subproblem
{
    MatrixXf lambda;
    MatrixXf Neighbors;
};

struct Array_w_idx
{
    MatrixXf Val_;
    MatrixXf idx_;
};

class empty_individual_class
{
public:
    Eigen::MatrixXf Position;
    Eigen::MatrixXf Cost;
    float g;
    bool IsDominated;

    empty_individual_class(int nVar, int nObj)
    {
        Position = Eigen::MatrixXf::Zero(1, nVar);
        Cost = Eigen::MatrixXf::Zero(1, nObj);
        g=0.0f;
        IsDominated=true;
    }
};

/** func. pre-declaration*/
Eigen::Array<empty_Subproblem, 1, Eigen::Dynamic> CreateSubProblems(int nObj, int nPop, int T);

Array_w_idx sort_idx(MatrixXf unsort_arr);

Eigen::MatrixXf pdist2(const Eigen::MatrixXf& A, const Eigen::MatrixXf& B);
float DecomposedCost(std::unique_ptr<empty_individual_class>&, MatrixXf, MatrixXf );
void MOEA_D(DTLZ1_Para_F MOP);
void DetermineDomination(Eigen::Matrix<std::unique_ptr<empty_individual_class>, 1, Eigen::Dynamic>&);
void SortDominatedPop(Eigen::Matrix<std::unique_ptr<empty_individual_class>, 1, Eigen::Dynamic>&,Eigen::Matrix<std::unique_ptr<empty_individual_class>, 1, Eigen::Dynamic>& );
void Crossover_Func(MatrixXf&, const Eigen::MatrixXf, const Eigen::MatrixXf);
void c_test(const Eigen::MatrixXf);

/**Main MOEAD function*/
int main()
{
//    DTLZ1_Para_F MOP;
    //temp initialization for debug-------------------
    MOP.MaxIt=300;
    MOP.nPop=100;
    MOP.nbox=25;
    MOP.objfun_dim=7;
    MOP.Position=MatrixXf::Constant(MOP.nPop, MOP.objfun_dim,0);

    for(int i=0; i<ceil(MOP.nPop); ++i)
    {
        MOP.Position.row(i)=unifrnd(0,1,7);
    }

    //-----------------------------------------------
    MOEA_D(MOP);
    return 0;
}


//MatrixXf EP=MOEA_D(DTLZ1_Para_F MOP) /** TBC*/
void MOEA_D(DTLZ1_Para_F MOP)
{

    // Problem Definition

    MatrixXf(*CostFunction)(MatrixXf)=DTLZ_1;
    Eigen::MatrixXf Lb(1,MOP.nobjfun), Ub(1,MOP.nobjfun);
    Lb=MOP.searchspaceLB;
    Ub=MOP.searchspaceUB;
    float VarMin = Lb(0);         // Decision Variables Lower Bound
    float VarMax = Ub(0);         // Decision Variables Upper Bound
//    float sigma=0.1*(VarMax-VarMin);//0.1 ??

    /** MOEA/D Settings*/

    Crossover_params.VarMin_cross=VarMin;
    Crossover_params.VarMax_cross=VarMax;
    Mutation_params.possibility=1/nVar;

    /** Initialization*/

    Eigen::Array<empty_Subproblem, 1, Eigen::Dynamic> sp=CreateSubProblems(nObj,MOP.nPop,CrsOver_T);// Create Sub-problems


    empty_individual_class empty_individual(nVar,nObj);// Empty Individual

    MatrixXf z_min_pt=MatrixXf::Zero(1,nObj);// Initialize Goal Point

    /** Create Initial Population*/
    //   Eigen::Matrix<empty_individual_class, Eigen::Dynamic, 1> pop(nPop,1);
    Eigen::Matrix<std::unique_ptr<empty_individual_class>, 1, Eigen::Dynamic> pop;
    pop.resize(1, MOP.nPop);

    // initial every element
    for(int i=0; i< MOP.nPop; ++i)
    {
        pop(i)=std::unique_ptr<empty_individual_class>(new empty_individual_class(nVar, nObj));
//            cout<<i<<","<<pop(0,i)->Position(0,0)<<endl;
    }


    for(int i=0; i<MOP.nPop; ++i)
    {
        pop(i)->Position=unifrnd(VarMin,VarMax,nVar);
//        cout<<i<<","<<pop(0,i)->Position<<endl;
    }


    for(int i=0; i<round(MOP.nPop*0.1); ++i)
    {
        pop(i)->Position=MOP.Position.row(i);//had been calculate in Matlab_DTLZ1
    }


    for(int i=0; i<MOP.nPop; ++i)
    {
        pop(i)->Cost=CostFunction(pop(i)->Position);
        z_min_pt=z_min_pt.eval().cwiseMin(pop(i)->Cost);
    }

    for(int i=0; i<MOP.nPop; ++i)
    {
        pop(i)->g=DecomposedCost(pop(i), z_min_pt, sp(i).lambda);
    }

    // Determine Population Domination Status
    DetermineDomination(pop);

    //Initialize Estimated Pareto Front
    Eigen::Matrix<std::unique_ptr<empty_individual_class>, 1, Eigen::Dynamic> Elite_Pop;
    SortDominatedPop(pop,Elite_Pop);

    // Main Loop
    int CrsOverRand, j0, j1;

//    std::unique_ptr<empty_individual_class> y(new empty_individual_class(nVar, nObj));
    empty_individual_class y(nVar, nObj);
    MatrixXf p0(1, nVar),p1(1, nVar);


    Matrix<int, 1, Eigen::Dynamic> vec(1,CrsOver_T);

    std::random_device rd;
    srand(time(NULL));

    for(int t=0; t<CrsOver_T; ++t)
    {
        vec(t)=t;
    }


    for (int it=0; it<MOP.MaxIt; ++it) //MOP.MaxIt
    {
        for(int i=0; i<10; ++i)//MOP.nPop
        {
            // Reproduction (Crossover)
            CrsOverRand=rand()%(Crossover_params.number)+1;  //(rand() % static_cast<int>(Crossover_params.number + 1));//0~20
//            cout<<"_CrsOver:"<<CrsOverRand<<endl;

            for (int nCrossover_iter=0; nCrossover_iter<10; ++nCrossover_iter)//CrsOverRand
            {

                std::shuffle(vec.begin(), vec.end(), std::default_random_engine(rd()));

//                cout<<"iter:"<<nCrossover_iter<<endl;
                j0=sp(i).Neighbors(vec(0))-1;
//                cout<<"j0 : "<<j0 <<endl;

                p0 = pop(j0)->Position;//memory issue??
//                cout<<", p0"<<p0<<endl;

                j1=sp(i).Neighbors(vec(1))-1;
//                cout<<"j1 : "<< j1 <<endl;
                p1 = pop(j1)->Position;
//                cout<<", p1"<<p1<<endl;
                Crossover_Func(y->Position, p0,p1);

//                cout<<"y Position:"<<y.Position<<endl;
                y->Cost=CostFunction(y->Position);
//                cout<<"y Cost:"<<y->Cost<<endl;
                z_min_pt=z_min_pt.eval().cwiseMin(y->Cost);
//              //value of z pt is TBC--------------------------------------------------------!!!
            }


            // Reproduction (Mutaion)

            int sp_N;
            for (int j=0; j<sp(0).Neighbors.cols(); ++j)//neighbor
            {

                sp_N=sp(i).Neighbors(j);
//                cout<<"Neighbor_it  "<<j<<", sp_N  "<<sp_N<<endl;

                y->g=DecomposedCost(y, z_min_pt,sp(sp_N).lambda);//the ith sp,the jth Neighbors
//                cout<<"g:"<<y->g<<endl;
                if(y->g <= pop(sp_N)->g)
                {
//                    pop(sp_N)=std::unique_ptr<empty_individual_class>(new empty_individual_class(y)));
                }
            }

            //ex. pop(i)->g=DecomposedCost(pop(i), z_min_pt, sp(i).lambda);
            /*
                        for (int nMutation_iter=0; nMutation_iter<mutation.number; ++nMutation_iter)
                        {
                            y->Position = Mutate_Func(pop(i).Position, mutation.possibility, sigma);
                            pop(i).Position=findlimits(pop(i).Position,Lb,Ub);
                            y->Cost = feval(CostFunction,y.Position);
                            z_min_pt=min(z_min_pt,y.Cost);
                        }

                        for (int j=0; j<sp(i).Neighbors; ++j)
                        {
                            y.g=DecomposedCost(y,z_min_pt,sp(j).lambda);
                            if(y.g<=pop(j).g)
                            {
                                pop(j)=y;
                            }
                        }
                        delete y;
                        */

        }


        // Determine Population Domination Status
        /*
                pop=DetermineDomination(pop);
                ndpop=pop(~[pop.IsDominated]);
                EP=[EP; ndpop];
                EP=DetermineDomination(EP);
                EP=EP(~[EP.IsDominated]);

                if(Elite_PopP.cols()>nArchive)
                {
                    Extra=numel(EP)-nArchive;
                    ToBeDeleted=randsample(numel(EP),Extra);
                    EP(ToBeDeleted)=[];
                }


                // Display Iteration Information
                cout<<"Iteration"<< it <<": Number of Pareto Solutions ="<<EP.cols()<<endl;
        */
        //for (int i = 0; i < 100; ++i) {
        //    pop(0, i).reset();//delete the storage of pop
        //}

    }


    //Reults
    /*
    EPC=[EP.Cost];
    for (int j=0;j<nObj;++j){
    cout<<"Objective #" num2str(j) ":"<<endl;
        disp<<"      Min = " << min(EPC(j,:))<<endl;
        disp<<"      Max = " << max(EPC(j,:))<<endl;
        disp<<"    Range = " <<(max(EPC(j,:))-min(EPC(j,:)))<<endl;
        disp<<"    St.D. = " <<std(EPC(j,:))<<endl;
        disp<<"    Mean = " <<mean(EPC(j,:))<<endl;
    }
    */

}

/**function define*/
//Create Subproblem
Eigen::Array<empty_Subproblem, 1, Eigen::Dynamic> CreateSubProblems(int nObj, int nPop, int T)
{

    //Create Subproblem with member lambda(dist) and Neighbor(closest pt.)
    Eigen::Array<empty_Subproblem, 1, Eigen::Dynamic> sp(1,nPop);


    MatrixXf lambda_;
    MatrixXf LAMBDA(nPop,nObj);
    MatrixXf Distance(nPop,nPop);
    //create pair of Val & indx
    for(int i=0; i<nPop; ++i)
    {
        lambda_=MatrixXf::Random(1, nObj).array() * 0.5 + 0.5;/**random type TBC*/

        lambda_=lambda_.array()/lambda_.norm();//Norm­

        sp[i].lambda=lambda_;
        LAMBDA.row(i)=lambda_;

    }

    Distance=pdist2(LAMBDA,LAMBDA);
    //Sort by idx & Val, then assign into "sp"
    for(int i=0; i<nPop; ++i)
    {
        sp[i].Neighbors.resize(1,T);
        Array_w_idx Sorted_Array=sort_idx(Distance.row(i));

        sp[i].Neighbors=Sorted_Array.idx_.leftCols(T);
        sp[i].Neighbors.array()-=1;

//        cout<<i<<" _Neighbor: "<<sp[i].Neighbors<<endl;
    }
    return sp;
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
        }
        std::sort(index.begin(),index.end());

        for(int i=0; i<unsort_arr.cols(); ++i)
        {
            Arr_sorted.Val_(j,i)=index[i].first;
            Arr_sorted.idx_(j,i)=index[i].second+1;
        }
    }
    return Arr_sorted;

}

Eigen::MatrixXf pdist2(const Eigen::MatrixXf& A, const Eigen::MatrixXf& B)
{
    if(B.cols()!=A.cols())cout<<"error: Matrix column are different in pdist2!";

    int A_r= A.rows();
    int B_r = B.rows();

    MatrixXf V_dist;

    Eigen::MatrixXf C(A_r, B_r);

    for (int i = 0; i < A_r; i++)
    {
        for (int j = 0; j < B_r; j++)
        {
            V_dist=(A.row(i)-B.row(j)).array().pow(2);
            C(i, j)=(V_dist.sum());
        }
    }
    return C;
}


//DecomposedCost

float DecomposedCost(std::unique_ptr<empty_individual_class>& individual, MatrixXf z, MatrixXf lambda)//pop(i).g=DecomposedCost(pop(i),z,sp(i).lambda);
{

    MatrixXf fx(1,lambda.cols());

    fx=individual->Cost;
    float g=(lambda.array()*abs(fx.array()-z.array())).matrix().maxCoeff();

    return g;

}



void DetermineDomination(Eigen::Matrix<std::unique_ptr<empty_individual_class>, 1, Eigen::Dynamic>& pop)
{
    // Determine Domination
    int nPop=pop.cols();
    MatrixXf Cost_mat(nPop-1,(pop(0)->Cost).cols()),  Cost_mat_origin(nPop,(pop(0)->Cost).cols());
    MatrixXf NPOP_mat;
    Matrix<bool,Eigen::Dynamic,3> compare1(nPop-1,3), compare2(nPop-1,3);
    Matrix<int,Eigen::Dynamic,Eigen::Dynamic>compare1_int, compare2_int,rowsum_c, Ones_set;
    Ones_set.resize(nPop-1,1);
    Ones_set.fill(1);

//---------------------------test case ---------------------------------
//    Cost_mat_origin=MatrixXf::Constant(nPop,3,2);//MatrixXf::Zero(nPop,3)
//    Cost_mat_origin.row(2)=MatrixXf::Constant(1,3,1);//1 1 1
//    Cost_mat_origin.row(3)<<0,1,1;                   //0 1 1
//    Cost_mat_origin.row(4)<<1,1,0;
//---------------------------test case ---------------------------------

    for(int i=0; i<nPop; ++i)
    {
        Cost_mat_origin.row(i)=pop(i)->Cost;
//        cout<<i<<"_"<<Cost_mat_origin.row(i)<<", ";

    }

    Cost_mat=Cost_mat_origin.bottomRows(nPop-1);

    for(int i=0; i<nPop; ++i)
    {
        if(i>0)
        {
            Cost_mat.topRows(i)=Cost_mat_origin.topRows(i);//Cost_mat.rows()=nPop-1=99
            Cost_mat.bottomRows(nPop-i-1)=Cost_mat_origin.bottomRows(nPop-i-1);
        }
        if(i==2)
        {
//                cout<<(i-1)<<"_minus_:"<<Cost_mat.row(i-1)<<", "<<endl;;
//        cout<<i<<" _i value_"<<Cost_mat.row(i)<<", "<<endl;;
//        cout<<(i+1)<<"_plus_"<<Cost_mat.row(i+1)<<", "<<endl;

        }


        NPOP_mat=MatrixXf::Constant(nPop-1,1,1);//reassign
        NPOP_mat=NPOP_mat.eval()*Cost_mat_origin.row(i);

        compare1=NPOP_mat.array() <= Cost_mat.array();//99x3
        compare1_int=compare1.cast<int>();
        rowsum_c= compare1_int.rowwise().sum();
        compare1_int.conservativeResize(nPop, 1);
        compare1_int=rowsum_c;
//        if(i==nPop-1){
//                cout<<"C1:"<<endl;
//        cout<<compare1_int;}

        compare2=NPOP_mat.array() < Cost_mat.array();
        compare2_int=compare2.cast<int>();
        rowsum_c= compare2_int.rowwise().sum();
        compare2_int.conservativeResize(nPop, 1);
        compare2_int=rowsum_c;
//if(i==nPop-1){cout<<"C2:"<<endl;cout<<compare2_int;}


//        Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> diff;//just for test & simplify
//        diff= compare1_int.array()==compare2_int.array();

        if( !(compare1_int.array() <Ones_set.array() ).any() )//why not compare2 only???
        {
            if( !(compare2_int.array() <Ones_set.array() ).any() )//if one pt sum is 0 means this pt is dominated
            {
                pop(i)->IsDominated=false;

            }

        }//compare 99x99 times

//        cout<<i<<"_"<<pop(i)->IsDominated<<endl;
        //To  do list: porduce random case to prove compare2 only is right

    }

}

void SortDominatedPop(Eigen::Matrix<std::unique_ptr<empty_individual_class>, 1, Eigen::Dynamic>& pop, Eigen::Matrix<std::unique_ptr<empty_individual_class>, 1, Eigen::Dynamic>& Elite_Pop)
{
    std::vector<std::unique_ptr<empty_individual_class>> filteredElements;


    for (int i = 0; i < pop.cols(); ++i)
    {
        if (!(pop(i)->IsDominated))
        {
//            auto& person=pop(i);
            filteredElements.push_back(std::unique_ptr<empty_individual_class>(new empty_individual_class(*pop(i))));
        }
    }

    int Element_size = static_cast<int>(filteredElements.size());
    Elite_Pop.resize(1, Element_size);

    for (int i = 0; i < Element_size; ++i)
    {
        Elite_Pop(i) = std::move(filteredElements[i]);
    }

}


void Crossover_Func(Eigen::MatrixXf& y_Position, const Eigen::MatrixXf x1, const Eigen::MatrixXf x2)
{
    srand(time(NULL));

    MatrixXf min_p=MatrixXf::Constant(1,x1.cols(),Crossover_params.VarMin_cross);
    MatrixXf max_p=MatrixXf::Constant(1,x1.cols(),Crossover_params.VarMax_cross);
    // alpha=unifrnd(-gamma,1+gamma,size(x1));
    MatrixXf alpha=MatrixXf::Random(x1.rows(),x1.cols()).array() * 0.5 + 0.5;//1 by n

    y_Position=alpha.array()*x1.array()+(1-alpha.array())*x2.array();
    y_Position=y_Position.eval().cwiseMax(min_p);
    y_Position=y_Position.eval().cwiseMin(max_p);
//    cout<<"Y_after:"<<y_Position<<endl;

}
void Mutate_Func(Eigen::MatrixXf& y_Position, Eigen::MatrixXf pop_Position)//y become pop first, then pick one to mutation.
{
    //sigma=0.1

    int size_j=ceil(Mutation_params.possibility*pop_Position.cols());
    int mu=ceil(Mutation_params.possibility*nVar);

    std::random_device rd;
    std::mt19937 generator( rd() );
    std::uniform_int_distribution<int> distribution(0, nVar-1);

    VectorXi j=VectorXi::Random(ceil(Mutation_params.possibility*pop_Position.cols()));//randsample(n,k)丟k個1~n均勻分布數值
    if(size_j==1){
        j(0)=distribution(generator);
    }
    else{
        for(int i=0;i<mu;++i){
            j(i) = distribution(generator);}
    }

    y_Position=pop_Position; //y become pop first, then pick one to mutation.
//    cout<<"before: "<<y_Position<<endl;

    float randi=sigma*distribution(generator);//randn(n): nxn
    y_Position(0,j)=pop_Position(0,j).array()+randi;//[1x7][7x7]

//    cout<<"after : "<<y_Position<<endl;

}
/*
function b=Dominates(x,y)
{
// Dominates
    if( isfield(x,'Cost') )
        x=x.Cost;

    if isfield(y,'Cost')
        y=y.Cost;

    b=all(x<=y) && any(x<y);
}





function ns=findlimits(ns,Lb,Ub)
{
    // Findlimits
    // Apply the lower bound
    ns_tmp=ns;
    I=ns_tmp<Lb;
    ns_tmp(I)=Lb(I);
    ns_tmp=ns_tmp+I.*rem(abs(Lb-ns),abs(Ub-Lb));

    // Apply the upper bounds
    J=ns_tmp>Ub;
    ns_tmp(J)=Ub(J);
    ns_tmp=ns_tmp-J.*rem(abs(ns-Ub),abs(Ub-Lb));

    // Update this new move
    ns=ns_tmp;
}
*/
