#include <iostream>
#include <fstream>
#include <random>
#include <algorithm> //for std::shuffle as randsample in matlab
#include <cstdlib> //rand()
#include <ctime> //time()
#include <Eigen/Dense>
#include <type_traits>//for isfield>>has_cost
#include <memory>//std::unique_ptr

#include "DTLZ1_Para.h"
#include "DTLZ1_function.h"
#include "unifrnd.h"
//#include "factorial.h"
#include "randsample.h"


using namespace std;
using namespace Eigen;
/**pre-define variable*/

DTLZ1_Para_F MOP;
int nVar=MOP.objfun_dim;          // Number of Decision Variables
//    Eigen::MatrixXf Lb(1,MOP.nobjfun), Ub(1,MOP.nobjfun);
//    auto Lb=MOP.searchspaceLB;
//    auto Ub=MOP.searchspaceUB;
float VarMin = 0;//MOP.searchspaceLB(0);         // Decision Variables Lower Bound
float VarMax = 1;//MOP.searchspaceUB(0);         // Decision Variables Upper Bound
int nObj=MOP.nobjfun;
float sigma=0.1*(VarMax-VarMin);//0.1 ??
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
float DecomposedCost(empty_individual_class, MatrixXf, MatrixXf );
void MOEA_D(DTLZ1_Para_F MOP);
void DetermineDomination(std::vector<empty_individual_class>&);
void SortDominatedPop(std::vector<empty_individual_class>&, std::vector<empty_individual_class>& );
void Crossover_Func(Eigen::MatrixXf&, const Eigen::MatrixXf, const Eigen::MatrixXf);
void c_test(const Eigen::MatrixXf);
void Mutate_Func(Eigen::MatrixXf& y_Position, Eigen::MatrixXf x);
void findlimits(Eigen::MatrixXf &pop_Position, Eigen::MatrixXf Lb, Eigen::MatrixXf Ub);

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

//    int nVar=MOP.objfun_dim;          // Number of Decision Variables
    Eigen::MatrixXf Lb(1,MOP.nobjfun), Ub(1,MOP.nobjfun);
    Lb=MOP.searchspaceLB;
    Ub=MOP.searchspaceUB;
    float VarMin = Lb(0);         // Decision Variables Lower Bound
    float VarMax = Ub(0);         // Decision Variables Upper Bound
//    int nObj=MOP.nobjfun;
//    float sigma=0.1*(VarMax-VarMin);//0.1 ??

    // MOEA/D Settings
    MOP.nRep=round(1.5*(pow(MOP.nbox, MOP.nobjfun-1))/(factorial(MOP.nobjfun-1)));//=469


    Crossover_params.VarMin_cross=VarMin;
    Crossover_params.VarMax_cross=VarMax;
    Mutation_params.possibility=1.0f/nVar;//1 must be float so that divide work properly.

    /** Initialization*/
    // Create Sub-problems
    Eigen::Array<empty_Subproblem, 1, Eigen::Dynamic> sp=CreateSubProblems(nObj,MOP.nPop,CrsOver_T);


    // Empty Individual
    empty_individual_class empty_individual(nVar,nObj);
    /** Initialize Goal Point*/
    MatrixXf z_min_pt=MatrixXf::Zero(1,nObj);//min. pt.

    // Create Initial Population
//   Eigen::Matrix<empty_individual_class, Eigen::Dynamic, 1> pop(nPop,1);
//    Eigen::Matrix<empty_individual_class, 1, Eigen::Dynamic> pop;
//pop.conservativeResize(1,MOP.nPop)
    std::vector<empty_individual_class> pop, non_dominate_pop;
    pop.resize(MOP.nPop, empty_individual_class(nVar, nObj));


    // initial every element
//    for(int i=0; i< MOP.nPop; ++i)
//    {
//        pop(i)=empty_individual_class(nVar, nObj);
//    }


    for(int i=0; i<MOP.nPop; ++i)
    {
        pop[i].Position=unifrnd(VarMin,VarMax,nVar);
//        cout<<i<<","<<pop[i].Position<<endl;
    }


    for(int i=0; i<round(MOP.nPop*0.1); ++i)
    {
        pop[i].Position=MOP.Position.row(i);//had been calculate in Matlab_DTLZ1
    }


    for(int i=0; i<MOP.nPop; ++i)
    {
        pop[i].Cost=CostFunction(pop[i].Position);
//        cout<<pop[i].Cost<<endl;
        z_min_pt=z_min_pt.eval().cwiseMin(pop[i].Cost);
    }
//        cout<<z_min_pt;


    for(int i=0; i<MOP.nPop; ++i)
    {
        pop[i].g=DecomposedCost(pop[i], z_min_pt, sp(i).lambda);

    }

    // Determine Population Domination Status
    DetermineDomination(pop);


    //Initialize Estimated Pareto Front
    std::vector<empty_individual_class> Elite_Pop;

    SortDominatedPop(pop,Elite_Pop);

    // Main Loop
    int CrsOverRand, j0, j1;

    empty_individual_class y(nVar, nObj);
    MatrixXf p0(1, nVar),p1(1, nVar);


    Matrix<int, 1, Eigen::Dynamic> vec(1,CrsOver_T);

    std::random_device rd;
    srand(time(NULL));

    for(int t=0; t<CrsOver_T; ++t)
    {
        vec(t)=t;
    }


    for (int it=0; it<20; ++it) //MOP.MaxIt
    {
        for(int i=0; i<MOP.nPop; ++i)//MOP.nPop
        {
            // Reproduction (Crossover)
//            cout<< i <<endl;
            CrsOverRand=rand()%(Crossover_params.number)+1;  //(rand() % static_cast<int>(Crossover_params.number + 1));//0~20
//            cout<<"_CrsOver:"<<CrsOverRand<<endl;
            y.Position=MatrixXf::Constant(1,nVar,0);
            for (int nCrossover_iter=0; nCrossover_iter<CrsOverRand; ++nCrossover_iter)//CrsOverRand
            {

                std::shuffle(vec.begin(), vec.end(), std::default_random_engine(rd()));

//                cout<<sp(i).Neighbors<<endl;

//                cout<<"iter:"<<nCrossover_iter<<endl;
                j0=sp(i).Neighbors(vec(0));
//                cout<<"j0 : "<<j0 <<endl;

                p0 = pop[j0].Position;//memory issue??
//                cout<<", p0:"<<p0<<endl;

                j1=sp(i).Neighbors(vec(1));
//                cout<<"j1 : "<< j1 <<endl;
                p1 = pop[j1].Position;
//                cout<<", p1:"<<p1<<endl;
                Crossover_Func(y.Position, p0,p1);

//                cout<<"y Position:"<<y.Position.rows()<<y.Position.cols()<<"_:"<<y.Position<<endl;
                y.Cost=CostFunction(y.Position);
//                cout<<"y Cost:"<<y.Cost<<endl;
                z_min_pt=z_min_pt.eval().cwiseMin(y.Cost);
//              //value of z pt is TBC--------------------------------------------------------!!!
            }


            // Reproduction (Mutaion)

            int sp_N;
            for (int j=0; j<sp(0).Neighbors.cols(); ++j)//neighbor
            {

                sp_N=sp(i).Neighbors(j);
                //                cout<<"Neighbor_it  "<<j<<", sp_N  "<<sp_N<<endl;

                y.g=DecomposedCost(y, z_min_pt,sp(sp_N).lambda);//the ith sp,the jth Neighbors
                //                cout<<"g:"<<y->g<<endl;
                if(y.g <= pop[sp_N].g)
                {
//                                pop(sp_N)=std::unique_ptr<empty_individual_class>(new empty_individual_class(y)));
                    pop[sp_N]=y;
                }
            }


            //ex. pop(i)->g=DecomposedCost(pop(i), z_min_pt, sp(i).lambda);

            for (int nMutation_iter=0; nMutation_iter<Mutation_params.number; ++nMutation_iter)//Mutation_params.number
            {

                Mutate_Func(y.Position, pop[i].Position);
                findlimits(pop[i].Position,Lb,Ub);
                y.Cost=CostFunction(y.Position);
                z_min_pt=z_min_pt.eval().cwiseMin(y.Cost);

            }

            //cout<<"sp(i).Neighbors;"<<sp(i).Neighbors;
            for (int j=0; j<sp(0).Neighbors.cols(); ++j)//neighbor
            {
                sp_N=sp(i).Neighbors(j);
                y.g=DecomposedCost(y,z_min_pt,sp(sp_N).lambda);
                if(y.g <= pop[sp_N].g)
                {
                    pop[sp_N]=y;
//                    cout<< "pop[sp_N]" << pop[sp_N]<<endl;
                }
            }

        }


        // Determine Population Domination Status

/** -------------------------------------------------cant define sixe of EP as 100,  it should be dynamic= =*/
        DetermineDomination(pop);
        for(unsigned int k=0; k<pop.size(); ++k)
        {
//            cout<<k<<" , IsDominated"<<pop[k].IsDominated<<endl;
            if(!pop[k].IsDominated)
            {
                non_dominate_pop.push_back(pop[k]);
                Elite_Pop.push_back(pop[k]);
            }
        }

        DetermineDomination(Elite_Pop);

        Elite_Pop.erase( remove_if( Elite_Pop.begin(), Elite_Pop.end(), [](const empty_individual_class& sub_pop)
        {
            return sub_pop.IsDominated;
        }), Elite_Pop.end());




        if(Elite_Pop.size() > MOP.nRep)//MOP.nRep
        {
            cout<<"Over Elite Size: "<< Elite_Pop.size()<<endl;
            int Extra=Elite_Pop.size() - MOP.nRep;//MOP.nRep
            cout<<"Elite Size: "<< Elite_Pop.size() <<", nRep: "<<MOP.nRep <<", Extra: "<<Extra<<endl;
            VectorXi ToBeDeleted;
            ToBeDeleted=randsample(Elite_Pop.size(), Extra);
            cout<<"ToBeDeleted"<< ToBeDeleted<<endl;
            Elite_Pop.erase(std::remove_if(Elite_Pop.begin(), Elite_Pop.end(), [&](const empty_individual_class& element) //remove_if(begin,end, lambda formula)
            {
                int index_E=&element - &Elite_Pop[0];//calculate the index
                return std::find(ToBeDeleted.begin(),ToBeDeleted.end(), index_E) !=ToBeDeleted.end();//confirm if indexE in ToBeDeleted
            }), Elite_Pop.end());

            cout<<"Final Elite Size: "<< Elite_Pop.size()<<endl;

        }

        // Display Iteration Information
        cout<<"Iteration"<< it <<": Number of Elite Pareto Solutions ="<<Elite_Pop.size()<<endl;

        //for (int i = 0; i < 100; ++i) {
        //    pop(0, i).reset();//delete the storage of pop
        //}

    }


    //Reults

    MatrixXf EP_Cost(Elite_Pop.size(),nObj);
    for(unsigned int ii=0; ii<Elite_Pop.size(); ++ii)
    {
        EP_Cost.row(ii)=Elite_Pop[ii].Cost;
    }

//    for (int j=0;j<nObj;++j){
//    cout<<"Objective #" num2str(j) ":"<<endl;
//        cout<<"      Min = " << min(EPC(j,:))<<endl;
//        cout<<"      Max = " << max(EPC(j,:))<<endl;
//        cout<<"    Range = " <<(max(EPC(j,:))-min(EP_Cost.row(j)))<<endl;
//        cout<<"    St.D. = " <<std(EPC(j,:))<<endl;
//        cout<<"    Mean = " <<mean(EPC(j,:))<<endl;
//    }


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
    for(int i=0; i<nPop; ++i)//nPop
    {
        lambda_=MatrixXf::Random(1, nObj).array() * 0.5 + 0.5;/**random type TBC*/
//        cout<<lambda_<<endl;

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

float DecomposedCost(empty_individual_class individual, MatrixXf z, MatrixXf lambda)//pop(i).g=DecomposedCost(pop(i),z,sp(i).lambda);
{

    MatrixXf fx(1,lambda.cols());

    fx=individual.Cost;
    float g=(lambda.array()*abs(fx.array()-z.array())).matrix().maxCoeff();

    return g;

}



void DetermineDomination(std::vector<empty_individual_class>& pop)
{
    // Determine Domination        Cost_mat.rows()=nPop-1
    MatrixXf Cost_mat(pop.size()-1,(pop[0].Cost).cols()),  Cost_mat_origin(pop.size(),(pop[0].Cost).cols());
    MatrixXf NPOP_mat;
    Matrix<bool,Eigen::Dynamic,3> compare1(pop.size()-1,3), compare2(pop.size()-1,3);
    Matrix<int,Eigen::Dynamic,Eigen::Dynamic>compare1_int, compare2_int,rowsum_c, Ones_set;
    Ones_set.resize(pop.size()-1,1);
    Ones_set.fill(1);

//---------------------------test case ---------------------------------
//    Cost_mat_origin=MatrixXf::Constant(nPop,3,2);//MatrixXf::Zero(nPop,3)
//    Cost_mat_origin.row(2)=MatrixXf::Constant(1,3,1);//1 1 1
//    Cost_mat_origin.row(3)<<0,1,1;                   //0 1 1
//    Cost_mat_origin.row(4)<<1,1,0;
//---------------------------test case ---------------------------------


    for(unsigned int i=0; i<pop.size(); ++i)
    {
        Cost_mat_origin.row(i)=pop[i].Cost;
//        cout<<i<<"_"<<Cost_mat_origin.row(i)<<", ";

    }

    Cost_mat=Cost_mat_origin.bottomRows(pop.size()-1);

    for(unsigned int i=0; i<pop.size(); ++i)
    {
        if(i>0)
        {
            Cost_mat.topRows(i)=Cost_mat_origin.topRows(i);//Cost_mat.rows()=nPop-1=99 ,Cost_mat saves cost except ith cost
            Cost_mat.bottomRows(pop.size()-i-1)=Cost_mat_origin.bottomRows(pop.size()-i-1);
        }
        if(i==2)
        {
//                cout<<(i-1)<<"_minus_:"<<Cost_mat.row(i-1)<<", "<<endl;;
//        cout<<i<<" _i value_"<<Cost_mat.row(i)<<", "<<endl;;
//        cout<<(i+1)<<"_plus_"<<Cost_mat.row(i+1)<<", "<<endl;

        }


        NPOP_mat=MatrixXf::Constant(pop.size()-1,1,1);//reassign
        NPOP_mat=NPOP_mat.eval()*Cost_mat_origin.row(i);

        compare1=NPOP_mat.array() <= Cost_mat.array();//99x3
        compare1_int=compare1.cast<int>();
        rowsum_c= compare1_int.rowwise().sum();
        compare1_int.conservativeResize(pop.size(), 1);
        compare1_int=rowsum_c;
//        if(i==nPop-1){
//                cout<<"C1:"<<endl;
//        cout<<compare1_int;}

        compare2=NPOP_mat.array() < Cost_mat.array();
        compare2_int=compare2.cast<int>();
        rowsum_c= compare2_int.rowwise().sum();
        compare2_int.conservativeResize(pop.size(), 1);
        compare2_int=rowsum_c;
//if(i==nPop-1){cout<<"C2:"<<endl;cout<<compare2_int;}


//        Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> diff;//just for test & simplify
//        diff= compare1_int.array()==compare2_int.array();

        if( !(compare1_int.array() <Ones_set.array() ).any() )//why not compare2 only???
        {
            if( !(compare2_int.array() <Ones_set.array() ).any() )//if one pt sum is 0 means this pt is dominated
            {
                pop[i].IsDominated=false;

            }

        }//compare 99x99 times

//        cout<<i<<"_"<<pop(i)->IsDominated<<endl;
        //To  do list: porduce random case to prove compare2 only is right

    }

}

void SortDominatedPop(std::vector<empty_individual_class>& pop, std::vector<empty_individual_class>& Elite_Pop)
{
    std::vector<empty_individual_class> filteredElements;


    for (int i = 0; i < MOP.nPop; ++i)
    {
        if (!(pop[i].IsDominated))
        {
//            auto& person=pop(i);
            filteredElements.push_back(pop[i]);
        }
    }
//
//    int Element_size = static_cast<int>(filteredElements.size());
//    Elite_Pop.resize(1, Element_size);
//
//    for (int i = 0; i < Element_size; ++i)
//    {
//        Elite_Pop(i) = std::move(filteredElements[i]);
//    }

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

    int size_j=ceil(Mutation_params.possibility*pop_Position.cols());//=1, y position only mutate 1 item
    int mu=ceil(Mutation_params.possibility*nVar);//=1

    std::random_device rd;
    std::mt19937 generator( rd() );
    std::uniform_int_distribution<int> distribution(0, nVar-1);

    VectorXi j=VectorXi::Random(size_j);//randsample(n,k)丟k個1~n均勻分布數值  randsample(nVar,nMu)
    if(size_j==1)
    {
        j(0)=distribution(generator);
    }
    else
    {
        for(int i=0; i<mu; ++i)
        {
            j(i) = distribution(generator);
        }
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



*/

void findlimits(Eigen::MatrixXf &pop_Position, Eigen::MatrixXf Lb, Eigen::MatrixXf Ub)
{
    // Findlimits
    // Apply the lower bound


    Eigen::MatrixXf pop_temp(1,nVar), L_float(1,nVar), U_float(1,nVar), B_remainder(1,nVar), B_Quotient(1,nVar);
    Eigen::MatrixXf A(1,nVar), B(1,nVar);

    pop_temp=pop_Position;
    Eigen::Matrix<bool, 1, 7> L_bool, U_bool;//=ns_tmp<Lb;



    L_bool=pop_temp.array() < Lb.array();
    for(int i=0; i<nVar; i++)
    {
        if(pop_temp(i)<Lb(i)) pop_temp(i) = Lb(i);
    }

    B_Quotient=abs(Lb.array()-pop_Position.array()) .cwiseQuotient(abs(Ub.array()-Lb.array()));//Quotient=A /B ;
    B_remainder=abs(Lb.array()-pop_Position.array()) - (B_Quotient.array() * (Ub.array()-Lb.array()));

    L_float=L_bool.cast<float>();

    pop_temp=pop_temp.array() + L_float.array() * B_remainder.array();
//    cout<<"after LB: "<<pop_temp<<endl;


    // Apply the upper bounds


    U_bool=pop_temp.array() < Ub.array();
    for(int i=0; i<nVar; i++)
    {
        if(pop_temp(i)>Ub(i)) pop_temp(i) = Ub(i);

    }


    B_Quotient=abs(Ub.array()-pop_Position.array()) .cwiseQuotient(abs(Ub.array()-Lb.array()));//Quotient=A /B ;
    B_remainder=abs(Ub.array()-pop_Position.array()) - (B_Quotient.array() * (Ub.array()-Lb.array()));

    U_float=U_bool.cast<float>();
    pop_temp=pop_temp.array() + L_float.array() * B_remainder.array();
//    cout<<"after UB: "<<pop_temp<<endl;


    // Update this new move

    pop_Position=pop_temp;
}

