
#include <iostream>
#include <fstream>
#include <conio.h>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;


void Read_file(float* p, int col, int row)
{

    ifstream file2("DTLZ1ReferencePoints.dat");
    if (!file2)
    {
        cout << "file can't open" <<endl;
    }
    else
    {
        for(int i=0; i<row; ++i)
        {
            for(int j=0; j<col; ++j)
            {
                file2>>*(p+(i*col)+j);
            }
//            printf("%d : %f, %f, %f  \n",i,*(p+(i*col)), *(p+(i*col)+1), *(p+(i*col)+2) );
        }

//                    cout<<ref_matrix[i][0]<<","<<ref_matrix[i][1]<<","<<ref_matrix[i][2]<<endl;


        cout<<"data read finished"<<endl;

        file2.close();
    }
}
int main()
{
    float ref_matrix[2500][3];
    float *p=&ref_matrix[0][0];
//    float* p[2500];
//    for(int i=0; i<2500; ++i)
//    {
//        p[i]=&ref_matrix[i];
//    }

    Read_file(p,3,3);
    cout<<"50,"<<ref_matrix[50][0]<<","<<ref_matrix[50][1]<<","<<ref_matrix[50][2]<<endl;
    const float a=ref_matrix[501][0];
    const float b=ref_matrix[501][1];
    const float c=ref_matrix[501][2];
    printf("a:%f, b:%f, c:%f\n",a,b,c);

    cout<<*(p+((501*3)+2))<<endl;//p+((i*col)+j)

    return 0;
//    getchar();

}





