/*
 *
 * author：X_Y 耀
 *
 * 2017/11/5
 *
 */

#include <QCoreApplication>
#include <iostream>
#include <iomanip>
#include <matrix.h>
#include <bandmatrix.h>
#include <fstream>
#include <assert.h>
#include <math.h>

using namespace std;


int main()
{
    double data[501*501]={0},b,c;
    int row,column,width,width_up;
    double lamda_1,lamda_501,lamda_max,lamda_s,mu[40],lamda_ik[40],cond_A2,det_A;
    bandmatrix A;
    matrix B;

    c=-0.064;
    b=0.16;

    ifstream fin("mat_A.txt");
    if(!fin)
    {
        ofstream fout("mat_A.txt");
        for(int i=2;i<501;i++)
            data[i]=c;

        for(int i=502;i<501*2;i++)
            data[i]=b;

        for(int i=1;i<=501;i++)
            data[i+2*501-1]=(1.64-0.024*i)*sin(0.2*i)-0.64*exp(0.1/i);

        for(int i=501*3;i<501*4-1;i++)
            data[i]=b;

        for(int i=501*4;i<501*5-2;i++)
            data[i]=c;

        A=bandmatrix(501,501,5,2,data);
        fout<<A;
    }
    else
    {
        assert(fin);
        fin>>row>>column>>width>>width_up;
        for(int i=0;i<501*5;i++)
            fin>>data[i];
        A=bandmatrix(row,column,width,width_up,data);
    }

    lamda_max=lamda_1=A.max_eig(1e-12,0);
    lamda_501=A.max_eig(1e-12,lamda_1)+lamda_1;
    if(lamda_1>lamda_501)
    {
        lamda_s=lamda_1;
        lamda_1=lamda_501;
        lamda_501=lamda_s;
    }
    lamda_s=A.min_eig(1e-12,0);

    for(int i=1;i<40;i++)
    {
        mu[i]=lamda_1+1.0*i*(lamda_501-lamda_1)/40;
        lamda_ik[i]=A.min_eig(1e-12,mu[i])+mu[i];
    }

    cond_A2=fabs(lamda_max/lamda_s);
    det_A=A.detmatrix();
    cout<< setiosflags(ios::uppercase)<<setiosflags(ios::scientific)<<setprecision(11);
    cout<<"lamda1: "<<lamda_1<<endl;
    cout<<"lamda501: "<<lamda_501<<endl;
    cout<<"lamda_s: "<<lamda_s<<endl;
    for(int i=1;i<40;i++)
    {
            cout<<"lamda_i"<<i<<": "<<lamda_ik[i];
            if(i%3==0)
            {
                cout<<endl;
            }
            else
            {
                cout<<'\t';
            }
    }
    cout<<"cond_A2: "<<cond_A2<<endl;
    cout<<"det A: "<<det_A<<endl;
}
