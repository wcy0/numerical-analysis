/*
 *
 * author：X_Y 耀
 *
 * 2017/11/5
 *
 */

#ifndef BANDMATRIX_H
#define BANDMATRIX_H

#include <iostream>
#include <cstring>
#include <math.h>
#include <matrix.h>

class bandmatrix:private matrix
{
public:
    bandmatrix();
    bandmatrix(int r,int c,int w,int u,double* p);
    virtual ~bandmatrix();
    bandmatrix(const bandmatrix&);
    bandmatrix& operator=(const bandmatrix&);
    matrix operator *(matrix &);
    friend ostream& operator<<(ostream&, bandmatrix&);
    double TakeElem(int i, int j);
    void SetElem(int,int,int);
    bandmatrix DoolittleLU();//Doolittle分解
    friend matrix DoolittleDecomposition(bandmatrix ,matrix); //带状doolittle分解法
    double max_eig(double, double);//按模最大特征值
    double min_eig(double, double);//按模最小特征值
    matrix band2mat();
    double detmatrix();

    int getRow(){return row;}
    int getColumn(){return column;}
    int getWidth(){return width;}
    int getWidthUp(){return width_up;}
    int getWidthBe(){return width_be;}
    void setRow(int i){row=i;}
    void setColumn(int i){column=i;}
    void setWidth(int i){width=i;}
    void setWidthUp(int i){width_up=i;}
    void setwidthBe(int i){width_be=i;}
    double* getData(){return  data;}



private:
    int row;
    int column;
    int width;
    int width_up;
    int width_be;
    double *data;
};

#endif // BANDMATRIX_H
