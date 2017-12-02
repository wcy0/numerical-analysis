/*
 *
 * author：X_Y 耀
 *
 * 2017/11/5
 *
 */

#ifndef matrix_H
#define matrix_H

#include <iostream>
#include <cstring>
#include <math.h>
using namespace std;

class matrix
{
public:
    matrix();
    matrix(int ,int ,double*);
    int getRow(){return row;};
    int getColumn(){return column;};
    matrix(const matrix&); //拷贝构造函数  //
    matrix& operator=(const matrix& m);
    matrix  operator*(const matrix&);
    matrix  operator*(double);
    friend ostream& operator<<(ostream&, const matrix&);  //
    inline double TakeElem(int i,int j){return data[i*column+j];}  //取矩阵中某个元素
    inline void SetElem(int i,int j,double b){this->data[i*column+j]=b;}//改变矩阵中某个元素的值
    double norm(char);//求矩阵范数‘1’==1范数，‘f’==F范数，‘i’==无穷范数
    ~matrix();
private:
    int row,column;
    double *data;
};



#endif // matrix_H
