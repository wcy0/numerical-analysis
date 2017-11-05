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
    matrix  operator+(const matrix&);
    matrix  operator-(const matrix&);
    matrix  operator*(const matrix&);
    matrix  operator*(double);
    matrix  operator/(matrix&);
    matrix& operator=(const matrix&);
    matrix& operator+=(const matrix&);
    matrix& operator-=(const matrix&);
    friend ostream& operator<<(ostream&, const matrix&);  //
    inline double TakeElem(int i,int j){return data[i*column+j];}  //取矩阵中某个元素
    inline void SetElem(int i,int j,double b){this->data[i*column+j]=b;}//改变矩阵中某个元素的值
    matrix operator[](int); //取某行，对操作符[]重载
    matrix operator()(int); //取某列，对操作符()重载
    matrix TakeSubBlock(int, int, int, int); //取子矩阵  //
    matrix Rotate(); //转秩
    matrix Inverse(); //逆矩阵
    friend double detmatrix(const matrix); //求方阵的行列式值
    friend matrix GaussEliminationSeq(const matrix,const matrix); //高斯消去法
    friend matrix GaussEliminationColPrin(const matrix,const matrix); //列主元高斯消去法
    friend matrix DoolittleDecomposition(const matrix ,const matrix); //doolittle分解法
    friend matrix DltDpsCP(const matrix,const matrix); //选主元doolittle分解法
    double norm(char);//求矩阵范数‘1’==1范数，‘f’==F范数，‘i’==无穷范数
    double max_eig(double, double);//按模最大特征值
    double min_eig(double, double);//按模最小特征值
    ~matrix();
private:
    int row,column;
    double *data;
};



#endif // matrix_H
