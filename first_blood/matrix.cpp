/*
 *
 * author：X_Y 耀
 *
 * 2017/11/5
 *
 */


#include "matrix.h"

matrix::matrix()
{
    row  = 0;  column = 0;  data = NULL;
}

matrix::matrix(int r, int c, double* p):row(r), column(c)
{
    data = new double[r*c];
    for(int i=0; i<r*c; i++)
        data[i] = p[i];
}

matrix::matrix(const matrix& m)
{
    row    = m.row;  column = m.column;
    data = new double[row * column];
    for(int i=0; i<row*column; i++)
        data[i] = m.data[i];
}

matrix matrix::operator*(const matrix& m)
{
    if(column != m.row)  exit(0);
    matrix temp;
    temp.data = new double[row*m.column];
    temp.row = row;
    temp.column = m.column;
    memset(temp.data, 0, sizeof(double)*row*m.column);
    for(int i=0; i<row; i++)
    {
        for(int j=0;j<m.column;j++)
        {
            for(int k=0;k<column;k++)
                temp.data[i*m.column+j]+=data[i*column+k]*m.data[k*temp.column+j];
        }
    }
    return temp;
}

matrix matrix::operator *(double m)
{
    matrix temp;
    temp.data = new double[row*column];
    temp.row = row;  temp.column = column;
    for(int i=0;i<temp.row*temp.column;i++)
        temp.data[i]=data[i]*m;
    return temp;
}

ostream& operator<<(ostream& o,const matrix& m)//overload the cout
{
    o<<"Print the matrix as follows:"<<endl;
    for(int i=0; i<m.row * m.column; i++)
    {
        if((i+1) % m.column == 0)
            o<<m.data[i]<<endl;
        else
            o<<m.data[i]<<"  ";
    }
    return o;
}

double matrix::norm(char a)
{
    if(a=='1')
    {
        double *t=new double[this->column];
        for(int i=0;i<this->column;i++)
        {
            t[i]=0;
            for(int j=0;j<this->row;j++)
            {
                this->data[j*this->column+i]>0?t[i]+=data[j*this->column+i]:t[i]-=data[j*this->column+i];
            }
        }
        double mx=t[0];
        for(int i=1;i<this->column;i++)
        {
            if(t[i]>mx)mx=t[i];
        }
        delete[] t;
        return mx;
    }
/*    else if(a=='2')
    {

    }*/
    else if(a=='F'||a=='f')
    {
        double m=0;
        for(int i=0;i<this->column*this->row;i++)
            m+=this->data[i]*this->data[i];
        return sqrt(m);
    }
    else if(a!='I'&&a!='i')
    {
        cout<<"use infinite norm as default"<<endl;
    }

    double *t=new double[this->row];
    for(int i=0;i<this->row;i++)
    {
        t[i]=0;
        for(int j=0;j<this->column;j++)
        {
            this->data[i*this->column+j]>0?t[i]+=data[i*this->column+j]:t[i]-=data[i*this->column+j];
        }
    }
    double mx=t[0];
    for(int i=1;i<this->column;i++)
    {
        if(t[i]>mx)mx=t[i];
    }
    delete[] t;
    return mx;

}

matrix& matrix::operator=(const matrix& m)
{
    row = m.row;
    column = m.column;
    data = new double[row*column];
    for(int i=0; i<row*column; i++)
        data[i] = m.data[i];  return *this;
}

matrix::~matrix()
{
    //cout<<"~matrix"<<endl;
    if(data != NULL)
    {
        delete[] data;
        data = NULL;
    }
}
