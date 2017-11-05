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

matrix matrix::operator+(const matrix& m)
{
    if(row != m.row || column != m.column)
        exit(0);
    matrix temp;
    temp.data = new double[row*column];
    temp.row = row;  temp.column = column;
    for(int i=0; i<row*column; i++)
        temp.data[i] = data[i] + m.data[i];
    return temp;
}

matrix matrix::operator-(const matrix& m)
{

    if(row != m.row || column != m.column)
        exit(0);  matrix temp;
    temp.data = new double[row*column];
    temp.row = row;  temp.column = column;
    for(int i=0; i<row*column; i++)
        temp.data[i] = data[i] - m.data[i];
    return temp;
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

matrix matrix::operator/(matrix& m)
{
    return (*this)* m.Inverse();
}

matrix& matrix::operator=(const matrix& m)
{
    row = m.row;
    column = m.column;
    data = new double[row*column];
    for(int i=0; i<row*column; i++)
        data[i] = m.data[i];  return *this;
}

matrix& matrix::operator+=(const matrix& m)
{
    if(row != m.row || column != m.column)
        exit(0);  for(int i=0; i<row*column; i++)
        data[i] += m.data[i];  return *this;
}

matrix& matrix::operator-=(const matrix& m)
{
    if(row != m.row || column != m.column)
        exit(0);
    for(int i=0; i<row*column; i++)
        data[i] -= m.data[i];  return *this;
}

matrix matrix::operator[](int i) //take one row
{
    int count = 0;
    matrix temp;
    temp.row = 1;
    temp.column = column;
    temp.data = new double[column];
    for(int j=(i-1)*column; j<i*column; j++)
        temp.data[count++] = data[j];
    return temp;
}

matrix matrix::operator()(int i) //take one column
{
    int count = 0;
    matrix temp;
    temp.row = row;
    temp.column = 1;
    temp.data = new double[row];
    for(int j=i-1; j<(row-1)*column+i; j+=column)
        temp.data[count++] = data[j];
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

matrix matrix::TakeSubBlock(int i1, int j1, int i2, int j2) //take one sub block
{
    if(i1 > i2 || j1 > j2)
        exit(0);  int count = 0;
    matrix temp;  temp.row = (i2-i1)+1;
    temp.column = (j2-j1)+1;
    temp.data = new double[temp.row * temp.column];
    for(int i=(i1-1)*column+j1-1; i<(i2-1)*column+j2; i++)
    {
        temp.data[count++] = data[i];
        if(count % temp.column == 0)
            i += (column+j1-j2-1);
    }
    return temp;
}

matrix matrix::Rotate()
{
    matrix temp;
    temp.row = column;
    temp.column = row;
    temp.data = new double[row*column];
    for(int i=0; i<row; i++)
        for(int j=0; j<column; j++)
            temp.data[j*row+i] = data[i*column+j];
    return temp;
}

double detmatrix(const matrix M)
{
    if(M.row != M.column)
        exit(0);
    int m, n, s, t, k=1;
    double c, x, sn, f=1;
    for(int i=0,j=0; i<M.row&&j<M.row; i++,j++)
    {
        if(M.data[i*M.column+j] == 0)
        {
            for(m=i; M.data[m*M.column+j]==0&&m<M.row; m++);
            if(m == M.row)
            {
                sn = 0;
                return sn;
            }
            else
                for(n=j; n<M.row; n++)
                {
                    c = M.data[i*M.column+n];
                    M.data[i*M.column+n] = M.data[m*M.column+n];
                    M.data[m*M.column+n] = c;
                }
            k *= (-1);
        }
        for(s=M.row-1; s>i; s--)
        {
            x = M.data[s*M.column+j];
            for(t=j; t<M.row; t++)
                M.data[s*M.column+t] -= M.data[i*M.column+t]*(x / M.data[i*M.column+j]);
        }
    }
    for(int i=0; i<M.row; i++)
        f *= M.data[i*M.column+i];
    sn = k * f;
    return sn;
}

matrix matrix::Inverse()
{
    int count = 0;
    double var = 0;
    double var1;
    if((var = detmatrix(*this)) == 0)
        exit(0);  matrix temp1, temp2;
    temp1.row = row;
    temp1.column = column;
    temp1.data = new double[row*column];
    temp2.row = row - 1;
    temp2.column = column - 1;
    temp2.data = new double[(row-1)*(column-1)];
    for(int i=0; i<row; i++)
        for(int j=0; j<column; j++)
        {
            int l;
            if((i+j)%2 == 0)
                l = 1;
            else
                l = -1;
            for(int k=0; k<(row-1)*(column-1);)
            {
                if(count/row == i || count%row == j)
                {
                    count ++;
                    continue;
                }
                temp2.data[k] = data[count];
                k ++;
                count ++;
            }
            count = 0;
            var1 = detmatrix(temp2)*l/var;
            temp1.data[j*row+i] = var1 ? var1:0;
        }
    return temp1;
}

matrix GaussEliminationSeq(const matrix ar,const matrix br)
{

    if(ar.row!=ar.column||ar.row!=br.row||br.column!=1)
        throw(1);

    int r=ar.row,c=ar.column;
    double *data = new double[r*c];
    for(int i=0; i<r*c; i++)
        data[i] = ar.data[i];
    double *b= new double[c];
    for(int i=0;i<c;i++)
        b[i]=br.data[i];
    double m=0;

    for(int k=0;k<r-1;k++)//消元过程
    {
        if(data[k*c+k]==0)throw(2);
        for(int i=k+1;i<r;i++)
        {
            m=data[i*c+k]/data[k*c+k];

            for(int j=k;j<c;j++)
            {
                data[i*c+j]=data[i*c+j]-data[k*c+j]*m;
            }
            b[i]=b[i]-b[k]*m;
        }
    }
/* 回代过程*/
    b[c-1]=b[c-1]/data[c*c-1];
    for(int i=c-2;i>=0;i--)
    {
        m=0;
        for(int j=c-1;j>i;j--)
        {
            m+=data[i*c+j]*b[j];
        }
        b[i]=(b[i]-m)/data[i*c+i];
    }
    matrix ans=matrix(r,1,b);
    delete[] data;
    delete[] b;
    return ans;
}

matrix GaussEliminationColPrin(const matrix ar,const matrix br)
{
    if(ar.row!=ar.column||ar.row!=br.row||br.column!=1)
        throw(1);

    int r=ar.row,c=ar.column;
    double *data = new double[r*c];
    for(int i=0; i<r*c; i++)
        data[i] = ar.data[i];
    double *b= new double[c];
    for(int i=0;i<c;i++)
        b[i]=br.data[i];
    double m=0;

    bool changable=0;
    int change_row=0;
    double max_num=0;
    double p;

    for(int k=0;k<r-1;k++)
    {
        changable=0;
        change_row=k;
        max_num=data[k*c+k];
        for(int i=k+1;i<r;i++)
        {
            if(max_num<data[i*c+k])
            {
                changable=1;
                change_row=i;
                max_num=data[i*c+k];
            }
        }

        if(changable)
        {
            for(int i=0;i<c;i++)
            {
                p=data[k*c+i];
                data[k*c+i]=data[change_row*c+i];
                data[change_row*c+i]=p;
            }
            p=b[k];
            b[k]=b[change_row];
            b[change_row]=p;
        }
    }

    for(int k=0;k<r-1;k++)//消元过程
    {
        if(data[k*c+k]==0)throw(2);
        for(int i=k+1;i<r;i++)
        {
            m=data[i*c+k]/data[k*c+k];

            for(int j=k;j<c;j++)
            {
                data[i*c+j]=data[i*c+j]-data[k*c+j]*m;
            }
            b[i]=b[i]-b[k]*m;
        }
    }
/* 回代过程*/
    b[c-1]=b[c-1]/data[c*c-1];
    for(int i=c-2;i>=0;i--)
    {
        m=0;
        for(int j=c-1;j>i;j--)
        {
            m+=data[i*c+j]*b[j];
        }
        b[i]=(b[i]-m)/data[i*c+i];
    }
    matrix ans=matrix(r,1,b);
    delete[] data;
    delete[] b;
    return ans;
}

matrix DoolittleDecomposition(const matrix ar,const matrix br)
{
    if(ar.row!=ar.column||ar.row!=br.row||br.column!=1)
        throw(1);

    int r=ar.row,c=ar.column;
    double *data = new double[r*c];
    for(int i=0; i<r*c; i++)
        data[i] = ar.data[i];
    double *b= new double[c];
    for(int i=0;i<c;i++)
        b[i]=br.data[i];
    double m;
    //分解过程
    for(int k=0;k<r;k++)
    {
        for(int j=k;j<r;j++)
        {
            m=0;
            for(int t=0;t<=k-1;t++)
            {
                m+=data[k*c+t]*data[t*c+j];
            }
            data[k*c+j]-=m;
        }
        for(int i=k+1;i<r;i++)
        {
            m=0;
            for(int t=0;t<=k-1;t++)
            {
                m+=data[i*c+t]*data[t*c+k];
            }
            data[i*c+k]=(data[i*c+k]-m)/data[k*c+k];
        }
    }
    matrix ans=matrix(r,c,data);
    cout<<ans;
    //求解过程
    for(int i=0;i<c;i++)
    {
        for(int j=0;j<i;j++)
            b[i]-=(data[i*c+j]*b[j]);
    }

    for(int i=c-1;i>=0;i--)
    {
        m=0;
        for(int j=i+1;j<c;j++)
            m+=data[i*c+j]*b[j];
        b[i]=(b[i]-m)/data[i*c+i];
    }
     ans=matrix(r,1,b);
    delete[] data;
    delete[] b;
    return ans;
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

double matrix::max_eig(double error, double trans)
{
    matrix a=*this;
    double eta,beta1=0,beta2=0;
    if(trans!=0)
    {
        for(int i=0;i<column;i++)
           a.data[i*column+i]=data[i*column+i]-trans;
    }
//    cout<<*this<<a;
    double *b=new double[row];
    memset(b,0,row*sizeof(double));
    b[0]=1;
    matrix u(row,1,b);
    eta=u.norm('f');
    matrix y;
    y=u*(1.0/eta);
    u=a*y;
    for(int i=0;i<row;i++)
        beta2+=y.TakeElem(i,0)*u.TakeElem(i,0);
    do
    {
        beta1=beta2;
        beta2=0;
        eta=u.norm('f');
        y=u*(1.0/eta);
        u=a*y;
        for(int i=0;i<y.row;i++)
            beta2+=y.TakeElem(i,0)*u.TakeElem(i,0);
        //cout<<"new beta"<<beta1<<" "<<beta2<<endl;
    }while(fabs((beta2-beta1)/beta2)>=error);
    delete[] b;
    return beta2;
}

double matrix::min_eig(double error,double trans)
{
    matrix a=*this;
    double eta,beta1=0,beta2=0;
    if(trans!=0)
    {
        for(int i=0;i<column;i++)
           a.data[i*column+i]=data[i*column+i]-trans;
    }
//    cout<<*this<<a;
    double *b=new double[row];
    memset(b,0,row*sizeof(double));
    b[0]=1;
    matrix u(row,1,b);
    eta=u.norm('f');
    matrix y;
    y=u*(1.0/eta);
    u=GaussEliminationColPrin(a,y);
    for(int i=0;i<row;i++)
        beta2+=y.data[i]*u.data[i];
    do
    {
        beta1=beta2;
        beta2=0;
        eta=u.norm('f');
        y=u*(1.0/eta);
        u=GaussEliminationColPrin(a,y);
        for(int i=0;i<y.row;i++)
            beta2+=y.TakeElem(i,0)*u.TakeElem(i,0);
        //cout<<"new beta"<<beta1<<" "<<beta2<<endl;
    }while(fabs((1.0-beta1/beta2))>=error);
    delete[] b;
    return 1.0/beta2;
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
