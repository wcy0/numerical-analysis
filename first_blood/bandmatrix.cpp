/*
 *
 * author：X_Y 耀
 *
 * 2017/11/5
 *
 */

#include "bandmatrix.h"

bandmatrix::bandmatrix()
{
    row  = 0;
    column = 0;
    data = NULL;
    width=0;
    width_be=0;
    width_up=0;
}

bandmatrix::bandmatrix(int r, int c, int w, int u, double* p):row(r),column(c),width(w),width_up(u),width_be(w-u-1)
{
    data = new double[r*w];
    for(int i=0; i<r*w; i++)
        data[i] = p[i];
}

bandmatrix::~bandmatrix()
{
     if(data != NULL)
    {
        delete[] data;
        data = NULL;
    }
}

bandmatrix::bandmatrix(const bandmatrix& m):matrix()
{
    row    = m.row;  column = m.column;
    width  = m.width; width_be = m.width_be; width_up=m.width_up;
    data = new double[row * width];
    for(int i=0; i<row*width; i++)
        data[i] = m.data[i];
}

bandmatrix& bandmatrix::operator=(const bandmatrix& m)
{
    row    = m.row;  column = m.column;
    width  = m.width; width_be = m.width_be; width_up=m.width_up;
    data = new double[row * width];
    for(int i=0; i<row*width; i++)
        data[i] = m.data[i];
    return *this;
}

ostream& operator<<(ostream& o,bandmatrix& m)//overload the cout
{
    double* p;
    p=m.getData();
    o<<m.getRow()<<" "<<m.getColumn()<<" "<<m.getWidth()<<" "<<m.getWidthUp()<<" "<<endl;
    for(int i=0;i<m.getColumn()*m.getWidth();i++,p++)
    {
        o<<*p<<" ";
        if(i%m.getColumn()==0)o<<endl;
    }
    return o;
}

double bandmatrix::TakeElem(int i, int j)
{
    if(i-j<=width_be&&j-i<=width_up)return data[(i+width_up-j)*column+j];
    else return 0;
}

void bandmatrix::SetElem(int i,int j, int b)
{
    if(i-j<=width_be&&j-i<=width_up)
        data[(i+width_up-j)*column+j]=b;
}

matrix bandmatrix::operator *(matrix &m)
{
    if(column != m.getRow())  exit(0);
    double *p=new double[row*m.getColumn()];
    memset(p, 0, sizeof(double)*row*m.getColumn());
    matrix temp(row,m.getColumn(),p);
    for(int i=0; i<row; i++)
    {
        for(int j=0;j<m.getColumn();j++)
        {
            for(int k=0;k<column;k++)
                temp.SetElem(i,j,temp.TakeElem(i,j)+this->TakeElem(i,k)*m.TakeElem(k,j));
        }
    }
    return temp;
}

bandmatrix bandmatrix::DoolittleLU()
{
    int n=column,s=width_up,r=width_be;
    double m;
    bandmatrix a=*this;
    for(int k=0;k<n;k++)
    {
        for(int j=k;j<n&&j<=k+s;j++)
        {
            m=0;
            int t=k-r>0?k-r:0;
            if(t<j-s)t=j-s;
            for(;t<=k-1;t++)
            {
                m+=a.data[(k+s-t)*n+t]*a.data[(t+s-j)*n+j];
            }
            a.data[(k+s-j)*n+j]-=m;
        }

        for(int i=k+1;i<n&&i<=k+r;i++)
        {
            m=0;
            int t=i-r>1?i-r:1;
            if(t<k-s)t=k-s;
            for(;t<=k-1;t++)
            {
                m+=a.data[(i-t+s)*n+t]*a.data[(t-k+s)*n+k];
            }
            a.data[(i-k+s)*n+k]=(a.data[(i-k+s)*n+k]-m)/a.data[s*n+k];
        }
    }
    return a;
}

matrix DoolittleDecomposition(bandmatrix ar,matrix br)
{
    int n=ar.column,s=ar.width_up,r=ar.width_be;
    double m,*b;
    b=new double[n];
    for(int i=0;i<n;i++)
        b[i]=br.TakeElem(i,0);
    matrix ans;
    bandmatrix lu=ar.DoolittleLU();
    for(int i=1;i<n;i++)
    {
        for(int t=i-r>0?i-r:0;t<=i-1;t++)
            b[i]-=(lu.data[(i-t+s)*n+t]*b[t]);
    }

    b[n-1]/=(lu.data[s*n+n-1]);
    for(int i=n-2;i>=0;i--)
    {
        m=0;
        for(int t=i+1;t<n&&t<=i+s;t++)
            m+=lu.data[(i-t+s)*n+t]*b[t];
        b[i]=(b[i]-m)/lu.data[s*n+i];
    }
    ans=matrix(n,1,b);
    delete[] b;
    return ans;
}

double bandmatrix::max_eig(double error, double trans)
{
    bandmatrix a=*this;
    int s=width_up;
    double eta,beta1=0,beta2=0;
    if(trans!=0)
    {
        for(int i=0;i<column;i++)
           a.data[s*column+i]=data[s*column+i]-trans;
    }
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
        for(int i=0;i<y.getRow();i++)
            beta2+=y.TakeElem(i,0)*u.TakeElem(i,0);
    }while(fabs((beta2-beta1)/beta2)>=error);
    delete[] b;
    return beta2;
}

double bandmatrix::min_eig(double error, double trans)
{
    bandmatrix a=*this;
    int s=width_up;
    double eta,beta1=0,beta2=0;
    if(trans!=0)
    {
        for(int i=0;i<column;i++)
           a.data[s*column+i]=data[s*column+i]-trans;
    }
    double *b=new double[row];
    memset(b,0,row*sizeof(double));
    b[0]=1;
    matrix u(row,1,b);
    eta=u.norm('f');
    matrix y;
    y=u*(1.0/eta);
    u=DoolittleDecomposition(a,y);
    for(int i=0;i<y.getRow();i++)
        beta2+=y.TakeElem(i,0)*u.TakeElem(i,0);
    do
    {
        beta1=beta2;
        beta2=0;
        eta=u.norm('f');
        y=u*(1.0/eta);
        u=DoolittleDecomposition(a,y);
        for(int i=0;i<y.getRow();i++)
            beta2+=y.TakeElem(i,0)*u.TakeElem(i,0);
    }while(fabs((1.0-beta1/beta2))>=error);
    delete[] b;
    return 1.0/beta2;
}

matrix bandmatrix::band2mat()
{
    double *a;
    a=new double[row*column];
    for(int i=0;i<row;i++)
        for(int j=0;j<column;j++)
            a[i*column+j]=TakeElem(i,j);
    matrix ans(row,column,a);
    delete[] a;
    return ans;
}

double bandmatrix::detmatrix()
{
    bandmatrix a= this->DoolittleLU();
    double ans=1,*p=NULL;
    p=a.data+width_up*column;
    for(int i=0;i<column;i++,p++)
        ans *= *p;
    return ans;
}
