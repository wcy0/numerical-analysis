#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<iomanip>
#include<time.h>

#define N 10               /*定义矩阵的阶数*/
#define E 1.0e-12           /*定义全局变量相对误差限*/
#define L 20               /*迭代次数*/
FILE *fp;
struct C                   /*定义结构体*/
{
  double R;  double I;

};

void shuchu(double a[N][N])         /*输出一个10*10的矩阵*/
{
  int i,j;
  double b=0;
  for(i=0;i<N;i++)
  {
   for(j=0;j<N;j++)
   {
    if(fabs(a[i][j])<E)




printf("%+.12e\t",b);
     else

printf("%+.12e\t",a[i][j]);

  }

  printf("\n");

 }

}


int sign(double x)                 /*符号函数*/
{
  if(x<0) return -1;
  if(x == 0) return 0;
return 1;
}


void assignment(double a[N][N])           /*输入矩阵A*/
{
  int i,j;

 for(i=0;i<N;i++)
  {
  for(j=0;j<N;j++)
   {

    if(i==j)

a[i][j]=1.52*cos((i+1)+1.2*(j+1));

   else

a[i][j]=sin(0.5*(i+1)+0.2*(j+1));

  }

 }

}

void Hessenberg(double a[N][N])          /*将矩阵A拟上三角化*/
{
  int r,i,j,b;
  double d,c,h,t;

 double p[N],q[N],u[N],w[N];


  for(r=0;r<N-2;r++)
  {

      b=0;

  for(i=r+2;i<N;i++)
      {
         if(fabs(a[i][r])>E)

b++;

  }

  if(b==0) continue;
   else
   {
    d=0;

   for(i=r+1;i<N;i++)

d=d+a[i][r]*a[i][r];
    d=sqrt(d);

   if(a[r+1][r]==0) c=d;
    else c=(-1)*sign(a[r+1][r])*d;
    h=c*c-c*a[r+1][r];
    for(i=0;i<=r;i++)

u[i]=0;

   u[r+1]=a[r+1][r]-c;
    for(i=r+2;i<N;i++)




u[i]=a[i][r];
    for(i=0;i<N;i++)
    {
     p[i]=0;

    for(j=r+1;j<N;j++)

p[i]=p[i]+a[j][i]*u[j];

    p[i]=p[i]/h;

   }

   for(i=0;i<N;i++)
    {
     q[i]=0;

    for(j=r+1;j<N;j++)

q[i]+=a[i][j]*u[j];
     q[i]=q[i]/h;

   }
    t=0;

   for (i=r+1;i<N;i++)

t+=p[i]*u[i];
    t=t/h;

   for(i=0;i<N;i++)

w[i]=q[i]-t*u[i];

   for(i=0;i<N;i++)
    for(j=0;j<N;j++)

a[i][j]=a[i][j]-w[i]*u[j]-u[i]*p[j];

  }

 }

}


void QRfenjie(double a[N][N],int m,int n)   /*将m维(m≤10)矩阵QR分解*/


{                                    /*当n=0时用于基本QR分解；n≠0时用于带双步位移的QR方法*/

 int r,i,j,k;
  double d,c,h,t,s,s1;

 double p[N],q[N],u[N],w[N],w1[N],v[N];
 double b[N][N]={{0}},g[N][N]={{0}};
  if(n==0) goto loop1;
  s=a[m-2][m-2]+a[m-1][m-1];

 s1=a[m-2][m-2]*a[m-1][m-1]-a[m-1][m-2]*a[m-2][m-1];
  for(i=0;i<m;i++)
  for(j=0;j<m;j++)
  {
   for(r=0;r<m;r++)

b[i][j]+=a[i][r]*a[r][j];

  b[i][j]=b[i][j]-s*a[i][j]+s1*(i==j);

 }


goto loop2;

loop1:  for(i=0;i<m;i++)
  for(j=0;j<m;j++)
  {
   b[i][j]=a[i][j];
   g[i][j]=(i==j);


}

loop2: for(r=0;r<m-1;r++)
 {
   k=0;

   for(i=r+1;i<m;i++)
    {
       if (fabs(b[i][r])>E)

k++;




}

    if(k==0) continue;
    d=0;

   for(i=r;i<m;i++)
    {
     d +=b[i][r]*b[i][r];

   }
    d=sqrt(d);
    if(b[r][r]==0) c=d;
    else c=(-1)*sign(b[r][r])*d;
    h=c*c-c*b[r][r];
    for(i=0;i<r;i++)

u[i]=0;

   u[r]=b[r][r]-c;
    for(i=r+1;i<m;i++)

u[i]=b[i][r];

   for(i=0;i<m;i++)
    {
     v[i]=0;
     for(j=r;j<m;j++)

v[i] +=b[j][i]*u[j];

    v[i]=v[i]/h;

   }

   for(i=0;i<m;i++)
    for(j=0;j<m;j++)

b[i][j]=b[i][j]-u[i]*v[j];
    for(i=0;i<m;i++)
    {
     p[i]=0;
     for(j=r;j<m;j++)





p[i] +=a[j][i]*u[j];
    p[i]=p[i]/h;

   }

   for(i=0;i<m;i++)
    {
     q[i]=0;
     for(j=r;j<m;j++)

q[i] +=a[i][j]*u[j];

    q[i]=q[i]/h;

   }
    t=0;

   for(i=r;i<m;i++)

t+=p[i]*u[i];
    t=t/h;

   for(i=0;i<m;i++)

w[i]=q[i]-t*u[i];

   for(i=0;i<m;i++)
    for(j=0;j<m;j++)

a[i][j]=a[i][j]-w[i]*u[j]-u[i]*p[j];

   if(n!=0) goto loop3;
    for(i=0;i<m;i++)
    {
        w1[i]=0;
     for(j=r;j<m;j++)

w1[i]+=g[i][j]*u[j];

   }

   for(i=0;i<m;i++)
    for(j=0;j<m;j++)


g[i][j]=g[i][j]-w1[i]*u[j]/h;

loop3:;

}

   if (n!=0) goto loop4;
   printf( "Q阵：\n");
   shuchu(g);

  printf("R阵：\n");
   shuchu(b);

  printf("RQ阵：\n");

shuchu(a);

 loop4:;

}

void Gauss(double d) /*用列主元素高斯消去法求属于特征值d的特征向量，并输出*/
{

  int i,j,k,t;
  double a[N][N],x[N]={0},b,m;
  assignment(a);
  for(i=0;i<N;i++)

a[i][i]-=d;

 for(k=0;k<N-1;k++)
  {
   t=k;

  for(i=k+1;i<N;i++)
   if(fabs(a[i][k])>fabs(a[t][k]))

t=i;

  for(j=k;j<N;j++)
   {
    b=a[k][j];
    a[k][j]=a[t][j];
    a[t][j]=b;

  }



for(i=k+1;i<N;i++)
{
    m=a[i][k]/a[k][k];
    for(j=k+1;j<N;j++)
    a[i][j]-=m*a[k][j];
  }
}
 x[N-1]=1;

 for(k=N-2;k>=0;k--)
  {
   for(j=N-1;j>k;j--)

   x[k]-=a[k][j]*x[j];

   x[k]=x[k]/a[k][k];

 }
  b=0;

  for(i=0;i<N;i++)                      /*将特征向量归一化*/
    b+=x[i]*x[i];
  b=sqrt(b);
  for(i=0;i<N;i++)
  x[i]=x[i]/b;

  printf("属于实特征值%+.12e的特征向量为:\n",d);
  for(i=0;i<N;i++)
    printf("%+.12e\n",x[i]);

}

void QRfa(double a[N][N])  /*带双步位移的QR方法求矩阵的所有特征值*/
{                          /*并输出所有实特征值的特征向量*/
  int i,k=1,m=N;
  double b,c,d;
  struct C y[N];

for (i=0;i<N;i++)

{
   y[i].R=0;
   y[i].I=0;


}

loop2: if(fabs(a[m-1][m-2])<=E)
{
    y[m-1].R=a[m-1][m-1];
    m-=1;


goto loop4;

}

else goto loop1;
  loop4: if(m==1)
  {
     y[0].R=a[0][0];


goto loop3;

 }
    else
    {
     if(m>1) goto loop2;
     else goto loop3;


}
 loop1: b=a[m-2][m-2]+a[m-1][m-1];
    c=a[m-2][m-2]*a[m-1][m-1]-a[m-1][m-2]*a[m-2][m-1];
    d=b*b-4*c;
    if (d>=0)
    {
     d=sqrt(d);
     y[m-1].R=(b+d)/2;
     y[m-2].R=(b-d)/2;




}
else

{
     d=sqrt(fabs(d));
     y[m-1].R=b/2;
     y[m-2].R=b/2;
     y[m-1].I=d/2;
     y[m-2].I=(-1)*d/2;
}

   if (m==2) goto loop3;
    if (fabs(a[m-2][m-3])<=E)
    {
     m=m-2;
     goto loop4;
    }
    if (k==L)
    {
     printf("未得到A的全部特征值\n");
     goto loop5;
    }

   QRfenjie(a,m,1);
    k++;

goto loop2;
loop3: printf("A的全部特征值已计算完毕\n");
    for (i=0;i<N;i++)

printf("%+.12e+i*(%+.12e)\n",y[i].R,y[i].I);

   for (i=0;i<10;i++)
    if (y[i].I==0)


Gauss(y[i].R);


loop5:;

}
int main()
{
  double A[N][N];
  assignment(A);
  Hessenberg(A);
  printf("拟上三角化后所得的矩阵A(n-1)：\n");
  shuchu(A);
  QRfenjie(A, N, 0);
  assignment(A);
  Hessenberg(A);
  QRfa(A);

}





