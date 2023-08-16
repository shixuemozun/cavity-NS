#include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
using namespace std;

double p[56][56],u[55][56],ru[55][56],ru2[55][56],v[56][55],rv[56][55],rv2[56][55];

double ru_[55][56],rv_[56][55];
int i,j;
double u_,un,v_,vn;

double rho = 1;
double niu = 1; 
double A,B;
double a,b,c,d;
double dx,dy,dt;

int main()
{
	
	dx = 1;

	dy = 1;

	dt = 0.01;
	rho = 1;
	niu = 1;
	
	
	
	
	
	for(i=0;i<=54;i++)
	{
		for(j=0;j<=55;j++)
		
	{
		u[i][j] = 0;
		ru[i][j] = u[i][j]*rho;
		ru2[i][j] = ru[i][j]*u[i][j];
		u[i][55] = 1;
		ru[i][55] = u[i][55]*rho;
		ru2[i][55] = ru[i][55]*u[i][55];
		
	}
		
	}
	
	for(i=0;i<=55;i++)
	{
		for(j=0;j<=54;j++)
		
	{
		v[i][j] = 0;
		rv[i][j] = v[i][j]*rho;
		rv2[i][j] = rv[i][j]*v[i][j];
		
	}
		
	}
	
	for(i=0;i<=55;i++)
	{
		for(j=0;j<=55;j++)
		{
			p[i][j] = 0;
			
		}
		
		
	}
	
	
	for(int t=0;t<1000;t++)
	{
	
	
	
		for(i=1;i<54;i++)
	{
		for(j=1;j<55;j++)
		{
		
		v_ = 1/2*(v[i][j]+v[i+1][j]);
		 vn = 1/2*(v[i][j-1]+v[i+1][j-1]);
		A = -1*(((ru2[i+1][j]-ru2[i-1][j])/(2*dx))+((ru[i][j+1]*v_-ru[i][j-1]*vn)/(2*dy)))+niu*((u[i+1][j]-2*u[i][j]+u[i-1][j])/(dx*dx)+(u[i][j+1]-2*u[i][j]+u[i][j-1])/(dy*dy));	
		ru_[i][j] = ru[i][j] + A*dt-(dt/dx)*(p[i+1][j]-p[i][j]);	//ru在n+1时刻的数值 
			
		}
				
	}
	
	//更新变量 Ux 
	
	for(i=1;i<54;i++)
	{
		for(j=1;j<55;j++)
		{
			ru[i][j] = ru_[i][j];
			u[i][j] = ru[i][j]/rho;
			ru2[i][j] = ru[i][j]*u[i][j];
			// v值还没有更新 
			
		}
		
		
	 } 
	
	//左右边界 
	for(j=1;j<55;j++)
	{
		ru_[0][j] = 0;
		ru_[54][j] = 0;
	 } 
	
	
	//上下边界
	for(i=0;i<=54;i++) 
	{
		ru_[i][0] = 0;
		ru_[i][55] = 1;
		
	}
	
	//========================================================================================//
	//计算 Uy 
	
	for(i=1;i<55;i++)
	{
		for(j=1;j<54;j++)
		{
		u_ = 1/2*(u[i][j+1]+u[i][j]); 
	    un = 1/2*(u[i-1][j+1]+u[i-1][j]);
	    
	    B = -1*((rv[i+1][j]*u_-rv[i-1][j]*un)/(2*dx)+(rv2[i][j+1]-rv2[i][j-1])/(2*dy))+niu*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dx*dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/(dy*dy));
		
		rv_[i][j] = rv[i][j]+B*dt-(dt/dy)*(p[i][j+1]-p[i][j]);
			
		}
		
		
	} 
	
	//更新变量 Uy
	
	 for(i=1;i<55;i++)
	 {
	 	for(j=1;j<54;j++)
	 	{
	 		rv[i][j] = rv_[i][j];
			v[i][j] = rv[i][j]/rho;
			rv2[i][j] = rv[i][j]*v[i][j];
	 		
	 		
		 }
	 	
	 	
	  } 
	
	//左右边界 
	
	for(j=1;j<54;j++)
	{
		rv_[0][j] = 0;
		rv_[55][j] = 0;
		
	 } 
	
	//上下边界
	for(i=0;i<=55;i++)
	{
	
	rv_[i][0] = 0;
	rv_[i][54] = 0; 
   }
	
}
	
	ostringstream name;
	  name<<"cavity_"<<6<<".dat";
	  ofstream out(name.str().c_str());
	  out<< "Title= \"LBM Lid Driven Flow\"\n" << "VARIABLES=\"X\",\"Y\",\"u\"\n" << "ZONE T=\"BOX\",I=" << 55 << ",J=" << 56 << ",F=POINT" << endl;
	  for(int j=0;j<=55;j++)
	     for(int i=0;i<=54;i++)
	     {
	     	out<<i*dx<<" "<<j*dx<<" "<<u[i][j]<<endl;
		 }
	
	
	
	
	
	return 0; 
 } 
