#include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
using namespace std;

double p[56][56],u[55][56],ru[55][56],ru2[55][56],v[56][55],rv[56][55],rv2[56][55];
double p_[56][56],p1[56][56]; 
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
	
	
	for(int t=0;t<500000;t++)
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
		ru[0][j] = ru_[0][j];
		u[0][j] = ru[0][j]/rho; 
		ru2[0][j] = ru[0][j]*u[0][j];	
	
		ru_[54][j] = 0;
		ru[54][j] = ru_[54][j];
		u[54][j] = ru[54][j]/rho;
		ru2[54][j] = ru[54][j]*u[54][j]; 
	
	
	
	
	
	 } 
	
	
	//上下边界
	for(i=0;i<=54;i++) 
	{
		ru_[i][0] = 0;
		ru[i][0] = ru_[i][0];
		u[i][0] = ru[i][0]/rho; 
		ru2[i][0] = ru[i][0]*u[i][0];
		
		
		
		ru_[i][55] = 1;
		ru[i][55] =ru_[i][55];
		u[i][55] = ru[i][55]/rho; 
		ru2[i][55] = ru[i][55]*u[i][55];
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
		rv[0][j] = rv_[0][j];
		v[0][j] = rv[0][j]/rho;
		rv2[0][j] = rv[0][j]*v[0][j];
		
		
		rv_[55][j] = 0;
		rv[55][j] =rv_[55][j];
		v[55][j] = rv[55][j]/rho; 
		rv2[55][j] = rv[55][j]*v[55][j];
	 } 
	
	//上下边界
	for(i=0;i<=55;i++)
	{
	
	rv_[i][0] = 0;
	rv[i][0] = rv_[i][0];
	v[i][0] = rv[i][0]/rho;
	rv2[i][0] = rv[i][0]*v[i][0];
	
	rv_[i][54] = 0; 
	rv[i][54] = rv_[i][54];
	v[i][54] = rv[i][54]/rho;
	rv2[i][54] = rv[i][54]*v[i][54];
	
	
	
   }
//=============================================================================//
//压力修正(内部流场点) 


for(i=1;i<55;i++)
{
for(j=1;j<55;j++)
{

a = 2*((dt)/(dx*dx)+(dt)/(dy*dy));
b = -1*(dt/(dx*dx));
c = -1*(dt/(dy*dy));
d = (1/dx)*(ru[i][j]-ru[i-1][j])+(1/dy)*(rv[i][j]-rv[i][j-1]);

p_[i][j] = (-1/a)*(b*p[i+1][j]+b*p[i-1][j]+c*p[i][j+1]+c*p[i][j-1]+d);
}
}

for(i=1;i<55;i++)
{
for(j=1;j<55;j++)
{

p1[i][j] = p[i][j]+1e-5*p_[i][j];
}
}

//边界上的压力

//左右边界

for(j=1;j<55;j++)
{
	p1[0][j] = p1[1][j];
	p1[55][j] = p1[54][j];
 } 

//上下边界

for(i=0;i<=55;i++)
{
	p1[i][0] = p1[i][1];
	p1[i][55] = p1[i][54];
	
 } 


//更新 p 值

for(i=0;i<=55;i++)
{
	
	for(j=0;j<=55;j++)
	{
		
		p[i][j] = p1[i][j];
		
	 } 
	
 } 

//利用更新的 p值重新计算 u、v; 

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
		ru[0][j] = ru_[0][j];
		u[0][j] = ru[0][j]/rho; 
		ru2[0][j] = ru[0][j]*u[0][j];	
	
		ru_[54][j] = 0;
		ru[54][j] = ru_[54][j];
		u[54][j] = ru[54][j]/rho;
		ru2[54][j] = ru[54][j]*u[54][j]; 
	
	
	
	
	
	 } 
	
	
	//上下边界
	for(i=0;i<=54;i++) 
	{
		ru_[i][0] = 0;
		ru[i][0] = ru_[i][0];
		u[i][0] = ru[i][0]/rho; 
		ru2[i][0] = ru[i][0]*u[i][0];
		
		
		
		ru_[i][55] = 1;
		ru[i][55] =ru_[i][55];
		u[i][55] = ru[i][55]/rho; 
		ru2[i][55] = ru[i][55]*u[i][55];
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
		rv[0][j] = rv_[0][j];
		v[0][j] = rv[0][j]/rho;
		rv2[0][j] = rv[0][j]*v[0][j];
		
		
		rv_[55][j] = 0;
		rv[55][j] =rv_[55][j];
		v[55][j] = rv[55][j]/rho; 
		rv2[55][j] = rv[55][j]*v[55][j];
	 } 
	
	//上下边界
	for(i=0;i<=55;i++)
	{
	
	rv_[i][0] = 0;
	rv[i][0] = rv_[i][0];
	v[i][0] = rv[i][0]/rho;
	rv2[i][0] = rv[i][0]*v[i][0];
	
	rv_[i][54] = 0; 
	rv[i][54] = rv_[i][54];
	v[i][54] = rv[i][54]/rho;
	rv2[i][54] = rv[i][54]*v[i][54];
	
	
	
   }


	
}
	
	
	ostringstream name;
	  name<<"cavity_"<<6<<".dat";
	  ofstream out(name.str().c_str());
	  out<< "Title= \"LBM Lid Driven Flow\"\n" << "VARIABLES=\"X\",\"Y\",\"v\"\n" << "ZONE T=\"BOX\",I=" << 56 << ",J=" << 56 << ",F=POINT" << endl;
	  for(int j=0;j<=55;j++)
	     for(int i=0;i<=55;i++)
	     {
	     	out<<i*dx<<" "<<j*dx<<" "<<v[i][j]<<endl;
		 }
	
	
	
	
	
	return 0; 
 } 
  
