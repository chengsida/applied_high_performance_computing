

#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;
void    f(double** vx, double** vy, double** h, double** kx, double** ky, double** kh);
void  matrixmult(double** k, double** m,  double c);
void matrixadd(double** k,double** h, double** m);
void final_matrixadd(double** a, double** b, double** c, double** d, double** e,double f);
void writeData(int k, double** h, double x_min, double y_min, double delta_x, double delta_y, int N_x, int N_y);



// Global variables
const double    x_min       =  0.00;
const double    x_max       = 100.00;
const double    t_min       =  0.00;
const double    t_max       = 100.00;
const double    y_min           =0.00;
const double    y_max           =100.00;
const double    delta_x     =1.00;
const double    delta_y         =1.00;
const double    delta_t     =0.2;
const int       N_x         =  (x_max-x_min)/delta_x+1;
const int       N_t         =  (t_max-t_min)/delta_t+1;
const int       N_y                     =  (y_max-y_min)/delta_y+1;
const double    g           =  9.81;

// Function declarations

int     main(int argc, char** argv)
{
	// Simulation parameters
	double t=0;
	int mm;
	int m;
	
	fstream     file;
	
	// Allocate arrays
	double** vx = new double* [N_x];
	vx[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		vx[m]=&vx[0][mm];
	}
	double** vy = new double* [N_x];
	vy[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		vy[m]=&vy[0][mm];
	}
	double** h = new double* [N_x];
	h[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		h[m]=&h[0][mm];
	}
	double** vxt = new double* [N_x];       //temp vx
	vxt[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		vxt[m]=&vxt[0][mm];
	}
	double** vyt = new double* [N_x];       //temp vy
	vyt[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		vyt[m]=&vyt[0][mm];
	}
	double** ht = new double* [N_x];        //temp h
	ht[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		ht[m]=&ht[0][mm];
	}
	double** kx1 = new double* [N_x];        // kx1 ky1 kh1
	kx1[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kx1[m]=&kx1[0][mm];
	}
	double** ky1 = new double* [N_x];
	ky1[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		ky1[m]=&ky1[0][mm];
	}
	double** kh1 = new double* [N_x];
	kh1[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kh1[m]=&kh1[0][mm];
	}
	double** kx2 = new double* [N_x];        // kx2 ky2 kh2
	kx2[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kx2[m]=&kx2[0][mm];
	}
	double** ky2 = new double* [N_x];        
	ky2[0]=new double [N_x*N_y];
	for( m=1,mm=N_y;m<N_x; m++,mm+=N_y)
	{
		ky2[m]=&ky2[0][mm];
	}
	double** kh2 = new double* [N_x];       
	kh2[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kh2[m]=&kh2[0][mm];
	}
	double** kx3 = new double* [N_x];         //kx3 ky3 kh3
	kx3[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kx3[m]=&kx3[0][mm];
	}
	double** ky3 = new double* [N_x];        
	ky3[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		ky3[m]=&ky3[0][mm];
	}
	double** kh3 = new double* [N_x];       
	kh3[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kh3[m]=&kh3[0][mm];
	}
	double** kx4 = new double* [N_x];               //kx4 ky4 kh4
	kx4[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kx4[m]=&kx4[0][mm];
	}
	double** ky4 = new double* [N_x];        
	ky4[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		ky4[m]=&ky4[0][mm];
	}
	double** kh4 = new double* [N_x];       
	kh4[0]=new double [N_x*N_y];
	for(  m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kh4[m]=&kh4[0][mm];
	}
	double** kxt = new double* [N_x];               //temp k values
	kxt[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kxt[m]=&kxt[0][mm];
	}
	double** kyt = new double* [N_x];        
	kyt[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kyt[m]=&kyt[0][mm];
	}
	double** kht = new double* [N_x];       
	kht[0]=new double [N_x*N_y];
	for( m=1, mm=N_y;m<N_x; m++,mm+=N_y)
	{
		kht[m]=&kht[0][mm];
	}
	// Set initial condition
	for(int i=0; i<N_x; i++)
	{
		for(int j=0; j<N_y;j++)
		{
			vx[i][j]=0;
			vy[i][j]=0;
			vxt[i][j]=0;
			vyt[i][j]=0;
			ht[i][j]=0;
			kx1[i][j]=0;
			ky1[i][j]=0;
			kh1[i][j]=0;
			kx2 [i][j]=0 ;
			ky2[i][j]=0 ;
			kh2[i][j]=0;
			kx3 [i][j]=0 ;
			ky3[i][j]=0 ;
			kh3[i][j]=0;
			kx4 [i][j]=0 ;
			ky4[i][j]=0 ;
			kh4[i][j]=0;
			kxt[i][j]=0 ;
			kyt[i][j]=0 ;
			kht[i][j]=0;
			double x=i*delta_x;
			double y=j*delta_y;
			h[i][j]=1+0.5*exp((-1/0.04)*(pow(x-30,2)+pow(y-30,2)));
		}
	}
	
	
	for(int l=0; l<N_t-1; l++)
	{
		t += delta_t;;
		f(vx,vy,h,kx1,ky1,kh1);
		matrixmult(kxt,kx1,delta_t/2);
		matrixmult(kyt,ky1,delta_t/2);
		matrixmult(kht,kh1,delta_t/2);
		matrixadd(vxt,vx,kxt);
		matrixadd(vyt,vy,kyt);
		matrixadd(ht,h,kht);
		f(vxt,vyt,ht,kx2,ky2,kh2);
		matrixmult(kxt,kx2,delta_t/2);
		matrixmult(kyt,ky2,delta_t/2);
		matrixmult(kht,kh2,delta_t/2);
		matrixadd(vxt,vx,kxt);
		matrixadd(vyt,vy,kyt);
		matrixadd(ht,h,kht);
		f(vxt,vyt,ht,kx3,ky3,kh3);
		matrixmult(kxt,kx3,delta_t);
		matrixmult(kyt,ky3,delta_t);
		matrixmult(kht,kh3,delta_t);
		matrixadd(vxt,vx,kxt);
		matrixadd(vyt,vy,kyt);
		matrixadd(ht,h,kht);
		f(vxt,vyt,ht,kx4,ky4,kh4);
		final_matrixadd(kxt,kx1,kx2,kx3,kx4,delta_t);
		final_matrixadd(kyt,ky1,ky2,ky3,ky4,delta_t);
		final_matrixadd(kht,kh1,kh2,kh3,kh4,delta_t);
		matrixadd(vx,vx,kxt);
		matrixadd(vy,vy,kyt);
		matrixadd(h,h,kht);
		
		
		
		
		
		
		
		
		if(l%50==0)
		{
			
			writeData(l, h, x_min, y_min, delta_x, delta_y,N_x, N_y);
			
		}			
		
		
		
	}
	file.close();
	
	
	
	// Deallocate arrays

	
	delete [] vx[0];
	delete [] vy[0];
	delete [] h[0];
	delete [] vxt[0];
	delete [] vyt[0];
	delete [] ht[0];
	delete [] kx1[0];
	delete [] ky1[0];
	delete [] kh1[0];
	delete [] kx2[0];
	delete [] ky2[0];
	delete [] kh2[0];
	delete [] kx3[0];
	delete [] ky3[0];
	delete [] kh3[0];
	delete [] kx4[0];
	delete [] ky4[0];
	delete [] kh4[0];
	delete [] kxt[0];
	delete [] kyt[0];
	delete [] kht[0];
	
	
	return 0;
	
}

void    f(double** vx, double** vy, double** h, double** kx, double** ky, double** kh)
{
	int x_a, x_b,x_c,x_d,y_e,y_f,y_g,y_h;
	
	for(int i=0; i<N_x; i++)
	{
		x_a=i-2;
		x_b=i-1;
		x_c=i+1;
		x_d=i+2;
		
		if(x_a==-2)
		{
			x_a=N_x-2;
		}
		else if (x_a==-1)
		{
			x_a=N_x-1;
		}
		else
		{
			x_a=i-2 ;
		}
		if(x_b==-1)
		{
			x_b=N_x-1;
		}
		if(x_c==N_x)
		{
			x_c=0;
		}
		if(x_d==N_x+1)
		{
			x_d=1;
		}
		else if(x_d==N_x)
		{
			x_d=0;
		}
		else
		{
			x_d=i+2;
		}
		
		for (int j=0;j<N_y;j++)
		{
			y_e=j-2;
			y_f=j-1;
			y_g=j+1;
			y_h=j+2;
			if(y_e==-2)
			{
				y_e=N_y-2;
			}
			else if (y_e==-1)
			{
				y_e=N_y-1;
			}
			if(y_f==-1)
			{
				y_f=N_y-1;
			}
			if(y_g==N_y)
			{
				y_g=0;
			}
			if(y_h==N_y+1)
			{
				y_h=1;
			}
			else if(y_h==N_y)
			{
				y_h=0;
			}
			else
			{
				y_h=j+2;
			}
			
			kx[i][j]=(((-1)*vx[i][j])/(12*delta_x))*(vx[x_a][j]-8*vx[x_b][j]+8*vx[x_c][j]-vx[x_d][j])+(((-1)*vy[i][j])/(12*delta_y))*(vx[i][y_e]-8*vx[i][y_f]+8*vx[i][y_g]-vx[i][y_h])+((-1)*g/(12*delta_x))*(h[x_a][j]-8*h[x_b][j]+8*h[x_c][j]-h[x_d][j]);
			ky[i][j]=(((-1)*vx[i][j])/(12*delta_x))*(vy[x_a][j]-8*vy[x_b][j]+8*vy[x_c][j]-vy[x_d][j])+(((-1)*vy[i][j])/(12*delta_y))*(vy[i][y_e]-8*vy[i][y_f]+8*vy[i][y_g]-vy[i][y_h])+((-1)*g/(12*delta_y))*(h[i][y_e]-8*h[i][y_f]+8*h[i][y_g]-h[i][y_h]);
			kh[i][j]=(-1/(12*delta_x))*(vx[x_a][j]*h[x_a][j]-8*vx[x_b][j]*h[x_b][j]+8*vx[x_c][j]*h[x_c][j]-vx[x_d][j]*h[x_d][j])+(-1/(12*delta_y))*(vy[i][y_e]*h[i][y_e]-8*vy[i][y_f]*h[i][y_f]+8*vy[i][y_g]*h[i][y_g]-vy[i][y_h]*h[i][y_h]);
			
			
			
		}
	}
	
	
	return;
}
void  matrixmult(double** k, double** m,  double c)
{
	for( int i=0;i<N_x;i++)
	{
		for (int j=0; j<N_y;j++)
		{
			k[i][j]=m[i][j]*c;
		}
	}
	
	
	return  ;
}
void matrixadd(double** k,double** h, double** m)
{
	for( int i=0;i<N_x;i++)
	{
		for (int j=0; j<N_y;j++)
		{
			k[i][j]=h[i][j]+m[i][j];
		}
	}
	return;
}



void writeData(int k, double** h, double x_min, double y_min, double delta_x, double delta_y, int N_x, int N_y)
{
	double x, y;
	char fileName[128];
	fstream file;
	sprintf(fileName, "waterequation_%05d.csv", k);
	file.open(fileName, ios::out);
	for(int i=0; i<N_x; i++)
	{
		x = x_min + i*delta_x;
		for(int j=0; j<N_y; j++)
		{
			y = y_min + j*delta_y;
			file << x << ",\t" << y << ",\t" << h[i][j] << endl;
		}
	}
	file.close();
	return;
}
void final_matrixadd(double** a, double** b, double** c, double** d, double** e, double f)
{
	for( int i=0; i<N_x;i++)
	{
		for (int j=0; j<N_y;j++)
		{
			a[i][j]=f*(b[i][j]/6+c[i][j]/3+d[i][j]/3+e[i][j]/6);
		}
	}
	return;
}
