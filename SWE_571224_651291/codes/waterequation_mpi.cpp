#include <string.h>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;
void    f(double** vx, double** vy, double** h, double** kx, double** ky, double** kh,double delta_x, double delta_y,int myN_x, int myN_y);
void  matrixmult(double** k, double** m,  double c,int myN_x, int myN_y);
void matrixadd(double** k,double** h, double** m,int myN_x, int myN_y);
void final_matrixadd(double** a, double** b, double** c, double** d, double** e, double f,int myN_x, int myN_y);
void  writeData(int k, double** h, double x_min, double y_min, double Delta_x, double Delta_y, int myN_x, int myN_y, int N_x, int N_y, int* myCoords, MPI_Comm& Comm2D);




// Global variables
int l;
const double    x_min       =  0.00;
const double    x_max       = 100.00;
const double    t_min       =  0.00;
const double    t_max       = 100.00;
const double    y_min           =0.00;
const double    y_max           =100.00;
const double    delta_t     =0.2;
const int       N_t         =  (t_max-t_min)/delta_t+1;
const double    g           =  9.81;
const int X=0;
const int Y=1;
const int numElementsPerBlock   = 1;

// Function declarations

int     main(int argc, char** argv)
{

    int myID;
    int N_Processers;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &N_Processers);
    MPI_Comm_rank(MPI_COMM_WORLD, &myID);
    
    // Simulation parameters
    int N=sqrt(N_Processers);
    int dimensions[2]={N,N};
    int isPeriodic[2]={1,1};
    int myCoords[2]={0,0};
    int reorder= 1; 
    int leftNeighbor = 0;
    int rightNeighbor = 0;
    int bottomNeighbor = 0;
    int topNeighbor = 0;
    double wtime=0;
    
    int mm;
    int m;
    int myN_x=20;
    int myN_y=20;
    int N_x=myN_x*dimensions[0];
    int N_y=myN_y*dimensions[1];
    double delta_x=(x_max-x_min)/(N_x-1);
    double delta_y=(y_max-y_min)/(N_y-1);

    MPI_Status status;
    MPI_Datatype strideType;
    MPI_Comm Comm2D;
    // Allocate arrays



    double** vx = new double* [myN_x+4];
    vx[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        vx[m]=&vx[0][mm];
    }
    double** vy = new double* [myN_x+4];
    vy[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        vy[m]=&vy[0][mm];
    }
    double** h = new double* [myN_x+4];
    h[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        h[m]=&h[0][mm];
    }
    double** vxt = new double* [myN_x+4];       //temp vx
    vxt[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        vxt[m]=&vxt[0][mm];
    }
    double** vyt = new double* [myN_x+4];       //temp vy
    vyt[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        vyt[m]=&vyt[0][mm];
    }
    double** ht = new double* [myN_x+4];        //temp h
    ht[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        ht[m]=&ht[0][mm];
    }
    double** kx1 = new double* [myN_x+4];        // kx1 ky1 kh1
    kx1[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kx1[m]=&kx1[0][mm];
    }
    double** ky1 = new double* [myN_x+4];
    ky1[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        ky1[m]=&ky1[0][mm];
    }
    double** kh1 = new double* [myN_x+4];
    kh1[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kh1[m]=&kh1[0][mm];
    }
    double** kx2 = new double* [myN_x+4];        // kx2 ky2 kh2
    kx2[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kx2[m]=&kx2[0][mm];
    }
    double** ky2 = new double* [myN_x+4];
    ky2[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1,mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        ky2[m]=&ky2[0][mm];
    }
    double** kh2 = new double* [myN_x+4];
    kh2[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kh2[m]=&kh2[0][mm];
    }
    double** kx3 = new double* [myN_x+4];         //kx3 ky3 kh3
    kx3[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kx3[m]=&kx3[0][mm];
    }
    double** ky3 = new double* [myN_x+4];
    ky3[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        ky3[m]=&ky3[0][mm];
    }
    double** kh3 = new double* [myN_x+4];
    kh3[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kh3[m]=&kh3[0][mm];
    }
    double** kx4 = new double* [myN_x+4];               //kx4 ky4 kh4
    kx4[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kx4[m]=&kx4[0][mm];
    }
    double** ky4 = new double* [myN_x+4];
    ky4[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        ky4[m]=&ky4[0][mm];
    }
    double** kh4 = new double* [myN_x+4];
    kh4[0]=new double [(myN_x+4)*(myN_y+4)];
    for(  m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kh4[m]=&kh4[0][mm];
    }
    double** kxt = new double* [myN_x+4];               //temp k values
    kxt[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kxt[m]=&kxt[0][mm];
    }
    double** kyt = new double* [myN_x+4];
    kyt[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kyt[m]=&kyt[0][mm];
    }
    double** kht = new double* [myN_x+4];
    kht[0]=new double [(myN_x+4)*(myN_y+4)];
    for( m=1, mm=myN_y+4;m<myN_x+4; m++,mm+=myN_y+4)
    {
        kht[m]=&kht[0][mm];
    }
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, isPeriodic, reorder, &Comm2D);
    MPI_Comm_rank(Comm2D, &myID);
    MPI_Cart_coords(Comm2D, myID, 2, myCoords);
    MPI_Cart_shift(Comm2D, X, 1, &leftNeighbor,   &rightNeighbor);
    MPI_Cart_shift(Comm2D, Y, 1, &bottomNeighbor, &topNeighbor);
    MPI_Barrier(Comm2D);




    if(myID==0)
    {
        wtime=MPI_Wtime();
    }

    // Set initial condition

    for(int i=2; i<myN_x+2; i++)
    {
        for(int j=2; j<myN_y+2;j++)
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
            double x=(myCoords[X]*myN_x + i - 2)*delta_x;
            double y=(myCoords[Y]*myN_y + i - 2)*delta_y;
            h[i][j]=1+0.5*exp((-1/0.04)*(pow(x-30,2)+pow(y-30,2)));
        }
    }


    MPI_Type_vector(myN_x, numElementsPerBlock, myN_y+4, MPI_DOUBLE, &strideType);            //DATATYPE for top and bottom neighbors
    MPI_Type_commit(&strideType);

    
    

    for(l=0; l<N_t-1; l++)
    {
       
        //for kx ky kh for the first values
        MPI_Sendrecv(&(vx[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vx[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vx[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vx[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vx[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vx[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vx[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vx[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vy[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vy[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vy[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vy[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vy[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vy[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vy[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vy[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(h[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(h[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(h[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(h[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(h[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(h[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(h[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(h[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vx[2][2]), 1, strideType,bottomNeighbor, 0, &(vx[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vx[2][3]), 1, strideType,bottomNeighbor, 0, &(vx[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vx[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(vx[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vx[2][myN_y]), 1, strideType,topNeighbor, 0, &(vx[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vy[2][2]), 1, strideType,bottomNeighbor, 0, &(vy[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vy[2][3]), 1, strideType,bottomNeighbor, 0, &(vy[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vy[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(vy[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vy[2][myN_y]), 1, strideType,topNeighbor, 0, &(vy[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(h[2][2]), 1, strideType,bottomNeighbor, 0, &(h[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(h[2][3]), 1, strideType,bottomNeighbor, 0, &(h[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(h[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(h[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(h[2][myN_y]), 1, strideType,topNeighbor, 0, &(h[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);

        f(vx,vy,h,kx1,ky1,kh1,delta_x,delta_y,myN_x,myN_y);
        // kx1 ky1 kh1 end here
        matrixmult(kxt,kx1,delta_t/2,myN_x,myN_y);
        matrixmult(kyt,ky1,delta_t/2,myN_x,myN_y);
        matrixmult(kht,kh1,delta_t/2,myN_x,myN_y);
        matrixadd(vxt,vx,kxt,myN_x,myN_y);
        matrixadd(vyt,vy,kyt,myN_x,myN_y);
        matrixadd(ht,h,kht,myN_x,myN_y);
        MPI_Sendrecv(&(vxt[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vxt[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vxt[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vxt[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vxt[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vyt[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vyt[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vyt[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vyt[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(ht[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(ht[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(ht[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(ht[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][2]), 1, strideType,bottomNeighbor, 0, &(vxt[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][3]), 1, strideType,bottomNeighbor, 0, &(vxt[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(vxt[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][myN_y]), 1, strideType,topNeighbor, 0, &(vxt[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][2]), 1, strideType,bottomNeighbor, 0, &(vyt[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][3]), 1, strideType,bottomNeighbor, 0, &(vyt[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(vyt[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][myN_y]), 1, strideType,topNeighbor, 0, &(vyt[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][2]), 1, strideType,bottomNeighbor, 0, &(ht[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][3]), 1, strideType,bottomNeighbor, 0, &(ht[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(ht[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][myN_y]), 1, strideType,topNeighbor, 0, &(ht[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);

        f(vxt,vyt,ht,kx2,ky2,kh2,delta_x,delta_y,myN_x,myN_y);
        //kx2 ky2 kh2 end here
        matrixmult(kxt,kx2,delta_t/2,myN_x,myN_y);
        matrixmult(kyt,ky2,delta_t/2,myN_x,myN_y);
        matrixmult(kht,kh2,delta_t/2,myN_x,myN_y);
        matrixadd(vxt,vx,kxt,myN_x,myN_y);
        matrixadd(vyt,vy,kyt,myN_x,myN_y);
        matrixadd(ht,h,kht,myN_x,myN_y);
        
        MPI_Sendrecv(&(vxt[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vxt[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vxt[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vxt[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vxt[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vyt[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vyt[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vyt[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vyt[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(ht[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(ht[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(ht[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(ht[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][2]), 1, strideType,bottomNeighbor, 0, &(vxt[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][3]), 1, strideType,bottomNeighbor, 0, &(vxt[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(vxt[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][myN_y]), 1, strideType,topNeighbor, 0, &(vxt[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][2]), 1, strideType,bottomNeighbor, 0, &(vyt[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][3]), 1, strideType,bottomNeighbor, 0, &(vyt[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(vyt[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][myN_y]), 1, strideType,topNeighbor, 0, &(vyt[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][2]), 1, strideType,bottomNeighbor, 0, &(ht[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][3]), 1, strideType,bottomNeighbor, 0, &(ht[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(ht[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][myN_y]), 1, strideType,topNeighbor, 0, &(ht[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        
        f(vxt,vyt,ht,kx3,ky3,kh3,delta_x,delta_y,myN_x,myN_y);
        //kx3 ky3 kh3 end here
        matrixmult(kxt,kx3,delta_t,myN_x,myN_y);
        matrixmult(kyt,ky3,delta_t,myN_x,myN_y);
        matrixmult(kht,kh3,delta_t,myN_x,myN_y);
        matrixadd(vxt,vx,kxt,myN_x,myN_y);
        matrixadd(vyt,vy,kyt,myN_x,myN_y);
        matrixadd(ht,h,kht,myN_x,myN_y);
        
        MPI_Sendrecv(&(vxt[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vxt[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vxt[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vxt[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vxt[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vyt[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(vyt[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vyt[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(vyt[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(ht[myN_x+2][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[3][2]), myN_y, MPI_DOUBLE,leftNeighbor, 0, &(ht[myN_x+3][2]), myN_y, MPI_DOUBLE, rightNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[myN_x+1][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(ht[1][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[myN_x][2]),myN_y, MPI_DOUBLE, rightNeighbor, 0, &(ht[0][2]), myN_y, MPI_DOUBLE, leftNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][2]), 1, strideType,bottomNeighbor, 0, &(vxt[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][3]), 1, strideType,bottomNeighbor, 0, &(vxt[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(vxt[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vxt[2][myN_y]), 1, strideType,topNeighbor, 0, &(vxt[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][2]), 1, strideType,bottomNeighbor, 0, &(vyt[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][3]), 1, strideType,bottomNeighbor, 0, &(vyt[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(vyt[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(vyt[2][myN_y]), 1, strideType,topNeighbor, 0, &(vyt[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][2]), 1, strideType,bottomNeighbor, 0, &(ht[2][myN_y+2]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][3]), 1, strideType,bottomNeighbor, 0, &(ht[2][myN_y+3]), 1, strideType, topNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][myN_y+1]), 1, strideType,topNeighbor, 0, &(ht[2][1]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);
        MPI_Sendrecv(&(ht[2][myN_y]), 1, strideType,topNeighbor, 0, &(ht[2][0]), 1, strideType, bottomNeighbor, 0, Comm2D, &status);

        
        f(vxt,vyt,ht,kx4,ky4,kh4,delta_x,delta_y,myN_x,myN_y);
        
        final_matrixadd(kxt,kx1,kx2,kx3,kx4,delta_t,myN_x,myN_y);
        final_matrixadd(kyt,ky1,ky2,ky3,ky4,delta_t,myN_x,myN_y);
        final_matrixadd(kht,kh1,kh2,kh3,kh4,delta_t,myN_x,myN_y);
        matrixadd(vx,vx,kxt,myN_x,myN_y);
        matrixadd(vy,vy,kyt,myN_x,myN_y);
        matrixadd(h,h,kht,myN_x,myN_y);
        
        
        
        
        
        
        
        
        if(l%50==0)
        {
            
        writeData(l,h,x_min,y_min,delta_x,delta_y,myN_x,myN_y,N_x,N_y,myCoords,Comm2D);
            
        }           
        
        
        
    }
    

    
    
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
    
    MPI_Finalize();
    return 0;
    
}

void    f(double** vx, double** vy, double** h, double** kx, double** ky, double** kh,double delta_x, double delta_y,int myN_x, int myN_y)
{
    
    
    for(int i=2; i<myN_x+2; i++)
    {
        for (int j=2;j<myN_y+2; j++)
        {

            kx[i][j]=(((-1)*vx[i][j])/(12*delta_x))*(vx[i-2][j]-8*vx[i-1][j]+8*vx[i+1][j]-vx[i+2][j])+(((-1)*vy[i][j])/(12*delta_y))*(vx[i][j-2]-8*vx[i][j-1]+8*vx[i][j+1]-vx[i][j+2])+((-1)*g/(12*delta_x))*(h[i-2][j]-8*h[i-1][j]+8*h[i+1][j]-h[i+2][j]);
            ky[i][j]=(((-1)*vx[i][j])/(12*delta_x))*(vy[i-2][j]-8*vy[i-1][j]+8*vy[i+1][j]-vy[i+2][j])+(((-1)*vy[i][j])/(12*delta_y))*(vy[i][j-2]-8*vy[i][j-1]+8*vy[i][j+1]-vy[i][j+2])+((-1)*g/(12*delta_y))*(h[i][j-2]-8*h[i][j-1]+8*h[i][j+1]-h[i][j+2]);
            kh[i][j]=(-1/(12*delta_x))*(vx[i-2][j]*h[i-2][j]-8*vx[i-1][j]*h[i-1][j]+8*vx[i+1][j]*h[i+1][j]-vx[i+2][j]*h[i+2][j])+(-1/(12*delta_y))*(vy[i][j-2]*h[i][j-2]-8*vy[i][j-1]*h[i][j-1]+8*vy[i][j+1]*h[i][j+1]-vy[i][j+2]*h[i][j+2]);
            
            
            
        }
    }
    
    
    return;
}
void  matrixmult(double** k, double** m,  double c,int myN_x, int myN_y)
{
    for( int i=2;i<myN_x+2;i++)
    {
        for (int j=2; j<myN_y+2;j++)
        {
            k[i][j]=m[i][j]*c;
        }
    }
    
    
    return  ;
}
void matrixadd(double** k,double** h, double** m,int myN_x, int myN_y)
{
    for( int i=2;i<myN_x+2;i++)
    {
        for (int j=2; j<myN_y+2;j++)
        {
            k[i][j]=h[i][j]+m[i][j];
        }
    }
    return;
}




void final_matrixadd(double** a, double** b, double** c, double** d, double** e, double f,int myN_x, int myN_y)
{
    for( int i=2; i<myN_x+2;i++)
    {
        for (int j=2; j<myN_y+2;j++)
        {
            a[i][j]=f*(b[i][j]/6+c[i][j]/3+d[i][j]/3+e[i][j]/6);
        }
    }
    return;
}
    
void writeData(int k, double** h, double x_min, double y_min, double Delta_x, double Delta_y, int myN_x, int myN_y, int N_x, int N_y, int* myCoords, MPI_Comm& Comm2D)
{
    double x, y;
    char fileName[128];
    MPI_File file;
    int offset;
    int N_CharPer_h = 24;
    int N_BytesPer_h = N_CharPer_h * sizeof(char);
    char buffer[N_BytesPer_h];
    sprintf(fileName, "waterequationMPI_%05d.csv", k);
    MPI_File_open(Comm2D, fileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    for(int i=2; i<myN_x+2; i++)
    {
        x = x_min + (myCoords[X]*myN_x + (i-2))*Delta_x;
        for(int j=2; j<myN_x+2; j++)
        {
            y = y_min + (myCoords[Y]*myN_y + (j-2))*Delta_y;
            offset = ( N_y*(myCoords[X]*myN_x + (i-1)) + myCoords[Y]*myN_y + (j-1) )*N_BytesPer_h;
            sprintf(buffer, "%01.4f,\t%01.4f,\t%+1.4f\n", x, y, h[i][j]);
            MPI_File_seek(file, offset, MPI_SEEK_SET);
            MPI_File_write(file, buffer, N_BytesPer_h, MPI_CHAR, MPI_STATUS_IGNORE);
        }
    }
    MPI_File_close(&file);
    return;
}
