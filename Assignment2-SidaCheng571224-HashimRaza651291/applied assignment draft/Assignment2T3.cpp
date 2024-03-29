/////////////////////////////////////////////////////////////////////////////
//
// Applied Numerical Methods
// with openMP and MPI
// 
//
// Problem:		Heat Equation \rho c \frac{\partial T}{\partial t} = (\nabla)^2 T
// Based on:    Example 15.3 
// Original Author: Dr. Steve Moore
// Re-write by :    Sida Cheng 571224
//                  Hashim Raza 651291
////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <math.h>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <mpi.h>
#include <omp.h>

using namespace std;

// Class definitions
class SparseMatrix {
public:
	SparseMatrix(int nrow, int nnzperrow) {
		// This constructor is called if we happen to know the number of rows
		// and an estimate of the number of nonzero entries per row.
		this->initialize(nrow, nnzperrow);
	}
	SparseMatrix() {
		// This constructor is called if we have no useful information
		N_row_ = 0;
		N_nz_ = 0;
		N_nz_rowmax_ = 0;
		allocSize_ = 0;
		val_ = NULL;
		col_ = NULL;
		row_ = NULL;
		nnzs_ = NULL;
	}
	~SparseMatrix() {
		if (val_)
			delete[] val_;
		if (col_)
			delete[] col_;
		if (row_)
			delete[] row_;
		if (nnzs_)
			delete[] nnzs_;
	}
	void initialize(int nrow, int nnzperrow) {
		N_row_ = nrow;
		N_nz_ = 0;
		N_nz_rowmax_ = nnzperrow;
		allocSize_ = N_row_ * N_nz_rowmax_;
		val_ = new double[allocSize_];
		col_ = new int[allocSize_];
		row_ = new int[N_row_ + 1];
		nnzs_ = new int[N_row_ + 1];

		memset(val_, 0, allocSize_ * sizeof(double));
		memset(col_, -1, allocSize_ * sizeof(int));
		memset(row_, 0, (N_row_ + 1) * sizeof(int));
		memset(nnzs_, 0, (N_row_ + 1) * sizeof(int));

		for (int k = 0, kk = 0; k < N_row_; k++, kk += N_nz_rowmax_) {
			row_[k] = kk;
		}
		return;
	}
	void finalize() {
		int minCol = 0;
		int insertPos = 0;
		int index = 0;

		// Now that the matrix is assembled we can set N_nz_rowmax_ explicitly by
		// taking the largest value in the nnzs_ array
		N_nz_rowmax_ = 0;
		for (int m = 0; m < N_row_; m++) {
			N_nz_rowmax_ = max(N_nz_rowmax_, nnzs_[m]);
		}

		double* tempVal = new double[N_nz_];
		int* tempCol = new int[N_nz_];
		int* tempRow = new int[N_row_ + 1];
		bool* isSorted = new bool[allocSize_]; // This array will help us sort the column indices

		memset(tempVal, 0, N_nz_ * sizeof(double));
		memset(tempCol, 0, N_nz_ * sizeof(int));
		memset(tempRow, 0, (N_row_ + 1) * sizeof(int));
		memset(isSorted, 0, allocSize_ * sizeof(bool));

		for (int m = 0; m < N_row_; m++) {
			for (int k = row_[m]; k < (row_[m] + nnzs_[m]); k++) {
				minCol = N_row_ + 1;
				for (int kk = row_[m]; kk < (row_[m] + nnzs_[m]); kk++) {
					if (!isSorted[kk] && col_[kk] < minCol) {
						index = kk;
						minCol = col_[index];
					}
				}
				tempVal[insertPos] = val_[index];
				tempCol[insertPos] = col_[index];
				isSorted[index] = true;
				insertPos++;
			}
			tempRow[m + 1] = tempRow[m] + nnzs_[m];
		}

		delete[] val_;
		delete[] col_;
		delete[] row_;
		delete[] nnzs_;
		delete[] isSorted;

		val_ = tempVal;
		col_ = tempCol;
		row_ = tempRow;
		nnzs_ = NULL;
		allocSize_ = N_nz_;

		return;
	}
	inline
	double& operator()(int m, int n) {
		int k = row_[m];
		bool foundEntry = false;

		// If the arrays are already full and inserting this entry would cause us to run off the end,
		// then we'll need to resize the arrays before inserting it
		if (nnzs_[m] >= N_nz_rowmax_) {
			this->reallocate();
		}
		// Search between row(m) and row(m+1) for col(k) = n (i.e. is the entry already in the matrix)
		while (k < (row_[m] + nnzs_[m]) && !foundEntry) {
			if (col_[k] == n) {
				foundEntry = true;
			}
			k++;
		}
		// If the entry is already in the matrix, then return a reference to it
		if (foundEntry) {
			return val_[k - 1];
		}
		// If the entry is not already in the matrix then we'll need to insert it
		else {
			N_nz_++;
			nnzs_[m]++;
			col_[k] = n;
			return val_[k];
		}
	}
	inline
	double& operator()(int k) {
		return val_[k];
	}
	void operator=(const SparseMatrix& A) {
		if (val_)
			delete[] val_;
		if (col_)
			delete[] col_;
		if (row_)
			delete[] row_;
		if (nnzs_)
			delete[] nnzs_;

		N_row_ = A.N_row_;
		N_nz_ = A.N_nz_;
		N_nz_rowmax_ = A.N_nz_rowmax_;
		allocSize_ = A.allocSize_;
		val_ = new double[allocSize_];
		col_ = new int[allocSize_];
		row_ = new int[N_row_ + 1];

		memcpy(val_, A.val_, N_nz_ * sizeof(double));
		memcpy(col_, A.col_, N_nz_ * sizeof(int));
		memcpy(row_, A.row_, (N_row_ + 1) * sizeof(int));
	}
	inline
	void multiply(double* u, double* v) {
		// Note: This function will perform a matrix vector multiplication with the input vector v, returning the output in u.
		for (int m = 0; m < N_row_; m++) {
			u[m] = 0.0;
			for (int k = row_[m]; k < row_[m + 1]; k++) {
				u[m] += val_[k] * v[col_[k]];
			}
		}
		return;
	}
	inline
	void multiply(double* u, double* v, bool* includerows, bool* includecols) {
		// Note: This function will perform a matrix vector multiplication on part of the matrix
		for (int m = 0; m < N_row_; m++) {
			u[m] = 0.0;
			if (includerows[m]) {
				for (int k = row_[m]; k < row_[m + 1]; k++) {

					if (includecols[col_[k]]) {
						u[m] += val_[k] * v[col_[k]];
					}
				}
			}
		}
		return;
	}
	inline
	void subtract(double u, SparseMatrix& A) {
		for (int k = 0; k < N_nz_; k++) {
			val_[k] -= (u * A.val_[k]);
		}
		return;
	}
	inline
	int getNnz() {
		return N_nz_;
	}
	inline
	int getNrow() {
		return N_row_;
	}
	void print(const char* name) {
		fstream matrix;
		cout << "Matrix " << name << " has " << N_row_ << " rows with " << N_nz_
				<< " non-zero entries - " << allocSize_ << " allocated."
				<< flush;
		matrix.open(name, ios::out);
		matrix << "Mat = [" << endl;
		for (int m = 0; m < N_row_; m++) {
			for (int n = row_[m]; n < row_[m + 1]; n++) {
				matrix << m + 1 << "\t" << col_[n] + 1 << "\t" << val_[n]
						<< endl;
			}
		}
		matrix << "];" << endl;
		matrix.close();
		cout << " Done." << flush << endl;
		return;
	}
protected:
	void reallocate() {
		// Double the memory allocation size
		N_nz_rowmax_ *= 2;

		allocSize_ = N_nz_rowmax_ * N_row_;

		// Create some temporary arrays of the new size
		double* tempVal = new double[allocSize_];
		int* tempCol = new int[allocSize_];

		memset(tempVal, 0, allocSize_ * sizeof(double));
		memset(tempCol, 0, allocSize_ * sizeof(int));

		for (int m = 0, mm = 0; m < N_row_; m++, mm += N_nz_rowmax_) {
			memcpy(&tempVal[mm], &val_[row_[m]], nnzs_[m] * sizeof(double));
			memcpy(&tempCol[mm], &col_[row_[m]], nnzs_[m] * sizeof(int));
			row_[m] = mm;
		}

		// Delete the memory used by the old arrays
		delete[] val_;
		delete[] col_;

		// Assign the addresses of the new arrays
		val_ = tempVal;
		col_ = tempCol;

		return;
	}
private:
	double* val_;		// [N_nz]    Stores the nonzero elements of the matrix
	int* col_;// [N_nz]    Stores the column indices of the elements in each row
	int* row_;		// [N_row+1] Stores the locations in val that start a row
	int* nnzs_;	// [N_row+1] Stores the number of nonzero entries per row during the assembly process
	int N_row_;			// The number of rows in the matric
	int N_nz_;	// The number of non-zero entries currently stored in the matrix
	int N_nz_rowmax_;// The maximum number of non-zero entries per row. This will be an estimate until the matrix is assembled
	int allocSize_;	// The number of non-zero entries currently allocated for in val_ and col_
};

class Boundary {
public:
	Boundary() {

	}
	string name_;
	string type_;
	int N_;
	int* indices_;
	double value_;
};

// Global variables
const double x_min = 0.00;
const double x_max = 0.02;
const double y_min = 0.00;
const double y_max = 0.02;
const double z_min = 0.00;
const double z_max = 0.02;
const double t_min = 0.00;
const double t_max = 100.00;
const double Delta_t = 0.1;
const double k_cond = 386.00;
const double rho = 8954.00;
const double C = 380.00;
const double T_air = 300.00;
const double h = 100.00;
const double Q_cpu = 40000.00;
const int N_t = static_cast<int>((t_max - t_min) / Delta_t + 1);
double* buffer = NULL;
int bufferSize = 0;

// Function declarations
void readData(char* filename, double**& Points, int**& Faces, int**& Elements, Boundary*& Boundaries, int& myN_p, int& myN_f, int& myN_e, int& myN_b, bool*& yourPoints, int myID);
void writeData(double* T, double**& Points, int**& Elements, int& myN_p,int& myN_e, int l, int myID);
void exchangeData(double* v, Boundary* Boundaries, int myN_b);
void assembleSystem(SparseMatrix& M, SparseMatrix& K, double* s, double* T, double** Points, int** Faces, int** Elements, Boundary* Boundaries, int myN_p, int myN_f, int myN_e, int myN_b, int myID);
void solve(SparseMatrix& A, double* T, double* b, Boundary* Boundaries, bool* yourPoints, int myN_b, int myID);
double computeInnerProduct(double* v1, double* v2, bool* yourPoints, int N_row);


int main(int argc, char** argv) {
	// Simulation parameters
	// Simulation parameters
    double**		Points		= NULL;
    int**			Faces		= NULL;
    int**			Elements	= NULL;
    Boundary*		Boundaries	= NULL;
    bool*			yourPoints	= NULL;
	double* 		buffer 		= NULL;
    int				myN_p		= 0;
    int				myN_f		= 0;
    int				myN_e		= 0;
    int				myN_b		= 0;
	int				myID		= 0;
    int				N_Processes	= 0;
	int				write		= 0;
	double			t			= 0;
	double wtime;
	int l;
	int m;
	// The following variable sets the number of threads for the program
	int N_Threads = omp_get_max_threads(); 	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myID);
	MPI_Comm_size(MPI_COMM_WORLD, &N_Processes);

	if (argc < 2) {
		if (myID == 0) {
			cerr << "No grid file specified" << endl;
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	} else {
		readData(argv[1], Points, Faces, Elements, Boundaries, myN_p, myN_f, myN_e, myN_b, yourPoints, myID);
	}

	// Allocate arrays
	double* T         = new double[myN_p];
	double* s 		  = new double[myN_p];
	double* b         = new double[myN_p];
	double* Z		  = new double[myN_p]; // Vector that's a multiple of Delta_T and s
	
	SparseMatrix M;
	SparseMatrix K;
	SparseMatrix A;

	if (myID == 0) {
		wtime = MPI_Wtime();
	}

	// Set initial condition
	for (m = 0; m < myN_p; m++) {
		T[m] = 300.0; // Which is T.air
	}
	
	assembleSystem(M, K, s, T, Points, Faces, Elements, Boundaries, myN_p, myN_f, myN_e, myN_b, myID);
	
	A = M;
	A.subtract(Delta_t, K); // At this point we have A = M-Delta_t*K

	// Compute the column vector to subtract from the right hand side to take account of fixed nodes
	// Use of OpenMP here is allowed 
	for (int m = 0; m < myN_p; m++) {
		Z[m] = Delta_t*s[m];
	}

	exchangeData(Z, Boundaries, myN_b);

	// Time marching loop
	for (l = 1; l < N_t; l++) {

		// Assemble b
		M.multiply(b, T);						// b = M*phi^l
		exchangeData(b, Boundaries, myN_b);
		for (int m = 0; m < myN_p; m++) {
			b[m] += Z[m];
		}	
		
		// Solve the linear system
		solve(A, T, b, Boundaries, yourPoints, myN_b, myID);

		// Write the solution
		if (!(l % 100)) {
			writeData(T, Points, Elements, myN_p, myN_e, l / 10, myID);
		}
	}

	if (myID == 0) {
		wtime = MPI_Wtime() - wtime;// Record the end time and calculate elapsed time
		cout << "Simulation took " << wtime << " seconds with " << N_Processes
				<< " processes " << " and " << N_Threads << " Threads per process..." << endl;
	}

	MPI_Buffer_detach(&buffer, &bufferSize);

	// Deallocate arrays - possible OpenMP
	for (int boundary = 0; boundary < myN_b; boundary++) {
		  delete[] Boundaries[boundary].indices_;
	}
	delete[] Points[0];
	delete[] Points;
	delete[] Faces[0];
	delete[] Faces;
	delete[] Elements[0];
	delete[] Elements;
	delete[] Boundaries;
	delete[] T;
	delete[] s;
	delete[] b;
	delete[] Z;
	delete[] buffer;

	MPI_Finalize();

	return 0;
}

void readData(char* filename, double**& Points, int**& Faces, int**& Elements, Boundary*& Boundaries, int& myN_p, int& myN_f, int& myN_e, int& myN_b, bool*& yourPoints, int myID) {
	
	int			myMaxN_sp	= 0;
	int			myMaxN_sb	= 0;
	int			maxN_sp		= 0;
	int			yourID		= 0;
	char		myFileName  [64];
    fstream		file;
	string      temp;

	
	
	if (myID == 0) {
		cout << "Reading " << filename << "'s... " << flush;
	}

	sprintf(myFileName, "%s%d", filename, myID);

	file.open(myFileName);
	if (!file.is_open()) {
		cerr << "Error opening file" << endl;
		exit(1);
	}

	file >> temp >> myN_p;
	file >> temp >> myN_f;
	file >> temp >> myN_e;
	file >> temp >> myN_b;

	Points			= new double*	[myN_p];
	Faces			= new int*		[myN_f];
	Elements		= new int*		[myN_e];
	Boundaries		= new Boundary  [myN_b];
    Points[0]       = new double    [myN_p*3];
    Faces[0]        = new int       [myN_f*3];
    Elements[0]     = new int       [myN_e*4];
	yourPoints		= new bool		[myN_p];
	
	
	for (int p = 1, pp = 3; p < myN_p; p++, pp += 3) {
		Points[p] = &Points[0][pp];
	}
	for (int f = 1, ff = 3; f < myN_f; f++, ff += 3) {
		Faces[f] = &Faces[0][ff];
	}
	for (int e = 1, ee = 4; e < myN_e; e++, ee += 4) {
		Elements[e] = &Elements[0][ee];
	}
	memset(yourPoints, false, myN_p * sizeof(bool));

	// Read in the points
	file >> temp;
	for (int p = 0; p < myN_p; p++) {
		file >> Points[p][0] >> Points[p][1] >> Points[p][2];
	}

	// Read in the faces
	file >> temp;
	for (int f = 0; f < myN_f; f++) {
		file >> Faces[f][0] >> Faces[f][1] >> Faces[f][2];
	}

	// Read in the elements
	file >> temp;
	for (int e = 0; e < myN_e; e++) {
		file >> Elements[e][0] >> Elements[e][1] >> Elements[e][2] >> Elements[e][3];
	}

	file >> temp;
	for (int b = 0; b < myN_b; b++) {
		file >> Boundaries[b].name_ >> Boundaries[b].type_ >> Boundaries[b].N_;
		Boundaries[b].indices_ = new int[Boundaries[b].N_];
		for (int n = 0; n < Boundaries[b].N_; n++) {
			file >> Boundaries[b].indices_[n];
		}
		file >> Boundaries[b].value_;
		if (Boundaries[b].type_ == "interprocess") {
			myMaxN_sb++;
			myMaxN_sp = max(myMaxN_sp, Boundaries[b].N_);
			yourID = static_cast<int>(Boundaries[b].value_);
			if (yourID > myID) {
				for (int p = 0; p < Boundaries[b].N_; p++) {
					yourPoints[Boundaries[b].indices_[p]] = true;
				}
			}
		}
	}

	MPI_Allreduce(&myMaxN_sp, &maxN_sp, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	buffer = new double[maxN_sp];
	bufferSize = (maxN_sp * sizeof(double) + MPI_BSEND_OVERHEAD) * myMaxN_sb;
	MPI_Buffer_attach(new char[bufferSize], bufferSize);

	file.close();

	if (myID == 0) {
		cout << "Done.\n" << flush;
	}
	return;
}

void writeData(double* T, double**& Points, int**& Elements, int& myN_p, int& myN_e, int l, int myID) {
	fstream file;
	char fileName[64];

	sprintf(fileName, "VTKBox/CPU_%02d_%04d.vtk", myID, l);

	file.open(fileName, ios::out);

	file << "# vtk DataFile Version 2.0" << endl;
	file << "Assignment 2, Created by Hash and Stan" << endl;
	file << "ASCII" << endl;
	file << "DATASET UNSTRUCTURED_GRID" << endl;

	file << "POINTS " << myN_p << " double" << endl;
	for (int p = 0; p < myN_p; p++) {
		file << setw(6) << setprecision(5) << fixed << Points[p][0] << "\t" << Points[p][1] << "\t" << Points[p][2] << endl;
	}

	file << "CELLS " << myN_e << " " << 5 * myN_e << endl;
	for (int e = 0; e < myN_e; e++) {
		file << "4\t" << Elements[e][0] << "\t" << Elements[e][1] << "\t" << Elements[e][2] << "\t" << Elements[e][3] << endl;
	}

	file << "CELL_TYPES " << myN_e << endl;
	for (int e = 0; e < myN_e; e++) {
		file << "10" << endl;
	}

	file << "POINT_DATA " << myN_p << endl;
	file << "SCALARS T double 1" << endl;
	file << "LOOKUP_TABLE default" << endl;
	// Set tempuerature values of the system
	for (int p = 0; p < myN_p; p++) {
		file << setprecision(5) << T[p] << endl;
	}

	file.close();

	return;
}

void exchangeData(double* v, Boundary* Boundaries, int myN_b) {
	int yourID = 0;
	int tag = 0;
	MPI_Status status;

	for (int b = 0; b < myN_b; b++) {
		if (Boundaries[b].type_ == "interprocess") {
			for (int p = 0; p < Boundaries[b].N_; p++) {
				buffer[p] = v[Boundaries[b].indices_[p]];
			}
			yourID = static_cast<int>(Boundaries[b].value_);
			MPI_Bsend(buffer, Boundaries[b].N_, MPI_DOUBLE, yourID, tag,
					MPI_COMM_WORLD);
		}
	}
	for (int b = 0; b < myN_b; b++) {
		if (Boundaries[b].type_ == "interprocess") {
			yourID = static_cast<int>(Boundaries[b].value_);
			MPI_Recv(buffer, Boundaries[b].N_, MPI_DOUBLE, yourID, tag,
					MPI_COMM_WORLD, &status);
			for (int p = 0; p < Boundaries[b].N_; p++) {
				v[Boundaries[b].indices_[p]] += buffer[p];
			}
		}
	}

	return;
}

void assembleSystem(SparseMatrix& M, SparseMatrix& K, double* s, double* T, double** Points, int** Faces, int** Elements, Boundary* Boundaries, int myN_p, int myN_f, int myN_e, int myN_b, int myID) {


	double x[4];
	double y[4];
	double z[4];
	double G[3][4];
	double Gp[3] = { 0.0, 0.0, 0.0 };
	double Gq[3] = { 0.0, 0.0, 0.0 };
	double M_e[4][4] = { { 2.0, 1.0, 1.0, 1.0 }, { 1.0, 2.0, 1.0, 1.0 }, { 1.0, 1.0, 2.0, 1.0 }, { 1.0, 1.0, 1.0, 2.0 } };
	double K_f[3][3] = { { 2.0, 1.0, 1.0 }, { 1.0, 2.0, 1.0 }, { 1.0, 1.0, 2.0 } };
	int Nodes[4] = { 0, 0, 0, 0 };
	double* Omega = new double[myN_e];
	double* Gamma = new double[myN_f];
	int m, n, o, p, q, r, b, f;

	// Assign all the indices to be free initially
#pragma omp parallel default(shared) private(m)
{	
	for (m = 0; m < myN_p; m++) {
		s[m] = 0.0;
	}

	// Calculate face lengths
	for (n = 0; n < myN_f; n++) {
		for (o = 0; o < 3; o++) {
			x[o] = Points[Faces[n][o]][0];
			y[o] = Points[Faces[n][o]][1];
			z[o] = Points[Faces[n][o]][2];
		}
			Gamma[n]	= sqrt(	  pow( (y[1]-y[0])*(z[2]-z[0]) - (z[1]-z[0])*(y[2]-y[0]) ,2.0) 
        	+ pow( (z[1]-z[0])*(x[2]-x[0]) - (x[1]-x[0])*(z[2]-z[0]) ,2.0)
       		+ pow( (x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0]) ,2.0) )/2.0;
	}

	// Calculate element areas
	for (q = 0; q < myN_e; q++) {
		for (p = 0; p < 4; p++) {
			x[p] = Points[Elements[q][p]][0];
			y[p] = Points[Elements[q][p]][1];
			z[p] = Points[Elements[q][p]][2];
		}
		Omega[q] = fabs(  x[0]*y[1]*z[2] - x[0]*y[2]*z[1] - x[1]*y[0]*z[2] 
                                         + x[1]*y[2]*z[0] + x[2]*y[0]*z[1] - x[2]*y[1]*z[0] 
                                         - x[0]*y[1]*z[3] + x[0]*y[3]*z[1] + x[1]*y[0]*z[3] 
                                         - x[1]*y[3]*z[0] - x[3]*y[0]*z[1] + x[3]*y[1]*z[0] 
                                         + x[0]*y[2]*z[3] - x[0]*y[3]*z[2] - x[2]*y[0]*z[3] 
                                         + x[2]*y[3]*z[0] + x[3]*y[0]*z[2] - x[3]*y[2]*z[0] 
                                         - x[1]*y[2]*z[3] + x[1]*y[3]*z[2] + x[2]*y[1]*z[3] 
                                         - x[2]*y[3]*z[1] - x[3]*y[1]*z[2] + x[3]*y[2]*z[1] ) /6; 
						
	}				 
	
}
	#pragma omp single
	// Assemble M, K, and s
	M.initialize(myN_p, 100);
	K.initialize(myN_p, 100);
	
	for (r = 0; r < myN_e; r++) {
		for (p = 0; p < 4; p++) {
			Nodes[p] = Elements[r][p];
			x[p] = Points[Nodes[p]][0];
			y[p] = Points[Nodes[p]][1];
			z[p] = Points[Nodes[p]][2];
		}

		G[0][0]	= (y[3]-y[1])*(z[2]-z[1])-(y[2]-y[1])*(z[3]-z[1]);	
		G[0][1]	= (y[2]-y[0])*(z[3]-z[2])-(y[2]-y[3])*(z[0]-z[2]);
		G[0][2]	= (y[1]-y[3])*(z[0]-z[3])-(y[0]-y[3])*(z[1]-z[3]);
		G[0][3] = (y[0]-y[2])*(z[1]-z[0])-(y[0]-y[1])*(z[2]-z[0]);
		G[1][0] = (x[2]-x[1])*(z[3]-z[1])-(x[3]-x[1])*(z[2]-z[1]);	
		G[1][1]	= (x[3]-x[2])*(z[2]-z[0])-(x[0]-x[2])*(z[2]-z[3]);
		G[1][2] = (x[0]-x[3])*(z[1]-z[3])-(x[1]-x[3])*(z[0]-z[3]);
		G[1][3] = (x[1]-x[0])*(z[0]-z[2])-(x[2]-x[0])*(z[0]-z[1]);
		G[2][0] = (x[3]-x[1])*(y[2]-y[1])-(x[2]-x[1])*(y[3]-y[1]);
		G[2][1] = (x[2]-x[0])*(y[3]-y[2])-(x[2]-x[3])*(y[0]-y[2]);
		G[2][2] = (x[1]-x[3])*(y[0]-y[3])-(x[0]-x[3])*(y[1]-y[3]);
		G[2][3] = (x[0]-x[2])*(y[1]-y[0])-(x[0]-x[1])*(y[2]-y[0]);

		// Outer loop over each node
		for (p = 0; p < 4; p++) {
			m = Nodes[p];
			Gp[0] = G[0][p];
			Gp[1] = G[1][p];
			Gp[2] = G[2][p];

			// Inner loop over each node
			for (q = 0; q < 4; q++) {
				n = Nodes[q];
				Gq[0] = G[0][q];
				Gq[1] = G[1][q];

				Gq[0] = G[0][q];
				Gq[1] = G[1][q];
				Gq[2] = G[2][q];
				M(m, n) += rho * C * M_e[p][q] * Omega[r] / 20;
				K(m, n) -= k_cond
						* (Gp[0] * Gq[0] + Gp[1] * Gq[1] + Gp[2] * Gq[2])
						/ (36 * Omega[r]);
			}
			s[m] += 0;
		}
	}

	// Apply boundary conditions
	for (b = 0; b < myN_b; b++) {
		if (Boundaries[b].type_ == "neumann") {
			for (f = 0; f < Boundaries[b].N_; f++) {
				for (p = 0; p < 3; p++) {
					Nodes[p] = Faces[Boundaries[b].indices_[f]][p];
					m = Nodes[p];
					s[m] += Q_cpu * Gamma[Boundaries[b].indices_[f]] / 3;
				}
			}
		} else if (Boundaries[b].type_ == "robin") {
			for (f = 0; f < Boundaries[b].N_; f++) {
				for (p = 0; p < 3; p++) {
					Nodes[p] = Faces[Boundaries[b].indices_[f]][p];
					m = Nodes[p];
					s[m] += h * T_air * Gamma[Boundaries[b].indices_[f]] / 3;
					for (int q = 0; q < 3; q++) {
						Nodes[q] = Faces[Boundaries[b].indices_[f]][q];
						n = Nodes[q];
						K(m, n) -= h * Gamma[Boundaries[b].indices_[f]]
								* K_f[p][q] / 12;
					}
				}
			}
		}
	}

	K.finalize();
	M.finalize();

	delete[] Gamma;
	delete[] Omega;

	MPI_Barrier (MPI_COMM_WORLD);

	return;
}

void solve(SparseMatrix& A, double* T, double* b, Boundary* Boundaries, bool* yourPoints, int myN_b, int myID) {
	
	int		N_row			= A.getNrow();
	double*	r_old			= new double [N_row];
	double*	r				= new double [N_row];
	double*	d				= new double [N_row];
	double*	Ad				= new double [N_row];
	double*	AT	      		= new double [N_row];
	double	alpha			= 0.0;
	double	beta			= 0.0;
	double	r_norm			= 0.0;
	double	tolerance		= 1e-8;
	double	maxIterations	= 1e+3;
	double	r_oldTr_old		= 0.0;
	double	rTr				= 0.0;
	double	dTAd			= 0.0;
	int		k				= 0;
	int		m				= 0;
	int		n				= 0;
	double tempT;

	memset(r_old, 0, N_row * sizeof(double));
	memset(r, 0, N_row * sizeof(double));
	memset(d, 0, N_row * sizeof(double));
	memset(Ad, 0, N_row * sizeof(double));

	// Compute the initial residual
	A.multiply(AT, T);
	exchangeData(AT, Boundaries, myN_b);
	#pragma omp parallel default(shared) private(m)
	{
		for (m = 0; m < N_row; m++) {
			r_old[m] = b[m] - AT[m];
			d[m] = r_old[m];
		}
		#pragma omp single
		{
			r_oldTr_old = computeInnerProduct(r_old, r_old, yourPoints, N_row);
			r_norm = sqrt(r_oldTr_old);
		}
// Conjugate Gradient iterative loop
		while (r_norm > tolerance && k < maxIterations) {
			#pragma omp single
			{
				A.multiply(Ad, d);
				exchangeData(Ad, Boundaries, myN_b);
				dTAd = computeInnerProduct(d, Ad, yourPoints, N_row);
			}
			alpha = r_oldTr_old / dTAd;
			#pragma omp for schedule(static) reduction(+ : tempT)
			for (m = 0; m < N_row; m++) {
				tempT = T[m];
				tempT += alpha * d[m];
				T[m] = tempT;
			}

			#pragma omp for schedule(static)
			for (m = 0; m < N_row; m++) {
				r[m] = r_old[m] - alpha * Ad[m];
			}

			#pragma omp single
			{
				rTr = computeInnerProduct(r, r, yourPoints, N_row);
			}

			beta = rTr / r_oldTr_old;
			#pragma omp for schedule(static)
			for (m = 0; m < N_row; m++) {
				d[m] = r[m] + beta * d[m];
			}
			#pragma omp for schedule(static)
			for (m = 0; m < N_row; m++) {
				r_old[m] = r[m];
			}

			#pragma omp single
			{
				r_oldTr_old = rTr;
				r_norm = sqrt(rTr);
				k++;
			}
			#pragma omp barrier
		}
	}


	delete[] r_old;
	delete[] r;
	delete[] d;
	delete[] Ad;
	delete[] AT;

	return;
}

double computeInnerProduct(double* v1, double* v2, bool* yourPoints, int N_row) {
	double myInnerProduct = 0.0;
	double innerProduct = 0.0;

	for (int m = 0; m < N_row; m++) {
		if (!yourPoints[m]) {
			myInnerProduct += v1[m] * v2[m];
		}
	}

	MPI_Allreduce(&myInnerProduct, &innerProduct, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

	return innerProduct;
}

