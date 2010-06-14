#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include "AVLTree.h"
#include <iostream>
#include <vector>
using namespace std;

class __declspec( dllexport ) SparseMatrixElement
{
public:

	SparseMatrixElement(){};
	SparseMatrixElement( int row, int col, double val ){ m_row = row; m_col = col; m_val = val;};
	~SparseMatrixElement(){};

	bool operator<( const SparseMatrixElement & elem )
	{
		if( m_col < elem.m_col ) return true;
		if( m_col > elem.m_col ) return false;
		if( m_row < elem.m_row ) return true;
		return false;
	};

	bool operator==( const SparseMatrixElement & elem )
	{
		return m_row == elem.m_row && m_col == elem.m_col;
	};

	int    & row() { return m_row; };
	int	   & col() { return m_col; };
	double & val() { return m_val; };

protected:
	int    m_row;
	int    m_col;
	double m_val;
};

struct StructMatrix
{
public:

	StructMatrix() { buffer = NULL; row_ind = NULL; col_cnt = NULL; };
	~StructMatrix()
	{ 
		if( buffer ) delete []buffer; 
		if( row_ind) delete []row_ind;
		if( col_cnt) delete []col_cnt;
	}
public:

	double * buffer;
	int    * row_ind;
	int    * col_cnt;
	int      nnz;
	int      nr;
};

typedef StructMatrix SM;

class __declspec( dllexport ) SparseMatrix
{
public:
	SparseMatrix(int rows=0, int cols=0)  { m_rows = rows; m_cols = cols; }

	~SparseMatrix();

	bool insert( int row, int col, double val );
	double element( int row, int col);
	void   set( int row, int col, double val );
	void   update( int row, int col, double val );
	int    size();

	int    rows() { return m_rows; };
	int    cols() { return m_cols; };

	void set_size(int rows, int cols) { m_rows = rows; m_cols = cols; }

	void operator *= (double s );
	SparseMatrix * Transpose();
	SparseMatrix * operator*( SparseMatrix & m );

	AVLTree<SparseMatrixElement>& get_elements() { return m_elements; }

	bool save(char* fname);
	bool load(char* fname);

	void convert( SM & sm );
	void output();

	bool IsPositiveDefinite();

	void save_to_mtl(char* fname);
	void save_to_matlab(char* fname);

protected:
	AVLTree<SparseMatrixElement> m_elements;

	int m_rows;
	int m_cols;
};

inline SparseMatrix* Add(SparseMatrix* A, SparseMatrix* B)
{
	if (A->rows() != B->rows() || A->cols() != B->cols())
	{
		cout << "A and B do not match!" << endl;
		return NULL;
	}

	SparseMatrix* C = new SparseMatrix(A->rows(), A->cols());

	int rows = A->rows();
	int cols = A->cols();
	for (int i = 0; i < rows; i ++)
	{
		for (int j = 0; j < cols; j ++)
		{
			double a = A->element(i, j);
			double b = B->element(i, j);

			if (a != 0 || b != 0)
			{
				C->insert(i, j, a+b);
			}
		}
	}

	return C;
}

__declspec(dllexport) void SuperLUSolve( SparseMatrix * A, double * b, int nrhs ); // AX = b solver, requires A is a square matrix
__declspec(dllexport) void MosekSolve( SparseMatrix * A, SparseMatrix * b ); // AX = b solver using MOSEK, A must be a symetric, semi-positive definite matrix
__declspec(dllexport) void CGSolve(SparseMatrix* A, double* b, int nrhs);
__declspec(dllexport) void LSolve(SparseMatrix* A, double* b, int nrhs);


// ne : number of requested eigenvalues, 
// nconv: number of converged eigenvalues
//__declspec(dllexport) void LaplaceEigen(SparseMatrix* A, int ne, double* evalue, double* evec, int& nconv); 


class __declspec( dllexport ) MyMatrix
{
public:
	vector< vector<double>* > elements;
	vector< vector<int>* > indices;
	int nr; // # of rows
};

inline void mult( const MyMatrix &A, const double *v, double *w ) 
{
	int n = A.nr;
	for (int i = 0; i < n; i ++)
	{
		w[i] = 0;

		for (int j = 0; j < (int)A.indices[i]->size(); j ++)
		{
			double f = A.elements[i]->at(j);
			int k = A.indices[i]->at(j);

			w[i] += f*v[k];
		}
//		cout << "w[" << i << "]=" << w[i] << endl;
	}
}

__declspec(dllexport) void LSolve(MyMatrix& A, double* b, int nrhs);

#endif
