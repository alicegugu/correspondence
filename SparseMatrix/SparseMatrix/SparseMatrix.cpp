#include "SparseMatrix.h"
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <fstream>
#include <assert.h>
#include "matrix/eigval.c"
//#include "matrix/solv.c"
#include "nr.h"
#include "nrutil.h"
#include <vector>

using namespace std;

bool SparseMatrix::insert( int row, int col, double val )
{
	SparseMatrixElement elem(row,col,0);
	if( m_elements.find( &elem ) != NULL ) return false;

	SparseMatrixElement * element = new SparseMatrixElement( row, col, val );

	m_elements.insert( element );
	return true;
}

void SparseMatrix::set( int row, int col, double val )
{
	SparseMatrixElement elem(row,col,0);
	SparseMatrixElement * e = m_elements.find( &elem );
	if( e != NULL ) 
		e->val() = val;
	else
		insert( row, col, val );
}


void SparseMatrix::update( int row, int col, double val )
{
	SparseMatrixElement elem(row,col,0);
	SparseMatrixElement * e = m_elements.find( &elem );
	if( e != NULL ) 
		e->val() += val;
	else
		insert( row, col, val );
}

SparseMatrix::~SparseMatrix()
{
	for( AVLTreeIterator<SparseMatrixElement> iter( m_elements ); !iter.end(); ++ iter )
	{
		SparseMatrixElement * elem = *iter;
		delete elem;
	}
}

int SparseMatrix::size() 
{
	int count = 0;
	for( AVLTreeIterator<SparseMatrixElement> iter( m_elements ); !iter.end(); ++ iter )
	{
		count ++;
	}
	return count;
}

double SparseMatrix::element(int row, int col )
{
	SparseMatrixElement elem( row, col, 0 );
	SparseMatrixElement * e = m_elements.find( &elem );
	if( e == NULL ) return 0;

	return e->val();
}


void SuperLUSolve( SparseMatrix * A, double* b, int nrhs )
{
	if( A->cols() != A->rows() ) 
	{
		cout << "fatal error in SuperLUSolve! A->cols() = " << A->cols() << ", A->rows() = " << A->rows() << endl;
		return;
	}

	cout << "SuperLUSolve: nrhs = " << nrhs << endl;
	cout << "A size " << A->cols() << endl;
	SM sm;

	A->convert( sm );

	cout << "calling SuperLU..." << endl;
	SuperLU( sm.buffer, sm.row_ind, sm.col_cnt, sm.nr, sm.nnz, nrhs, b);
}

void SparseMatrix::operator *=( double s )
{
	for( AVLTreeIterator<SparseMatrixElement> eiter(m_elements); !eiter.end(); ++ eiter )
	{
		SparseMatrixElement * elem = *eiter;
		elem->val() *= s;
	}
}

SparseMatrix * SparseMatrix::Transpose()
{
	SparseMatrix * t = new SparseMatrix( m_cols, m_rows );

	for( AVLTreeIterator<SparseMatrixElement> iter( m_elements ); !iter.end(); ++ iter )
	{
		SparseMatrixElement * ele = *iter;
		t->insert( ele->col(), ele->row(), ele->val() );
	}
	return t;
}

SparseMatrix * SparseMatrix::operator*( SparseMatrix & m )
{
	if( m_cols != m.rows() ) return NULL;

	int rows = m_rows;
	int cols = m.cols();

	double  * buffer = new double[rows*cols];
	memset( buffer, 0, sizeof(double)*rows*cols);

	SparseMatrix * result = new SparseMatrix(rows,cols);

	SparseMatrix * mt = m.Transpose();

	SM sm;
	convert( sm );

	SM smt;
	mt->convert( smt );


	for( int c = 0 ; c < m_cols; c ++ )
	{
		for( int k = sm.col_cnt[c]; k < sm.col_cnt[c+1]; k ++ )
		{
			double v = sm.buffer[k];
			int    r = sm.row_ind[k];
			
			for( int kt = smt.col_cnt[c]; kt < smt.col_cnt[c+1]; kt ++ )
			{
				double vt = smt.buffer[kt];
				int    rt = smt.row_ind[kt];
				buffer[r * cols + rt] += v * vt;
			}
		}
	}

	for( int i = 0; i < rows; i ++ )
	for( int j = 0; j < cols; j ++ )
	{
		if( buffer[i*cols+j] != 0 )
		result->insert(i,j,buffer[i*cols+j]);
	}

	delete mt;

	return result;
}

void SparseMatrix::convert(  SM & sm )
{

	sm.nnz	    = size();
	sm.buffer   = new double[sm.nnz];
	sm.row_ind  = new int[sm.nnz];
	sm.col_cnt  = new int[cols()+1];
	sm.nr = rows();


	int k = 0;

	int cur_col = 0;

	sm.col_cnt[0] = 0;
	
	for( AVLTreeIterator<SparseMatrixElement> iter( m_elements ); !iter.end(); ++ iter )
	{
		SparseMatrixElement * element = * iter;

		int r = element->row();
		int c = element->col();	
		double v = element->val();

		sm.row_ind[k] = r;
		sm.buffer[k] = v;

		if( c != cur_col )
		{
			sm.col_cnt[c] = k;

			cur_col ++;

			for( ; cur_col < c; cur_col ++ )
			{
				sm.col_cnt[cur_col] = sm.col_cnt[c];
			}
		}

		k++;
	}

	for(cur_col++ ; cur_col < cols() +1; cur_col ++ )
	{
		sm.col_cnt[cur_col] = sm.nnz;
	}
}

void SparseMatrix::save_to_mtl(char* fname)
{
	ofstream ofs(fname, ios::out);
	int count = size();

	ofs << rows() << " " << cols() << endl;
	ofs << count << endl;

	for( AVLTreeIterator<SparseMatrixElement> eiter(m_elements); !eiter.end(); ++ eiter )
	{
		SparseMatrixElement * elem = *eiter;
		
		ofs << elem->row() << " " << elem->col() << " " << elem->val() << endl;
	}

	ofs.close();
}

void SparseMatrix::save_to_matlab(char* fname)
{
	ofstream ofs(fname, ios::out);
	ofs.precision(10);

	int count = size();

	for( AVLTreeIterator<SparseMatrixElement> eiter(m_elements); !eiter.end(); ++ eiter )
	{
		SparseMatrixElement * elem = *eiter;
		
		ofs << elem->row()+1 << "\t\t" << elem->col()+1 << "\t\t" << elem->val() << endl;
	}

	ofs.close();
}


void SparseMatrix::output()
{
	for( int i = 0; i < m_rows; i ++ )
	{
		for( int j = 0; j < m_cols; j ++ )
		{
			double v = element(i,j);
			//if( v!= 0 )
			//	fprintf(stderr, "(%d,%d) %f ", i,j, v );
			fprintf(stderr,"%.3f ",v);
		}
		fprintf(stderr,"\n");
	}
}

bool SparseMatrix::save(char* fname)
{
	ofstream ofs(fname, ios::out);
	ofs.precision(20);
	ofs << rows() << " " << cols() << endl;
	ofs << size() << endl;

	for (AVLTreeIterator<SparseMatrixElement> iter(m_elements); !iter.end(); ++ iter)
	{
		SparseMatrixElement * ele = *iter;
		ofs << ele->row() << " " << ele->col() << " " << ele->val() << endl;
	}

	ofs.close();

	return true;
}

bool SparseMatrix::IsPositiveDefinite()
{
	if (cols() != rows())
	{
		cout << "it is not a square matrix!" << endl;
		return false;
	}

	int n = cols();

	double* QQ = new double[n*n]; // test eigenvalues
	for (int qq = 0; qq < n*n; qq ++)
		QQ[qq] = 0;

	cout << "\nTest: is Q positive definite?" << endl;
	double* ev = new double[n];

	eigval(QQ, ev, n);
	
	bool positive_definite(true);
	for (int i = 0; i < n; i ++)
	{
		if (ev[i] < 0)
		{
			cout << "eigenvalues[" << i << "] = " << ev[i] << endl;
			positive_definite = false;
		}

		if (ev[i] < 0.001)
		{
			cout << "eigenvalues[" << i << "] is near zero!" << ev[i] << endl;
		}
	}

	delete[] ev;
	delete[] QQ;

	return positive_definite;
}

bool SparseMatrix::load(char* fname)
{
	ifstream ifs(fname, ios::in);
	if (!ifs.is_open() || ifs.fail())
	{
		cout << "failed to load " << fname << endl;
		return false;
	}

	cout << "load " << fname << endl;

	int nrows, ncols;
	ifs >> nrows >> ncols;
	if (nrows <= 0 || ncols <= 0)
	{
		cout << "invalid # of rows or columns! nrows = " << nrows << ", ncols = " << ncols << endl;
		return false;
	}

	m_rows = nrows;
	m_cols = ncols;

	int nz;
	ifs >> nz;

	int count(0);
	for (int i = 0; i < nz; i ++)
	{
		int row, col;
		double f;

		ifs >> row >> col >> f;

		if (row < 0 || row >= nrows || col < 0 || col >= ncols)
		{
			cout << "invalid element! row = " << row << ", col = " << col << ", f = " << f << endl;
		}
		else
		{
			insert(row, col, f);
			count ++;
		}
	}

	cout << "nrows = " << nrows << ", ncols = " << ncols << endl;
	cout << "load " << count << " non-zero elements" << endl;
	ifs.close();

	return true;
}

//////////////////////////////////////////// MOSEK /////////////////////////////////////////////
//

#include "mosek.h"

static void MSKAPI printstr(void *handle,
                            char str[])
{
  printf("%s",str);
} /* printstr */

// quadratic programming
bool mosek_qp(
		int NUMCON,
		int NUMANZ,
		int NUMVAR,
		int NUMQNZ,
		int* qsubi, // Q
		int* qsubj,
		double* qval,
		double* c,	// c
		double fix_term,
		//int* bkc,
		MSKboundkeye* bkc,
		double* blc,
		double* buc,
		int* ptrb,
		int* ptre,
		int* sub,
		double* val,
		//int* bkx,   // upper & lower bound
		MSKboundkeye* bkx,
		double* blx,
		double* bux,
		double* xx, // x
		double& objval,
		bool display_splash)
{
	clock_t start, finish;
	start = clock();

	int r;
	MSKenv_t  env;
	MSKtask_t task;

	/* Make the mosek environment. */
	r = MSK_makeenv(&env,NULL,NULL,NULL,NULL);

	/* Check whether the return code is ok. */
	if (r != MSK_RES_OK)
	{
		cout << "fail to make mosek enviroment" << endl;
		return false;
	}
	MSK_linkfunctoenvstream(env,MSK_STREAM_LOG,NULL,printstr);

	  /* Initialize the environment. */   
	  r = MSK_initenv(env);
	  if (r != MSK_RES_OK)
	  {
		  cout << "fail to init mosek enviroment" << endl;
		  return false;
	  }

	/* Directs the log stream to the user specified procedure 'printstr'. */
	

	/* Make the optimization task. */
	r = MSK_maketask(env,NUMCON,NUMVAR,&task);

	if (r != MSK_RES_OK)
	{
		cout << "fail to make mosek task" << endl;
		return false;
	}

	// input the linear part of an optimization task. 
	r = MSK_inputdata(task,
		                  NUMCON,NUMVAR,
                          NUMCON,NUMVAR,
                          c,fix_term,
                          ptrb,
                          ptre,
                          sub,
                          val,
                          bkc,
                          blc,
                          buc,
                          bkx,
                          blx,
                          bux);
	if (r != MSK_RES_OK)
	{
		cout << "fail to input data to MOSEK! returned value = " << r << endl;
		return false;
	}

	if (display_splash)
		cout << "calling MSK_putqobj" << endl;
	/* Input the Q for the objective. */
	r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval);

	if (r != MSK_RES_OK)
	{
		cout << "fail to input the Q for the objective" << endl;
		return false;
	}

	if (display_splash)
		cout << "calling MSK_optimize" << endl;
	r = MSK_optimize(task);
	
	if (r != MSK_RES_OK)
	{
		cout << "optimization failed! mosek returned value = " << r << endl;
		return false;
	}

	if (display_splash)
		cout << "calling MSK_getsolutionslice" << endl;
	MSK_getsolutionslice(task,
                         MSK_SOL_ITR, // basic/interior/integer solution
                         MSK_SOL_ITEM_XX,
                         0,
                         NUMVAR,
                         xx);
	
	MSK_deletetask(&task);
	MSK_deleteenv(&env);

	if (display_splash)
	{
		cout << "MSK released the task and enviroment" << endl;
		finish = clock();
		double duration = (double)(finish - start) / CLOCKS_PER_SEC;
		cout.precision(10);
		cout << "\ttime for mosek = " << duration << endl;
	}

	return true;
}


void MosekSolve( SparseMatrix * A, double * b )
{
	/*
		suppose A is a symetric, semi-positive definite matrix, 
		Solve AX=b in a convex quadratic programming 
			min 0.5*X^{T}*A*X - b*X
	*/
	if (A->rows() != A->cols())
	{
		cout << "fatal error in MosekSolve! A must be a squared matrix!" << endl;
		return;
	}

	for( int r = 0; r < A->rows(); r ++ )
	{
		b[r] = -b[r];
	}

	

	int NUMCON = 0;			// Number of constraints.             
	int NUMVAR = A->rows();	// Number of variables.               
	int NUMANZ = 0;  
	int NUMQNZ = 0; // Number of nonzeros in Q.  
	for( AVLTreeIterator<SparseMatrixElement> iter0( A->get_elements() ); !iter0.end(); ++ iter0 )
	{
		SparseMatrixElement * elem = *iter0;

		int row = elem->row();
		int col = elem->col();

		if (row >= col)
		{
			NUMQNZ ++;
		}
	}

	fprintf(stderr, "there are %d non-zeros in the lower triangle of A\n", NUMQNZ);
	
	int* qsubi = new int[NUMQNZ];
	int* qsubj = new int[NUMQNZ];
	double* qval = new double[NUMQNZ];

	//int* bkx = new int[NUMVAR];
	MSKboundkeye* bkx = new MSKboundkeye[NUMVAR];
	double* blx = new double[NUMVAR];
	double* bux = new double[NUMVAR];

	double* x = new double[NUMVAR];

	for (int i = 0; i < NUMVAR; i ++)
	{
		blx[i] = -MSK_INFINITY;
		bux[i] = MSK_INFINITY;
		bkx[i] = MSK_BK_FR;
	}

	int count(0);
	for( AVLTreeIterator<SparseMatrixElement> iter( A->get_elements() ); !iter.end(); ++ iter )
	{
		SparseMatrixElement * elem = *iter;

		int row = elem->row();
		int col = elem->col();

		if (row >= col)
		{
			qsubi[count] = row;
			qsubj[count] = col;
			qval[count] = elem->val();
				
			count ++;
		}
	}
	if (count != NUMQNZ)
	{
		fprintf(stderr, "fatal error in setup mosek!\n");
	}

	double objval, c_fix(0);

	if (!mosek_qp(0, 0, NUMVAR, NUMQNZ, qsubi, qsubj, qval, b, c_fix,
				NULL, NULL, NULL, NULL, NULL, NULL, NULL, // linear constraints
				bkx, blx, bux,   // upper & lower bound
				x, objval, true))
	{
		fprintf(stderr, "failed to solve linear equations!\n");
	}


	for(int r = 0; r < A->rows(); r ++ )
	{
		b[r] = x[r];
	}

	delete[] qsubi;
	delete[] qsubj;
	delete[] qval;
	delete[] bkx;
	delete[] blx;
	delete[] bux;
	delete[] x;
	delete[] b;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
void CGSolve(SparseMatrix* A, double* bb, int nrhs)
{
	cout << "CGSolve " << endl;
	int ITOL = 1;
	double TOL = 1e-12;
	int ITMAX = 500;

	unsigned int NP = A->rows();
	cout << "NP = " << NP << endl;

	unsigned nmax = NP*NP/20; // 5 percent 
	cout << "nmax = " << nmax << endl;

	unsigned long* ija = new unsigned long[nmax];
	double* sa = new double[nmax];
/*
	for (j=1;j<=n;j++) sa[j]=a[j][j];
	ija[1]=n+2;
	k=n+1;
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) {
			if (fabs(a[i][j]) >= thresh && i != j) {
				if (++k > nmax) nrerror("sprsin: nmax too small");
				sa[k]=a[i][j];
				ija[k]=j;
			}
		}
		ija[i+1]=k+1;
	}
*/
	sa[0] = 0;;
	ija[0] = 0;
	for (unsigned int j = 1; j <= NP; j ++) 
		sa[j] = A->element(j-1, j-1);

	ija[1] = NP + 2;
	unsigned int k = NP + 1;

	for (unsigned int i = 1; i <= NP;i ++)
	{
		if (i % 10 == 0 || i == NP)
		{
			fprintf(stderr, "\r formating... row %d %4f percent", i, 100*double(i)/NP);
			fflush(stderr);     
		}

		for (unsigned int j = 1; j <= NP; j ++) 
		{
			double e = A->element(i-1,j-1);

			if (fabs(e) >= 1e-6 && i != j) 
			{
				if (++k > nmax) 
				{
					nrerror("sprsin: nmax too small");
				}

				sa[k] = e;
				ija[k] = j;
			}
		}
		ija[i+1] = k+1;
	}
	cout << endl;

	cout << "sa.size() = " << k << endl;

	int iter;
	double *b,*x,err;
	b=dvector(1,NP);
	x=dvector(1,NP);

	cout << "nrhs = " << nrhs << endl;

	for (int j = 0; j < nrhs; j ++)
	{
		cout << "solving j = " << j << endl;

		for (unsigned int i = 1;i <= NP;i ++) 
		{
			x[i] = 1.0;
			b[i] = bb[j*NP+(i-1)];
		}

		linbcg(NP, b, x, ITOL, TOL, ITMAX, &iter, &err, ija, sa);

		printf("%s %15e\n","Estimated error:",err);
		printf("%s %6d\n","Iterations needed:",iter);

		for (unsigned int i = 1; i <= NP; i ++)
		{
			bb[j*NP+(i-1)] = x[i];
		}
	}

	free_dvector(x,1,NP);
	free_dvector(b,1,NP);
	
	delete[] ija;
	delete[] sa;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

#include "cghs.h"
#include "bicgsq.h"
#include "bicgstab.h"
#include "gmres.h"

void LSolve(SparseMatrix* A, double* bb, int nrhs)
{
	cout << "LSolve " << endl;

	int NP = A->rows();

	MyMatrix mm;
	mm.nr = NP;

	int count(0);
	for (int i = 0; i < NP;i ++)
	{
		if (i % 10 == 0 || i == NP-1)
		{
			fprintf(stderr, "\r formating... row %d %4f percent", i, 100*double(i+1)/NP);
			fflush(stderr);     
		}

		vector<double>* ee = new vector<double>;
		vector<int>* ii = new vector<int>;

		for (int j = 0; j < NP; j ++) 
		{
			double e = A->element(i,j);

			if (fabs(e) >= 1e-6) 
			{
				ee->push_back(e);
				ii->push_back(j);
				count ++;
			}
		}

		mm.elements.push_back(ee);
		mm.indices.push_back(ii);
	}
	cout << endl;

	cout << "nnz = " << count << endl;

	double* b = new double[NP];
	double* x = new double[NP];
	double* w = new double[NP];
	double* v = new double[NP];

	double EPS(1e-14);

	for (int j = 0; j < nrhs; j ++)
	{
		cout << "solving j = " << j << endl;

		for (int i = 0;i < NP; i ++) 
		{
			x[i] = 1.0;
			b[i] = bb[j*NP+i];
		}

		int its = bicgsq(NP, mm, b, x, EPS, true);
		cout << "iteration = " << its << endl;

		double max_error(0);
	    mult(mm, x, v);
		for (int i = 0; i < NP; i ++)
		{
			w[i] = fabs(v[i] - b[i]);
			if (w[i] > max_error)
			{
				max_error = w[i];
			}
		}
		cout << "max_error = " << max_error << endl;

		for (int i = 0; i < NP; i ++)
		{
			bb[j*NP+i] = x[i];
		}
	}


	delete[] b;
	delete[] x;

	for (int i = 0; i < NP; i ++)
	{
		delete mm.elements[i];
		delete mm.indices[i];
	}
}


void LSolve(MyMatrix& mm, double* bb, int nrhs)
{
	cout << "LSolve " << endl;

	int NP = mm.nr;

	double* b = new double[NP];
	double* x = new double[NP];
	double* w = new double[NP];
	double* v = new double[NP];

	double EPS(1e-14);

	for (int j = 0; j < nrhs; j ++)
	{
		cout << "=============================== solving equation " << j << " ==============================" << endl;

		for (int i = 0;i < NP; i ++) 
		{
			x[i] = 1.0;
			b[i] = bb[j*NP+i];
		}

		int its = bicgsq(NP, mm, b, x, EPS, true);
		cout << "iteration = " << its << endl;

		double max_error(0);
	    mult(mm, x, v);
		for (int i = 0; i < NP; i ++)
		{
			w[i] = fabs(v[i] - b[i]);
			if (w[i] > max_error)
			{
				max_error = w[i];
			}
		}
		cout << "max_error = " << max_error << endl;

		for (int i = 0; i < NP; i ++)
		{
			bb[j*NP+i] = x[i];
		}
	}


	delete[] b;
	delete[] x;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

/*
void LaplaceEigen(SparseMatrix* A, int ne, double* evalue, double* evec, int& nconv) // ne : # of eigenvalues
{
	if( A->cols() != A->rows() ) 
	{
		cout << "fatal error in SuperLUSolve! A->cols() = " << A->cols() << ", A->rows() = " << A->rows() << endl;
		return;
	}

	int n = A->cols();

	SM sm;

	A->convert( sm );

	ARluNonSymMatrix<double> matrix(sm.nr, sm.nnz, sm.buffer, sm.row_ind, sm.col_cnt);

	// Defining what we need: the four eigenvectors of A with largest magnitude.

	ARluNonSymStdEig<double> dprob(ne, matrix);

	// Finding eigenvalues and eigenvectors.

	dprob.FindEigenvectors();

	nconv = Prob.ConvergedEigenvalues();

	cout << "Dimension of the system            : " << Prob.GetN()    << endl;
	cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << endl;
	cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
	cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << endl;
	cout << "Number of iterations taken         : " << Prob.GetIter() << endl;
	cout << endl;

	if (Prob.EigenvaluesFound()) 
	{
		for (int i = 0; i < nconv; i ++) 
		{
			evalue[i] = Prob.EigenvalueReal(i);
			if (Prob.EigenvalueImag(i)!=0.0) 
			{
				cout << "no. " << i << " eigenvalue is complex! " << Prob.EigenvalueReal(i) << "+" << Prob.EigenvalueImag(i) << "i" << endl;
			}
			
			double* ev = Prob.RawEigenvector(i);
			for (int j = 0; j < n; j ++)
			{
				evec[i*n+j] = ev[j];
			}
		}
    }
}*/