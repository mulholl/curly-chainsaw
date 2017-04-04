#include "AListClass.hpp"

using namespace std;
using namespace Eigen;

alist_matrix::alist_matrix(MatrixXi &H)
{
	N = H.cols();
	M = H.rows();

	mlist = new int *[M];
	nlist = new int *[N];

	num_mlist = new int[M]{0};
	num_nlist = new int[N]{0};

	int count;

	for (int m = 0; m < M; m++)
	{
		for (int n = 0; n < N; n++)
		{		
			if (H(m, n))	
			{
				num_mlist[m] += H(m, n);
				num_nlist[n] += H(m, n);
			}
		}
	}

	biggest_num_m = 0;

	for (int m = 0; m < M; m++)
	{
		biggest_num_m = (num_mlist[m] > biggest_num_m) ? num_mlist[m] : biggest_num_m;
	}

	biggest_num_n = 0;

	for (int n = 0; n < N; n++)
	{
		biggest_num_n = (num_nlist[n] > biggest_num_n) ? num_nlist[n] : biggest_num_n;
	}

	for (int m = 0; m < M; m++)
	{
		mlist[m] = new int[biggest_num_m]{0};

		count = 0;

		for (int n = 0; n < N; n++)
		{
			if (H(m, n))
			{
				cout << "mlist[" << m << "][" << count << "] = " << n + 1 << endl;
				mlist[m][count++] = n + 1;
				if (count >= num_mlist[m])
				{
					cout << "Breaking: count = " << count << endl;
					break;
				}
			}
		}
	}

	for (int n = 0; n < N; n++)
	{
		nlist[n] = new int[biggest_num_n]{0};

		count = 0;

		for (int m = 0; m < M; m++)
		{
			if (H(m, n))
			{
				nlist[n][count++] = m + 1;
				if (count >= num_nlist[n])
				{
					break;
				}
			}
		}
	}	

	biggest_num_m_alloc = biggest_num_m;
	biggest_num_n_alloc = biggest_num_n;

	same_length = 1;
}

// alist_matrix::~alist_matrix()
// {
// 	for (int n = 0; n < N; n++)
// 	{
// 		delete nlist[n];
// 	}
// 	delete nlist;
// 	for (int m = 0; m < M; m++)
// 	{
// 		delete mlist[m];
// 	}
// 	delete mlist;
// 	delete num_nlist;
// 	delete num_mlist;
// }

void write_alist ( FILE *fp , alist_matrix *a ) {
  /* this assumes that mlist and nlist have the form of a rectangular
     matrix in the file; if lists have unequal lengths, then the 
     entries should be present (eg zero values) but are ignored
     */
  int N = a->N , M = a->M ;

  fprintf ( fp , "%d %d\n" , N , M ) ; 
  fprintf ( fp , "%d %d\n" , a->biggest_num_n , a->biggest_num_m ) ; 
  write_ivector ( fp , a->num_nlist , 0 , N-1 ) ; 
  write_ivector ( fp , a->num_mlist , 0 , M-1 ) ; 
  write_imatrix ( fp , a->nlist , 0 , N-1 , 0 , a->biggest_num_n - 1) ; 
  write_imatrix ( fp , a->mlist , 0 , M-1 , 0 , a->biggest_num_m - 1) ; 
}  

void write_ivector( FILE *fp, int *m, int l1, int h1)
{
  int i;
  
  for (i=l1; i<=h1; i++)
  	fprintf( fp , "%d " , m[i] );
  fprintf( fp , "\n");
}

void write_imatrix( FILE *fp ,  int **m, int l1, int h1, int l2, int h2)
{
  int i,j;
  
  for (i=l1; i<=h1; i++){
  	fprintf(fp, "i = %d: ", i);
    for (j=l2; j<=h2; j++)
    {
    	fprintf( fp , "%d ",m[i][j]);
  	}
    fprintf(fp , "\n");
  }                   
}