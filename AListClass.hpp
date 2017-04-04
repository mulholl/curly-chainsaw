#ifndef ALIST_CLASS_HPP
#define ALIST_CLASS_HPP

#include <Eigen/Dense>
#include <iostream>

class alist_matrix
{
	public:
		alist_matrix();
		alist_matrix(Eigen::MatrixXi &H);
		// ~alist_matrix();
		int N , M ;      /* size of the matrix */
		int **mlist;     /* list of integer coordinates in the m direction where the non-zero entries are */
		int **nlist;     /* list of integer coordinates in the n direction where the non-zero entries are */
		int *num_mlist;  /* weight of each row, m */
		int *num_nlist;  /* weight of each column n */
		int *l_up_to ;
		int *u_up_to ;
		int *norder ; /* Doesn't appear to be used at all - Ian */
		int biggest_num_m ;       /* actual biggest sizes */
		int biggest_num_n ; 
		int biggest_num_m_alloc ; /* sizes used for memory allocation */
		int biggest_num_n_alloc ; 
		int tot ;
		int same_length ;  /* whether all vectors in mlist and nlist have same length */
};


void write_alist(FILE *fp, alist_matrix *a );
void write_ivector( FILE *fp, int *m, int l1, int h1);
void write_imatrix( FILE *fp ,  int **m, int l1, int h1, int l2, int h2);

#endif