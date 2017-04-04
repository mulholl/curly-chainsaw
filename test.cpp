#include <iostream>
#include <Eigen/Dense>
#include "ProtographClass.hpp"
#include "AListClass.hpp"
#include "PEXIT_AWGN.hpp"

using namespace std;

int main(void)
{
	Protograph P;

	Eigen::MatrixXi B(3, 6);
	Eigen::MatrixXi H;

	FILE *fp;

	fp = fopen("alist.txt", "w");

	// B << 1 , 2 , 2 , 1 , 1 , 0 , 2 , 1 , 2 , 1 , 1 , 0 , 1 , 1 , 0 , 0 , 0 , 1;

	// cout << "B = \n" << B << endl;

	P = Protograph("test4.txt");
	cout << "B = \n" << P.BaseMatrix() << endl;
	// exit(0);

	for (unsigned int i = 0; i < 1; i++)
	{
		try
		{
			P.lift(10, 1, H);
			cout << "H.rows() = " << H.rows() << "\t\tH.cols() = " << H.cols() << endl;
		}
		catch (exception &e)
		{
			std::cout << "Exception: " << e.what() << std::endl << "i = " << i << std::endl;;
			exit(0);
		}
	}

	alist_matrix A(H);

	write_alist(fp, &A);

	// double threshold = PEXIT_Threshold_AWGN(P, dBtoLin(0.0), dBtoLin(5.0), 10);

	// cout << "threshold = " << lintodB(threshold) << " dB (" << threshold << ") (" << dBtoLin(lintodB(threshold)) << ")" << endl;

	// double threshold2 = PEXIT_Sweep_AWGN(P, dBtoLin(0.0), dBtoLin(5.0), 0.02, 10);

	// cout << "threshold2 = " << lintodB(threshold2) << " dB (" << threshold2 << ") (" << dBtoLin(lintodB(threshold2)) << ")" << endl;
	
	return 0;
}