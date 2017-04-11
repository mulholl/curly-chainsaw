#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "ProtographClass.hpp"
#include "AListClass.hpp"
#include "PEXIT_AWGN.hpp"

#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

int main(void)
{
	const int max_iters = 10;

	// Protograph P;

	// Eigen::MatrixXi B(3, 6);
	// Eigen::MatrixXi H;

	// FILE *fp;
	// B << 1 , 2 , 2 , 1 , 1 , 0 , 2 , 1 , 2 , 1 , 1 , 0 , 1 , 1 , 0 , 0 , 0 , 1;

	// cout << "B = \n" << B << endl;

	// P = Protograph("test4.txt");
	// cout << "B = \n" << P.BaseMatrix() << endl;
	// // exit(0);

	// for (unsigned int i = 0; i < 1; i++)
	// {
	// 	try
	// 	{
	// 		P.lift(1250, 1, H);
	// 		cout << "H.rows() = " << H.rows() << "\t\tH.cols() = " << H.cols() << endl;
	// 	}
	// 	catch (exception &e)
	// 	{
	// 		std::cout << "Exception: " << e.what() << std::endl << "i = " << i << std::endl;;
	// 		exit(0);
	// 	}
	// }


	// fp = fopen("alist.txt", "w");

	// alist_matrix A(H);

	// write_alist(fp, &A);
	// fclose(fp);

	FILE *fp1;

	fp1 = fopen("BER.txt", "w");

	itpp::LDPC_Parity mat;
	mat.load_alist("alist.txt");

	LDPC_Code code(&mat, 0, false);
	code.set_exit_conditions(max_iters, true, true);

	cout << "Code rate is " << code.get_rate() << endl;


	//Vectors
	bvec bits_in, encodedBits, dec_bits;
	vec symbols, rec;
	//Classes
	BPSK bpsk;  //The BPSK modulator/debodulator class
	BERC berc;  //The Bit Error Rate Counter class
	//Randomize the random number generator
	RNG_randomize();

	float SNRdB, SNR_lin; /* 	Eb/N0, the ratio of Eb, the energy per bit, to N0, the noise
								power spectral density, in decibal and linear form,
								respectively */

	
	int N;
	double N0; /* Noise power spectral density */
	double sigma, sigma_squared; /* Noise standard deviation and variance, respectively */

	//Init
	N = mat.get_nvar(); //The number of bits to simulate
	bits_in = zeros_b(N);
	//Do the BPSK modulation
	bpsk.modulate_bits(bits_in, symbols);

	int CWs_to_do = 100000;
	int max_CW_errs = 1000;
	int CW_errs;
	int CW;

	float R = code.get_rate();
	R = 1; /* If we want to simulate uncoded transmission, set the code rate to 1 */

	float CER = 0.0;

	QLLRvec LLRin, LLRout;
	vec softbits;

	double errs_prev, errs_curr = 0;

	vector<float> SNRdB_vals = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

	for (vector<float>::iterator it = SNRdB_vals.begin(); it < SNRdB_vals.end(); it++)
	// for (int SNRdB_int = 0; SNRdB_int <= 50; SNRdB_int += 5)
	{
		SNRdB = *it;//(float)SNRdB_int / 10;
		cout << "\n\nSNRdB = " << SNRdB << endl;

		SNR_lin = pow(10, SNRdB/10);

		N0 = 1 / (R * SNR_lin);
		sigma_squared = N0 / 2;
		sigma = sqrt(sigma_squared);
		cout << "N0 = " << N0 << endl;

		berc.clear();

		errs_prev = 0;
		errs_curr = 0;
		CW_errs = 0;

		for (CW = 0; CW < CWs_to_do; CW++)
		{
			cout << CW << " " << flush;
			//Add the AWGN
			rec = symbols + sigma * randn(N);
			//Decode the received bits
			// softbits = bpsk.demodulate_soft_bits(rec, N0);
			// code.bp_decode(code.get_llrcalc().to_qllr(softbits), LLRout);
			// dec_bits = LLRout < 0;

			bpsk.demodulate_bits(rec, dec_bits); // Uncoded

			//Count the number of errors
			berc.count(bits_in, dec_bits);

			errs_curr = berc.get_errors();	

			if (errs_curr > errs_prev)
			{
				CW_errs++;
				errs_prev = errs_curr;
				CER = (float)CW_errs / (CW+1);
				if (CW_errs >= max_CW_errs){
					break;
				}
			}

			cout << "(" << CW_errs << ") " << flush;
		}

		//Print the results
		cout << "\n\nCW = " << CW << ", CW_errs = " << CW_errs << endl;
		cout << "There were " << berc.get_errors() << " received bits in error." << endl;
		cout << "There were " << berc.get_corrects() << " correctly received bits." << endl;
		cout << "The error probability was " << berc.get_errorrate() << endl;
		cout << "The theoretical error probability is " << 0.5*erfc(sqrt(SNR_lin)) << endl;
		fprintf(fp1, "%5.2f %10e %10e\n", SNRdB, berc.get_errorrate(), CER);
		fflush(fp1);
	}

	fclose(fp1);
		
	return 0;
}