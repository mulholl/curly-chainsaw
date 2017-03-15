#include "PEXIT_AWGN.hpp"

/* Use PEXIT analysis to find the decoding threshold of for a protograph over the AWGN
	Arguments:
		Input:
			const Protograph &P 			-	The protograph to be used
			const double &EbN0_linear_min 	-	The lower boundary of the Eb/N0 range
			const double &EbN0_linear_max 	-	The upper boundary of the Eb/N0 range
			const unsigned int &max_its 	-	The maximum number of iterations allowed
	Return value:
		double 	-	-1.0 if no threshold was found within the range, otherwise the threshold 
					value is returned
 */
double PEXIT_Threshold_AWGN(const Protograph &P, const double &EbN0_linear_min, const double &EbN0_linear_max, const unsigned int &max_its)
{
	unsigned int its_taken;
	double MI_reached;

	double EbN0_linear_upper, EbN0_linear_lower, EbN0_linear_mid, EbN0_linear_mid_old, EbN0_linear_best;	

	if (EbN0_linear_max >= EbN0_linear_min)
	{
		EbN0_linear_upper = EbN0_linear_max;
		EbN0_linear_lower = EbN0_linear_min;
	}
	else
	{
		EbN0_linear_upper = EbN0_linear_min;
		EbN0_linear_lower = EbN0_linear_max;
	}

	/* If convergence occurs at the minimum Eb/N0 value, then we don't need to do the 
	search as we won't find a better Eb/N0 value within the range */
	if (ret = PEXIT_AWGN(P, EbN0_linear_lower, max_its, its_taken, MI_reached))
	{
		return EbN0_linear_lower;
	}

	/* If convergence fails to occur at the maximum Eb/N0 value, then we don't need to
	do the search as we won't find any Eb/N0 value within the range for which convergence
	occurs */
	if (!(PEXIT_AWGN(P, EbN0_linear_upper, max_its, its_taken, MI_reached)))
	{
		return -1.0;
	}

	EbN0_linear_best = EbN0_linear_upper;
	EbN0_linear_mid_old = EbN0_linear_best;
	EbN0_linear_mid = (EbN0_linear_lower + EbN0_linear_upper) / 2;

	while (fabs(EbN0_linear_best - EbN0_linear_mid) > PEXIT_AWGN_TOLERANCE)
	{
		if (PEXIT_AWGN(P, EbN0_linear_mid, max_its, its_taken, MI_reached))
		{
			EbN0_linear_upper = EbN0_linear_mid;
			EbN0_linear_best = EbN0_linear_mid;
		}
		else
		{
			EbN0_linear_lower = EbN0_linear_mid;
		}

		EbN0_linear_mid = (EbN0_linear_lower + EbN0_linear_upper) / 2;
	}
	
	return EbN0_linear_best;
}

/* Use PEXIT analysis to determine whether a code lifted from a given protograph will 
	converge for the specified Eb/N0 value
	Arguments:
		Input:
			const Protograph &P 			-	The protograph to be used
			const double &EbN0_linear 		-	The Eb/N0 value to be used (linear)
			const unsigned int &max_its 	-	The maximum number of iterations allowed
		Output:
			unsigned int &its_taken 		-	The actual number of iterations required for
												convergence at this Eb/N0 value
			double &MI_reached 				-	The average value of I_APP after the 
												final iteration at this Eb/N0 value
	Return value:
		bool 	- 	TRUE if the code converges for this Eb/N0 value, FALSE otherwise
 */
bool PEXIT_AWGN(const Protograph &P, const double &EbN0_linear, const unsigned int &max_its, unsigned int &its_taken, double &MI_reached)
{
	Eigen::MatrixXi<double> IAC, IAV, IEC, IEV;
	Eigen::VectorXi<double> sigma_ch_squared, IAPP;
	double IAPP_min;

	PEXIT_AWGN_Initialize(P, sigma_ch_squared, IAV, IEV, IAC, IEC, EbN0_linear);

	for (its_taken = 1; current_it <= max_its; current_it++)
	{
		PEXIT_AWGN_IEV(P, sigma_ch_squared, IAV, IEV);

		IAC = IEV;

		PEXIT_AWGN_IEC(P, IAC, IEC);

		IAV = IEC;

		if (PEXIT_AWGN_IAPP(P, sigma_ch_squared, IAV, IAPP, MI_reached, IAPP_min))
		{
			return TRUE;
		}
	}
	
}

/* Use PEXIT analysis to determine whether a code lifted from a given protograph will 
	converge for the specified Eb/N0 value
	Arguments:
		Input:
			const Protograph &P 			-	The protograph to be used
			const double &EbN0_linear 		-	The Eb/N0 value to be used (linear)
			const unsigned int &max_its 	-	The maximum number of iterations allowed
	Return value:
		bool 	- 	TRUE if the code converges for this Eb/N0 value, FALSE otherwise
 */
bool PEXIT_AWGN(const Protograph &P, const double &EbN0_linear, const unsigned int &max_its)
{
	unsigned int *dummy_i;
	double *dummy_d;

	return PEXIT_AWGN(P, EbN0_linear, max_its, dummy_i, dummy_d);
}

/* Initialize the variables that will be used in the PEXIT analysis
	Arguments:
		Input:
			const Protograph &P 						-	The protograph to be used
			const Eigen::VectorXi<double> EbN0_linear 	-	Vector containing the Eb/N0 values
															for each bit
			const double R 								- 	The code rate
			Eigen::VectorXi<double> &sigma_ch_squared	-	Vector containing noise variance
															values for each bit
			Eigen::MatrixXi<double> &IAV 				-	Matrix to store a priori MI values 													into the VNs
			Eigen::MatrixXi<double> &IEV 				-	Matrix to store a posteriori 														extrinsic MI values from the VNs
			Eigen::MatrixXi<double> &IAC 				-	Matrix to store a priori MI values 													into the CNs
			Eigen::MatrixXi<double> &IEC 				-	Matrix to store a posteriori 														extrinsic MI values from the CNs
			Eigen::VectorXi<double> &IAPP 				-	Vector to store the MI between 
															the a posteriori probability 
															log-likelihood ratios and 
															the associated codeword bits		
 */
void PEXIT_AWGN_Initialize(const Protograph &P, const Eigen::VectorXi<double> EbN0_linear, const double R, Eigen::VectorXi<double> &sigma_ch_squared, Eigen::MatrixXi<double> &IAV, Eigen::MatrixXi<double> &IEV, Eigen::MatrixXi<double> &IAC, Eigen::MatrixXi<double> &IEC, Eigen::VectorXi<double> &IAPP)
{
	unsigned int N = P.VNs;
	unsigned int M = P.CNs;

	sigma_ch_squared.resize(N);
	for (unsigned int j = 0; j <= N-1; j++)
	{
		Ich(j) = 8 * R * EbN0_linear(j);
	}
	IAV.resize(M, N);
	IEV.resize(M, N);
	IAC.resize(M, N);
	IEC.resize(M, N);
	IAPP.resize(N);
}

/* Initialize the variables that will be used in the PEXIT analysis
	Arguments:
		Input:
			const Protograph &P 						-	The protograph to be used
			const Eigen::VectorXi<double> EbN0_linear 	-	Eb/N0 value to be used for all bits
															(linear)
			const double R 								- 	The code rate
			Eigen::VectorXi<double> &sigma_ch_squared	-	Vector containing noise variance
															values for each bit
			Eigen::MatrixXi<double> &IAV 				-	Matrix to store a priori MI values
															into the VNs
			Eigen::MatrixXi<double> &IEV 				-	Matrix to store a posteriori 														extrinsic MI values from the VNs
			Eigen::MatrixXi<double> &IAC 				-	Matrix to store a priori MI values
															into the CNs
			Eigen::MatrixXi<double> &IEC 				-	Matrix to store a posteriori 														extrinsic MI values from the CNs
			Eigen::VectorXi<double> &IAPP 				-	Vector to store the MI between 
															the a posteriori probability 
															log-likelihood ratios and 
															the associated codeword bits		
 */
void PEXIT_AWGN_Initialize(const Protograph &P, const double EbN0_linear, const double R, Eigen::VectorXi<double> &sigma_ch_squared, Eigen::MatrixXi<double> &IAV, Eigen::MatrixXi<double> &IEV, Eigen::MatrixXi<double> &IAC, Eigen::MatrixXi<double> &IEC, Eigen::VectorXi<double> &IAPP)
{
	unsigned int N = P.VNs;
	unsigned int M = P.CNs;

	Ich.resize(N);
	sigma_ch_squared = Eigen::VectorXi.Constant(8 * R * EbN0_linear);
	IAV.resize(M, N);
	IEV.resize(M, N);
	IAC.resize(M, N);
	IEC.resize(M, N);
}

/* Calculate the extrinsic MI values from the VN side
	Arguments:
		Input:
			const Protograph &P 						-	The protograph to be used
			Eigen::VectorXi<double> &sigma_ch_squared	-	Vector containing noise variance
															values for each bit
			Eigen::MatrixXi<double> &IAV 				-	Matrix storing a priori MI values
															into the VNs
		Output:
			Eigen::MatrixXi<double> &IEV 				-	Matrix to store a posteriori 														extrinsic MI values from the VNs
*/
void PEXIT_AWGN_IEV(const Protograph &P, const Eigen::VectorXi<double> &sigma_ch_squared, const Eigen::MatrixXi<double> &IAV, Eigen::MatrixXi<double> &IEV)
{
	unsigned int N = P.VNs;
	unsigned int M = P.CNs;

	double term1, term2;
	double b_ij;

	for (unsigned int j = 0; j <= N-1; j++)
	{
		for (unsigned int i = 0; i <= M-1; i++)
		{
			b_ij = P.BaseMatrix(i, j);

			if (b_ij != 0)
			{
				term1 = 0.0;

				for (unsigned int s = 0; s <= M-1; s++)
				{
					if (s != i)
					{
						term1 += P.BaseMatrix(s, j) * pow(Jinv(IAV(s, j), 2);
					}
					else if (b_ij > 1)
					{
						term1 += (b_ij - 1) * pow(Jinv(IAV(i, j), 2);
					}
				}

				term2 = sigma_ch_squared(j);

				IEV(i, j) = J(sqrt(term1 + term2));
			}
			else
			{
				IEV(i, j) = 0.0;
			}
		}
	}
}

/* Calculate the extrinsic MI values from the CN side
	Arguments:
		Input:
			const Protograph &P 						-	The protograph to be used
			Eigen::MatrixXi<double> &IAC 				-	Matrix storing a priori MI values
															into the CNs
		Output:
			Eigen::MatrixXi<double> &IEC 				-	Matrix to store a posteriori 														extrinsic MI values from the CNs
*/
void PEXIT_AWGN_IEC(const Protograph &P, const Eigen::MatrixXi<double> &IAC, Eigen::MatrixXi<double> &IEC)
{
	unsigned int N = P.VNs;
	unsigned int M = P.CNs;

	double term1;
	double b_ij;

	for (unsigned int j = 0; j <= N-1; j++)
	{
		for (unsigned int i = 0; i <= M-1; i++)
		{
			b_ij = P.BaseMatrix(i, j);

			if (b_ij != 0)
			{
				term1 = 0.0;

				for (unsigned int s = 0; s <= N-1; s++)
				{
					if (s != j)
					{
						term1 += P.BaseMatrix(i, s) * pow(Jinv(1 - IAC(i, s), 2);
					}
					else if (b_ij > 1)
					{
						term1 += (b_ij - 1) * pow(Jinv(1 - IAC(i, j), 2);
					}
				}

				IEC(i, j) = 1 - J(sqrt(term1));
			}
			else
			{
				IEC(i, j) = 0.0;
			}
		}
	}
}

/* Calculate the MI between the a posteriori probability log-likelihood ratios and the 
associated codeword bits
	Arguments:
		Input:
			const Protograph &P 						-	The protograph to be used
			Eigen::VectorXi<double> &sigma_ch_squared	-	Vector containing noise variance
															values for each bit
			Eigen::MatrixXi<double> &IAV 				-	Matrix storing a priori MI values
															into the VNs
		Output:
			Eigen::VectorXi<double> &IAPP 				-	Vector to store the MI between 
															the a posteriori probability 
															log-likelihood ratios and 
															the associated codeword bits
			double &IAPP_avg							-	The average value of the vector IAPP
			double &IAPP_min							-	The minimum value of the vector IAPP
	Return value:
		bool	-	TRUE if IAPP converges to 1 (to within the tolerance value given by
		 			PEXIT_AWGN_TOLERANCE), FALSE otherwise
*/
bool PEXIT_AWGN_IAPP(const Protograph &P, const Eigen::VectorXi<double> &sigma_ch_squared, const Eigen::MatrixXi<double> &IAV, Eigen::VectorXi<double> &IAPP, double &IAPP_avg, double &IAPP_min)
{
	double term1, term2;
	double IAPP_j;
	unsigned int N = P.VNs;
	bool ret = TRUE;
	IAPP_avg = 0.0;

	IAPP.resize(N);
	for (unsigned int j = 0; j <= N-1; j++)
	{
		float term2 = sigma_ch_squared;

		term1 = 0.0;

		for (unsigned int s = 0; s <= P.CNs-1; s++)
		{
			term1 += P.BaseMatrix(s, j) * pow(Jinv((IAV(s, j))), 2;
		}

		IAPP(j) = sqrt(term1 + term2);
		IAPP_avg += IAPP(j);

		ret = (IAPP_j < (1 - PEXIT_AWGN_TOLERANCE)) ? FALSE : ret;
	}

	IAPP_avg /= N;

	return ret;
}