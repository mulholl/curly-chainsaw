#include "PEXIT_AWGN.hpp"

std::vector<double> JLUT, JinvLUT;

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

	double EbN0_linear_upper, EbN0_linear_lower, EbN0_linear_mid, EbN0_linear_best;	

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
	if (PEXIT_AWGN(P, EbN0_linear_lower, max_its, its_taken, MI_reached))
	{
		return EbN0_linear_lower;
	}

	/* If convergence fails to occur at the maximum Eb/N0 value, then we don't need to
	do the search as we won't find any Eb/N0 value within the range for which convergence
	occurs */
	if (!(PEXIT_AWGN(P, EbN0_linear_upper, max_its, its_taken, MI_reached)))
	{
		std::cout << "Exiting threshold search early" << std::endl;
		return -1.0;
	}

	EbN0_linear_best = EbN0_linear_upper;
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

double PEXIT_Sweep_AWGN(const Protograph &P, const double &EbN0_linear_min, const double &EbN0_linear_max, const double &step, const unsigned int &max_its)
{
	unsigned int its_taken;
	double MI_reached;

	double EbN0_linear_upper, EbN0_linear_lower, EbN0_linear_best;	

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

	// /* If convergence occurs at the minimum Eb/N0 value, then we don't need to do the 
	// search as we won't find a better Eb/N0 value within the range */
	// if (PEXIT_AWGN(P, EbN0_linear_lower, max_its, its_taken, MI_reached))
	// {
	// 	return EbN0_linear_lower;
	// }

	//  If convergence fails to occur at the maximum Eb/N0 value, then we don't need to
	// do the search as we won't find any Eb/N0 value within the range for which convergence
	// occurs 
	// if (!(PEXIT_AWGN(P, EbN0_linear_upper, max_its, its_taken, MI_reached)))
	// {
	// 	std::cout << "Exiting threshold search early" << std::endl;
	// 	return -1.0;
	// }

	bool found_best = false;
	bool ret;

	unsigned int numVals = (EbN0_linear_upper - EbN0_linear_lower) / step;

	std::vector<float> MI_vals(numVals, 0.0);
	unsigned int count = 0;

	EbN0_linear_best = -1.0;

	std::cout << "step = " << step << " ( = " << lintodB(step) << " dB)" << std::endl;
	std::cout << "EbN0_linear_lower = " << EbN0_linear_lower << " ( = " << lintodB(EbN0_linear_lower) << " dB)" << std::endl;
	std::cout << "EbN0_linear_upper = " << EbN0_linear_upper << " ( = " << lintodB(EbN0_linear_upper) << " dB)" << std::endl;

	for (double current_EbN0_linear = EbN0_linear_lower; current_EbN0_linear < EbN0_linear_upper; current_EbN0_linear += step)
	{
		std::cout << "current_EbN0_linear = " << current_EbN0_linear << " ( = " << lintodB(current_EbN0_linear) << " dB)";
		ret = PEXIT_AWGN(P, current_EbN0_linear, max_its, its_taken, MI_reached);
		if (ret)
		{
			std::cout << "\t*";
		}
		else
		{
			std::cout << "\t\t\t" << MI_reached;
		}
		MI_vals[count++] = MI_reached;
		if (ret && !found_best)
		{
			found_best = true;
			EbN0_linear_best = current_EbN0_linear;
		}
		std::cout << std::endl;
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
		bool 	- 	true if the code converges for this Eb/N0 value, false otherwise
 */
bool PEXIT_AWGN(const Protograph &P, const double &EbN0_linear, const unsigned int &max_its, unsigned int &its_taken, double &MI_reached)
{
	Eigen::MatrixXd IAC, IAV, IEC, IEV;
	Eigen::VectorXd sigma_ch_squared, IAPP;
	double IAPP_min;
	double R = P.Rate();

	PEXIT_AWGN_Initialize(P, EbN0_linear, R, sigma_ch_squared, IAV, IEV, IAC, IEC, IAPP);

	for (its_taken = 1; its_taken <= max_its; its_taken++)
	{
		// std::cout << "\nIAV = \n" << IAV << std::endl;

		PEXIT_AWGN_IEV(P, sigma_ch_squared, IAV, IEV);
		
		// std::cout << "\nIEV = \n" << IEV << std::endl;

		IAC = IEV;

		// std::cout << "\nIAC = \n" << IAC << std::endl;

		PEXIT_AWGN_IEC(P, IAC, IEC);

		// std::cout << "\nIEC = \n" << IEC << std::endl;

		IAV = IEC;

		if (PEXIT_AWGN_IAPP(P, sigma_ch_squared, IAV, IAPP, MI_reached, IAPP_min))
		{
			return true;
		}
	}

	return false;	
}

/* Use PEXIT analysis to determine whether a code lifted from a given protograph will 
	converge for the specified Eb/N0 value
	Arguments:
		Input:
			const Protograph &P 			-	The protograph to be used
			const double &EbN0_linear 		-	The Eb/N0 value to be used (linear)
			const unsigned int &max_its 	-	The maximum number of iterations allowed
	Return value:
		bool 	- 	true if the code converges for this Eb/N0 value, false otherwise
 */
bool PEXIT_AWGN(const Protograph &P, const double &EbN0_linear, const unsigned int &max_its)
{
	unsigned int dummy_i;
	double dummy_d;

	return PEXIT_AWGN(P, EbN0_linear, max_its, dummy_i, dummy_d);
}

/* Initialize the variables that will be used in the PEXIT analysis
	Arguments:
		Input:
			const Protograph &P 						-	The protograph to be used
			const Eigen::VectorXd EbN0_linear 	-	Vector containing the Eb/N0 values
															for each bit
			const double R 								- 	The code rate
			Eigen::VectorXd &sigma_ch_squared	-	Vector containing noise variance
															values for each bit
			Eigen::MatrixXd &IAV 				-	Matrix to store a priori MI values 													into the VNs
			Eigen::MatrixXd &IEV 				-	Matrix to store a posteriori 														extrinsic MI values from the VNs
			Eigen::MatrixXd &IAC 				-	Matrix to store a priori MI values 													into the CNs
			Eigen::MatrixXd &IEC 				-	Matrix to store a posteriori 														extrinsic MI values from the CNs
			Eigen::VectorXd &IAPP 				-	Vector to store the MI between 
															the a posteriori probability 
															log-likelihood ratios and 
															the associated codeword bits		
 */
void PEXIT_AWGN_Initialize(const Protograph &P, const Eigen::VectorXd EbN0_linear, const double R, Eigen::VectorXd &sigma_ch_squared, Eigen::MatrixXd &IAV, Eigen::MatrixXd &IEV, Eigen::MatrixXd &IAC, Eigen::MatrixXd &IEC, Eigen::VectorXd &IAPP)
{
	unsigned int N = P.VNs();
	unsigned int M = P.CNs();

	sigma_ch_squared.resize(N);
	for (unsigned int j = 0; j <= N-1; j++)
	{
		sigma_ch_squared(j) = 8 * R * EbN0_linear(j);
	}
	IAV = Eigen::MatrixXd::Zero(M, N);
	IEV = Eigen::MatrixXd::Zero(M, N);
	IAC = Eigen::MatrixXd::Zero(M, N);
	IEC = Eigen::MatrixXd::Zero(M, N);
	IAPP = Eigen::VectorXd::Zero(N);

	createJLUT(JLUT);
	createJinvLUT(JinvLUT);
}

/* Initialize the variables that will be used in the PEXIT analysis
	Arguments:
		Input:
			const Protograph &P 				-	The protograph to be used
			const Eigen::VectorXd EbN0_linear 	-	Eb/N0 value to be used for all bits
													(linear)
			const double R 						- 	The code rate
			Eigen::VectorXd &sigma_ch_squared	-	Vector containing noise variance
													values for each bit
			Eigen::MatrixXd &IAV 				-	Matrix to store a priori MI values
													into the VNs
			Eigen::MatrixXd &IEV 				-	Matrix to store a posteriori 														extrinsic MI values from the VNs
			Eigen::MatrixXd &IAC 				-	Matrix to store a priori MI values
													into the CNs
			Eigen::MatrixXd &IEC 				-	Matrix to store a posteriori 														extrinsic MI values from the CNs
			Eigen::VectorXd &IAPP 				-	Vector to store the MI between 
													the a posteriori probability 				log-likelihood ratios and the 
													associated codeword bits		
 */
void PEXIT_AWGN_Initialize(const Protograph &P, const double EbN0_linear, const double R, Eigen::VectorXd &sigma_ch_squared, Eigen::MatrixXd &IAV, Eigen::MatrixXd &IEV, Eigen::MatrixXd &IAC, Eigen::MatrixXd &IEC, Eigen::VectorXd &IAPP)
{
	unsigned int N = P.VNs();
	unsigned int M = P.CNs();

	sigma_ch_squared.resize(N);
	sigma_ch_squared = Eigen::VectorXd::Constant(N, 8 * R * EbN0_linear);
	IAV = Eigen::MatrixXd::Zero(M, N);
	IEV = Eigen::MatrixXd::Zero(M, N);
	IAC = Eigen::MatrixXd::Zero(M, N);
	IEC = Eigen::MatrixXd::Zero(M, N);
	IAPP = Eigen::VectorXd::Zero(N);

	createJLUT(JLUT);
	createJinvLUT(JinvLUT);
}

/* Calculate the extrinsic MI values from the VN side
	Arguments:
		Input:
			const Protograph &P 						-	The protograph to be used
			Eigen::VectorXd &sigma_ch_squared	-	Vector containing noise variance
															values for each bit
			Eigen::MatrixXd &IAV 				-	Matrix storing a priori MI values
															into the VNs
		Output:
			Eigen::MatrixXd &IEV 				-	Matrix to store a posteriori 														extrinsic MI values from the VNs
*/
void PEXIT_AWGN_IEV(const Protograph &P, const Eigen::VectorXd &sigma_ch_squared, const Eigen::MatrixXd &IAV, Eigen::MatrixXd &IEV)
{
	unsigned int N = P.VNs();
	unsigned int M = P.CNs();

	double term1, term2;
	unsigned int b_ij;

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
						term1 += P.BaseMatrix(s, j) * pow(getJinv(IAV(s, j)), 2);
					}
					else if (b_ij > 1)
					{
						term1 += (b_ij - 1) * pow(getJinv(IAV(i, j)), 2);
					}
				}

				term2 = sigma_ch_squared(j);

				IEV(i, j) = getJ(sqrt(term1 + term2));
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
			Eigen::MatrixXd &IAC 				-	Matrix storing a priori MI values
															into the CNs
		Output:
			Eigen::MatrixXd &IEC 				-	Matrix to store a posteriori 														extrinsic MI values from the CNs
*/
void PEXIT_AWGN_IEC(const Protograph &P, const Eigen::MatrixXd &IAC, Eigen::MatrixXd &IEC)
{
	unsigned int N = P.VNs();
	unsigned int M = P.CNs();

	double term1;
	unsigned int b_ij;

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
						term1 += P.BaseMatrix(i, s) * pow(getJinv(1 - IAC(i, s)), 2);
					}
					else if (b_ij > 1)
					{
						term1 += (b_ij - 1) * pow(getJinv(1 - IAC(i, j)), 2);
					}
				}

				IEC(i, j) = 1 - getJ(sqrt(term1));
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
			Eigen::VectorXd &sigma_ch_squared	-	Vector containing noise variance
															values for each bit
			Eigen::MatrixXd &IAV 				-	Matrix storing a priori MI values
															into the VNs
		Output:
			Eigen::VectorXd &IAPP 				-	Vector to store the MI between 
															the a posteriori probability 
															log-likelihood ratios and 
															the associated codeword bits
			double &IAPP_avg							-	The average value of the vector IAPP
			double &IAPP_min							-	The minimum value of the vector IAPP
	Return value:
		bool	-	true if IAPP converges to 1 (to within the tolerance value given by
		 			PEXIT_AWGN_TOLERANCE), false otherwise
*/
bool PEXIT_AWGN_IAPP(const Protograph &P, const Eigen::VectorXd &sigma_ch_squared, const Eigen::MatrixXd &IAV, Eigen::VectorXd &IAPP, double &IAPP_avg, double &IAPP_min)
{
	double term1;
	unsigned int N = P.VNs();
	bool ret = true;
	IAPP_avg = 0.0;

	IAPP.resize(N);
	for (unsigned int j = 0; j <= N-1; j++)
	{
		double term2 = sigma_ch_squared[j];

		term1 = 0.0;

		for (unsigned int s = 0; s <= P.CNs()-1; s++)
		{
			term1 += P.BaseMatrix(s, j) * pow(getJinv(IAV(s, j)), 2);
		}

		IAPP(j) = getJ(sqrt(term1 + term2));
		IAPP_avg += IAPP(j);

		ret = (IAPP(j) >= (1 - PEXIT_AWGN_TOLERANCE)) ? ret : false;
	}

	IAPP_avg /= N;

	return ret;
}

double dBtoLin(const double &dB)
{
	return pow(10, dB / 10);
}

double lintodB(const double &lin)
{
	return 10 * log10(lin);
}

void createJLUT(std::vector<double> &JLUT)
{
	unsigned long int numels = 10 * (unsigned long int)JLUT_MULTIPLIER;
	JLUT.resize(numels);

	unsigned long int stop_ind = SIGMA_STAR * JLUT_MULTIPLIER;

	const double a_J1 = -0.0421061;
	const double b_J1 = 0.209252;
	const double c_J1 = -0.00640081;

	const double a_J2 = 0.00181491;
	const double b_J2 = -0.142675;
	const double c_J2 = -0.0822054;
	const double d_J2 = 0.0549608;

	for (unsigned long int ind = 0; ind < stop_ind; ind++)
	{
		double sigma = (double)ind / JLUT_MULTIPLIER;
		JLUT[ind] = a_J1 * pow(sigma, 3) + b_J1 * pow(sigma, 2) + c_J1 * sigma;
	}

	for (unsigned long int ind = stop_ind; ind < numels; ind++)
	{
		double sigma = (double)ind / JLUT_MULTIPLIER;
		JLUT[ind] = 1 - exp(a_J2 * pow(sigma, 3) + b_J2 * pow(sigma, 2) + c_J2 * sigma + d_J2);
	}
}

double getJ(const double &sigma, const std::vector<double> &JLUT)
{
	if (sigma >= 10.0)
	{
		return 1.0;
	}
	else
	{
		unsigned long int ind = sigma * JLUT_MULTIPLIER;

		return JLUT[ind];
	}
}

double getJ(const double &sigma)
{
	if ((sigma >= 0) && (sigma <= SIGMA_STAR))
	{
		return A_J1 * pow(sigma, 3) + B_J1 * pow(sigma, 2) + C_J1 * sigma;
	}
	else if (sigma < 10)
	{
		return 1.0 - exp(A_J2 * pow(sigma, 3) + B_J2 * pow(sigma, 2) + C_J2 * sigma + D_J2);
	}
	else if (sigma >= 10)
	{
		return 1.0;
	}
	else
	{
		return -1.0;
	}
}

void createJinvLUT(std::vector<double> &JinvLUT)
{
	unsigned int numels = JINVLUT_MULTIPLIER;
	JinvLUT.resize(numels);

	unsigned int stop_ind = I_STAR * JINVLUT_MULTIPLIER;

	const double a_s1 = 1.09542;
	const double b_s1 = 0.214217;
	const double c_s1 = 2.33727;

	const double a_s2 = 0.706692;
	const double b_s2 = 0.386013;
	const double c_s2 = -1.75017;

	for (unsigned int ind = 0; ind < stop_ind; ind++)
	{
		double I = (double)ind / JINVLUT_MULTIPLIER;
		JinvLUT[ind] = a_s1 * pow(I, 2) + b_s1 * I + c_s1 * pow(I, 2);
	}

	for (unsigned int ind = stop_ind; ind < numels; ind++)
	{
		double I = (double)ind / JINVLUT_MULTIPLIER;
		JinvLUT[ind] = -a_s2 * log(b_s2 * (1 - I)) - c_s2 * I;
	}
}

double getJinv(const double &I, const std::vector<double> &JinvLUT)
{
	unsigned long int ind = I * JINVLUT_MULTIPLIER;
	return JLUT[ind];
}

double getJinv(const double &I)
{
	if ((I >= 0) && (I <= I_STAR))
	{
		return (A_S1 * pow(I, 2)) + (B_S1 * I) + (C_S1 * sqrt(I));
	}
	else if (I < 1)
	{
		return (-A_S2 * log(B_S2 * (1 - I))) - (C_S2 * I);
	}
	else
	{
		return -1.0;
	}
}