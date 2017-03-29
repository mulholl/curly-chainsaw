#include "ProtographClass.hpp"

using namespace std;
using namespace Eigen;

/* 	Default constructor */
Protograph::Protograph() : BM(MatrixXi(0,0)), numVNs(0), numCNs(0){
	symbolsToTransmit.clear();
	puncturedVNs.clear();

	R = 0;
}

/* 	Constructor
	Arguments:
		Input:
			const Eigen::MatrixXi &B 	-	The base matrix of the protograph
*/
Protograph::Protograph(const MatrixXi &B) : BM(B), numVNs(B.cols()), numCNs(B.rows()) {
	symbolsToTransmit.clear();
	puncturedVNs.clear();
	numPunctured = 0;
	calcRate();
	calcEdges();
}

/* 	Constructor
	Arguments:
		Input:
			const std::string &fn 	- 	The name of the file containing the protograph
*/
Protograph::Protograph(const string &fn)
{
	readFile(fn);
}

void Protograph::readFile(const string &fn)
{
	ifstream ifs(fn, ios::in);
	istringstream iss;
	string str;

	if (ifs.fail())
	{
		throw PGFileFailure(fn);
	}

	unsigned int cols_this_row = 0;
	unsigned int cols = 0;
	unsigned int rows = 0;
	unsigned int val;

	unsigned int CN = 0, VN = 0;

	while (getline(ifs, str))
	{
		cols_this_row = 0;

		iss.str(str);

		while (iss >> val)
		{
			cols_this_row++;
		}

		cols = (cols == 0) ? cols_this_row : cols;
		if (cols)
		{
			if (cols_this_row == cols) 
			{
				cols = cols_this_row;
				rows++;
			}
			else if (cols_this_row)
			{
				throw PGFileFailure(fn);
			}
		}

		iss.clear();
	}

	ifs.clear();
	ifs.seekg(0, ios::beg);

	numCNs = rows;
	numVNs = cols;
	BM.resize(numCNs, numVNs);

	while (ifs >> BM(CN, VN))
	{
		if (++VN == numVNs)
		{
			VN = 0;
			if (++CN == numCNs)
			{
				break;
			}
		}
	}

	ifs.close();
}

void Protograph::saveFile(const string &fn)
{
	ofstream ofs(fn, ios::out);

	for (unsigned int CN = 0; CN < numCNs; CN++)
	{
		for (unsigned int VN = 0; VN < numVNs-1; VN++)
		{
			ofs << BM(CN, VN) << " ";
		}
		if (numVNs)
		{
			ofs << BM(CN, numVNs-1) << "\n";
		}
	}

	ofs.flush();
	ofs.close();
}

/* 	Constructor
	Arguments:
		Input:
			const Eigen::MatrixXi &B 		-	The base matrix of the protograph
			const std::vector<bool> &T 		-	A vector containing one boolean value for 
												each VN of the protograph whose value is TRUE
												if the bit associated with that is transmitted
												or FALSE if that bit is punctured
*/
Protograph::Protograph(const MatrixXi &B, const vector<bool> &T) : BM(B), numVNs(B.cols()), numCNs(B.rows()){
	symbolsToTransmit = T;

	size_t max = symbolsToTransmit.size();
	numPunctured = 0;

	puncturedVNs.resize(max);

	for (unsigned int i = 0; i < max; i++)
	{
		if (symbolsToTransmit[i])
		{
			puncturedVNs[numPunctured++] = i;
		}
	}

	puncturedVNs.resize(numPunctured);

	calcRate();
	calcEdges();
}

/* 	Constructor
	Arguments:
		Input:
			const Eigen::MatrixXi &B 		-	The base matrix of the protograph
			const std::vector<unsigned int> &P 	-	A vector containing the indices of each
													of the VNs of the protograph that are
													punctured
*/
Protograph::Protograph(const MatrixXi &B, const vector<unsigned int> &P) : BM(B), numVNs(B.cols()), numCNs(B.rows()) {
	vector<bool> T(B.cols(), true);

	puncturedVNs = P;
	numPunctured = 0;

	for (vector<unsigned int>::const_iterator it = P.begin(); it != P.end(); ++it){
		T[*it] = false;
		numPunctured++;
	}

	symbolsToTransmit = T;

	calcRate();
	calcEdges();
}

/*	Specify which VNs of the protograph are punctured
	Arguments:
		Input:
			const std::vector<bool> &T 		-	A vector containing one boolean value for 
												each VN of the protograph whose value is TRUE
												if the bit associated with that is transmitted
												or FALSE if that bit is punctured
*/
void Protograph::setPuncturedVNs(const std::vector<bool> &T)
{
	symbolsToTransmit = T;

	size_t max = symbolsToTransmit.size();
	numPunctured = 0;

	puncturedVNs.resize(max);

	for (unsigned int i = 0; i < max; i++)
	{
		if (symbolsToTransmit[i])
		{
			puncturedVNs[numPunctured++] = i;
		}
	}

	puncturedVNs.resize(numPunctured);

	/* The rate may have changed, so recalculate it */
	calcRate();
}

/*	Specify which VNs of the protograph are punctured
	Arguments:
		Input:
			const std::vector<unsigned int> &P 	-	A vector containing the indices of each
													of the VNs of the protograph that are
													punctured
*/
void Protograph::setPuncturedVNs(const std::vector<unsigned int> &P)
{
	vector<bool> T(BM.cols(), true);

	puncturedVNs = P;
	numPunctured = 0;

	for (vector<unsigned int>::const_iterator it = P.begin(); it != P.end(); ++it){
		T[*it] = false;
		numPunctured++;
	}

	symbolsToTransmit = T;

	/* The rate may have changed, so recalculate it */
	calcRate();
}

/* Calculate the rate of the Protograph */
void Protograph::calcRate()
{
	unsigned int N = numVNs;
	unsigned int K = N - numCNs;

	R = (double) K / (N - numPunctured);
}

/* 	Return the rate of the protograph
	Return value:
		double 	-	The rate of the protograph
*/
double Protograph::Rate()
{
	return R;
}

/* Calculate the number of edges in the Protograph */
void Protograph::calcEdges()
{
	E = 0;

	for (unsigned int it1 = 0; it1 < numVNs; it1++)
	{
		for (unsigned int it2 = 0; it2 < numCNs; it2++)
		{
			E += BM(it2, it1);
		}
	}
}

/* 	Return the number of edges in the protograph
	Return value:
		unsigned int 	-	The number of edges in the protograph
*/
unsigned int Protograph::Edges()
{
	return E;
}

/* 	Return a referene to the base matrix of the protograph
	Return value:
		const MatrixXi& 	- Reference to the protograph base matrix
*/
const MatrixXi& Protograph::BaseMatrix() const
{
	return BM;
}

/* 	Return the number of VNs in the protograph
	Return value:
		unsigned int 	-	The number of VNs in the protograph
*/
unsigned int Protograph::VNs()
{
	return numVNs;
}


/* 	Return the number of CNs in the protograph
	Return value:
		unsigned int 	-	The number of CNs in the protograph
*/
unsigned int Protograph::CNs()
{
	return numCNs;
}

/* 	Generate an LDPC code by randomly lifting the protograph
	Arguments:
		Input:
			const unsigned int &numCopies 	-	The number of copies of the protograph
												to be lifted
			const unsigned int &seed		-	Seed to be used with the pseudo-RNG
		Output:
			Eigen::MatrixXi &H				-	The parity-check matrix of the generated 
												LDPC code
*/
void Protograph::lift(const unsigned int &numCopies, const unsigned int &seed, MatrixXi &H)
{
	unsigned int M = numCNs * numCopies;
	unsigned int N = numVNs * numCopies;

	H.resize(M, N);
	H.setZero();

	minstd_rand RNG (seed);

	vector< vector<unsigned int> > minConnectedCNs;
	vector<unsigned int> minConnectedCNAmts, minConnectedCNEdges;

	liftSetup(numCopies, minConnectedCNs, minConnectedCNAmts, minConnectedCNEdges);

	unsigned int EdgeType = 0;
	unsigned int init_EdgeType;
	// unsigned int H_VN;

	for (unsigned int PG_VN = 0; PG_VN < numVNs; PG_VN++)
	{
		init_EdgeType = EdgeType;

		for (unsigned int H_VN = PG_VN * numCopies; H_VN < ((PG_VN+1) * numCopies); H_VN++)
		{
			EdgeType = init_EdgeType;

			for (unsigned int PG_CN = 0; PG_CN < numCNs; PG_CN++)
			{
				// cout << "BM(" << PG_CN << ", " << PG_VN << ") = " << BM(PG_CN, PG_VN) << ", EdgeType = " << EdgeType << endl;
				for (unsigned int counter = BM(PG_CN, PG_VN); counter > 0; counter--)
				{
					AddEdge(H_VN, EdgeType, numCopies, RNG, minConnectedCNs, minConnectedCNAmts, minConnectedCNEdges, H);
				}
				if (BM(PG_CN, PG_VN))
				{
					EdgeType++;
				}
			}
		}
	}

	if (!ValidateLiftedH(H))
	{
		throw liftexc();
	}
}

/* 	Generate an LDPC code by randomly lifting the protograph - this version uses a time
	seed, rather than a seed specified by the calling function
	Arguments:
		Input:
			const unsigned int &numCopies 	-	The number of copies of the protograph
												to be lifted
		Output:
			Eigen::MatrixXi &H 				-	The parity-check matrix of the generated 
												LDPC code
*/
void Protograph::lift(const unsigned int &numCopies, MatrixXi &H)
{
	unsigned int seed = chrono::system_clock::now().time_since_epoch().count();

	lift(numCopies, seed, H);
}

/* Find the number of non-zero elements in the protograph base matrix
*/
unsigned int Protograph::NNZ()
{
	unsigned int nonZeros = 0;

	for (unsigned int CN = 0; CN < numCNs; CN++)
	{
		for (unsigned int VN = 0; VN < numVNs; VN++)
		{
			if (BM(CN, VN))
			{
				nonZeros++;
			}
		}
	}

	return nonZeros;
}

/*	Perform the setup necessary for the first edge to be added to a parity-check matrix
	from a lifted protograph.

	In the following, a "minimally-connected" CN for a given edge type is one which may be
	connected to edges of that type and which, among all such CNs, has the fewest edges of
	that type currently connected to it. For example, if four CNs are connected to 2, 2, 3
	and 4 edges of the same type, the CNs that have only two such connections are 
	"minimally-connected" by this definition.
	Arguments:
		Input:
			const unsigned int &numCopies 	-	The number of copies of the protograph
												to be lifted
			std::vector< std::vector<unsigned int> > &minConnectedCNs
											-	Each element of this vector is itself a
												vector containig the row-indices of each
												of the "minimally-connected" CNs for a
												given edge type.
			std::vector<unsigned int> &minConnectedCNAmts
											-	Each element of this vector represents
												the number of "minimally-connected" CNs
												for a given edge type, i.e., the size of
												each of the vectors in minConnectedCNs
			std::vector<unsigned int> &minConnectedCNEdges
											-	Each element of this vector represents
												the number of edges currently connecte
												to each of the "minimally-connected" CNs
												for a given edge type, i.e., the number 
												of edges connected to the CNs listed in
												each of the vectors in minConnectedCNs
												(every CN in each of these vectors is,
												by definition, connected to the same
												number of edges of the given type)
*/			
void Protograph::liftSetup(const unsigned int &numCopies, vector< vector<unsigned int> > &minConnectedCNs, vector<unsigned int> &minConnectedCNAmts, vector<unsigned int> &minConnectedCNEdges)
{
	unsigned int nonZeros = NNZ();
	minConnectedCNs.resize(nonZeros);

	minConnectedCNAmts.assign(nonZeros, numCopies);
	minConnectedCNEdges.assign(nonZeros, 0);

	for (vector< vector<unsigned int> >::iterator it_1 = minConnectedCNs.begin(); it_1 < minConnectedCNs.end(); it_1++)
	{
		(*it_1).resize(numCopies);
		
		CNsForEdgeType(distance(minConnectedCNs.begin(), it_1), numCopies, *it_1);
	}
}

/*	Add an edge to a parity-check matrix being constructed by lifting the protograph
	Arguments:
		Input:
			const unsigned int &VN  		-	The index of the VN
			const unsigned int &EdgeType 	-	A number representing the edge type
			const unsigned int &numCopies 	-	The number of copies of the protograph
												being made in the lifting process
			minstd_rand &RNG 				-	The pseudo-random number generator to
												be used in choosing which CN to connect
												the edge to
		Output
			std::vector< std::vector<unsigned int> > &minConnectedCNs
											-	Each element of this vector is itself a
												vector containig the row-indices of each
												of the "minimally-connected" CNs for a
												given edge type.			
			std::vector<unsigned int> &minConnectedCNAmts
											-	Each element of this vector represents
												the number of "minimally-connected" CNs
												for a given edge type, i.e., the size of
												each of the vectors in minConnectedCNs
			std::vector<unsigned int> &minConnectedCNEdges
											-	Each element of this vector represents
												the number of edges currently connecte
												to each of the "minimally-connected" CNs
												for a given edge type, i.e., the number 
												of edges connected to the CNs listed in
												each of the vectors in minConnectedCNs
												(every CN in each of these vectors is,
												by definition, connected to the same
												number of edges of the given type)
			Eigen::MatrixXi &H 				-	The parity-check matrix to which the edge
												is to be added
*/
void Protograph::AddEdge(const unsigned int &VN, const unsigned int &EdgeType, const unsigned int &numCopies, minstd_rand &RNG, vector< vector<unsigned int> > &minConnectedCNs, vector<unsigned int> &minConnectedCNAmts, vector<unsigned int> &minConnectedCNEdges, MatrixXi &H)
{
	const unsigned int max_tries = 10; /* The maximum number of attempts */
	unsigned int try_num; /* Used to keep count of the number of failed attempts */

	unsigned int max = minConnectedCNAmts[EdgeType]; // The number of minimally connected CNs for this edge type

	unsigned int minConnectedCNs_ind; /* Used to find the index in the parity-check matrix of the chosen CN */

	unsigned int temp; /* Used later for swapping */

	for (try_num = 0; try_num < max_tries; try_num++)
	{
		/* 	If max > 1, there is more than one possible CN to connect to, so we use the
			RNG to pick it
		*/
		minConnectedCNs_ind = (max > 1) ? (RNG() % (max - 1)) : 0;

		unsigned int CN = (minConnectedCNs[EdgeType])[minConnectedCNs_ind]; /* The index in the parity-check matrix of the chosen CN */

		if (!H(CN, VN))
		{
			/* First, connect the edge! */
			H(CN, VN) = 1;

			minConnectedCNAmts[EdgeType]--; // Reduce the count of the number of minimally connected CNs for this edge type

			/*  Check if we've used up the last of the minimally connected CNs for edges of this
				type. If so, we need to generate a new vector of minimally connected CNs for this
				edge type
			*/
			if (!minConnectedCNAmts[EdgeType])
			{
				minConnectedCNEdges[EdgeType]++;
				minConnectedCNAmts[EdgeType] = numCopies;

				/*  Generate the new list of minimally connected CNs - its just a list of
					CNs that can connect to edges of this type, as they should all
					currently be connected to the same number of edges
				*/
				CNsForEdgeType(EdgeType, numCopies, minConnectedCNs[EdgeType]);
			}
			/* Otherwise, swap the index for the CN we just connected with the last one in the array

			The index for the CN we just connected could easily be discarded, but preserving it
			may be useful if debugging is needed... */
			else
			{
				temp = (minConnectedCNs[EdgeType])[minConnectedCNs_ind];
				(minConnectedCNs[EdgeType])[minConnectedCNs_ind] = (minConnectedCNs[EdgeType])[max - 1];
				(minConnectedCNs[EdgeType])[max - 1] = temp; 
			}

			/* All done, so exit the function */
			return;
		}
	}

	/* Check whether we reached the maximum number of attempts - if so, throw an exception */
	if (try_num == max_tries)
	{
		throw liftexc();
	}

}

/*	Generate a list of CNs that can be connected to a given edge type
	Arguments:
		Input:
			const unsigned int &EdgeType 	-	Number representing the edge type
			const unsigned int &numCopies 	-	The number of copies of the protograph
												to be lifted
		Output:
			std::vector<unsigned int> &ConnectableCNs 
											-	The row-indices of all CNs in the lifted
												parity-check matrix that can potentially
												be connected to an edge of the type
												specified by EdgeType
*/
void Protograph::CNsForEdgeType(const unsigned int &EdgeType, const unsigned int &numCopies, vector<unsigned int> &ConnectableCNs)
{
	ConnectableCNs.resize(numCopies);

	unsigned int nextEdgeType = 0;

	for (unsigned int VN = 0; VN < numVNs; VN++)
	{
		for (unsigned int CN = 0; CN < numCNs; CN++)
		{
			if (BM(CN, VN))
			{
				nextEdgeType++;
				if ((nextEdgeType - 1) == EdgeType)
				{
					for (unsigned int ind = 0; ind < numCopies; ind++)
					{
						ConnectableCNs[ind] = (CN * numCopies) + ind;
					}
				}
			}
		}
	}
}

/* 	Verify that a parity-check matrix H was correctly lifted from the base matrix
	Arguments:
		Input:
			cost Eigen::MatrixXi &H 	-	The parity-check matrix to be tested
	Output:
		bool 	-	TRUE if H was correctly lifted
					FALSE otherwise
*/
bool Protograph::ValidateLiftedH(const MatrixXi &H)
{
	unsigned int N = H.cols();
	unsigned int M = H.rows();

	unsigned int numCopies = N / numVNs;

	unsigned int VN, PG_CN;

	MatrixXi test(numCNs, numVNs);

	for (unsigned int copy = 0; copy < numCopies; copy++)
	{
		test.setZero();

		for (unsigned int PG_VN = 0; PG_VN < numVNs; PG_VN++)
		{
			VN = (PG_VN * numCopies) + copy;
			for (unsigned int CN = 0; CN < M; CN++)
			{
				PG_CN = CN / numCopies;
				if (H(CN, VN))
				{
					test(PG_CN, PG_VN)++;
				}
			}
		}
		if (test != BM)
		{
			return false;
		}
	}

	return true;
}