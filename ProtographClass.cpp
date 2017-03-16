#include "ProtographClass.hpp"

using namespace std;

/* 	Default constructor */
Protograph::Protograph() : BM(Eigen::MatrixXi(0,0)), numVNs(0), numCNs(0){
	symbolsToTransmit.clear();
	puncturedVNs.clear();

	R = 0;
}

/* 	Constructor
	Arguments:
		Input:
			const Eigen::MatrixXi &B 	-	The base matrix of the protograph
*/
Protograph::Protograph(const Eigen::MatrixXi &B) : BM(B), numVNs(B.cols()), numCNs(B.rows()) {
	symbolsToTransmit.clear();
	puncturedVNs.clear();
	numPunctured = 0;
	calcRate();
}

/* 	Constructor
	Arguments:
		Input:
			const Eigen::MatrixXi &B 	-	The base matrix of the protograph
			const vector<bool> &T 		-	A vector containing one boolean value for 
											each VN of the protograph whose value is TRUE
											if the bit associated with that is transmitted
											or FALSE if that bit is punctured
*/
Protograph::Protograph(const Eigen::MatrixXi &B, const vector<bool> &T) : BM(B), numVNs(B.cols()), numCNs(B.rows()){
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
}

/* 	Constructor
	Arguments:
		Input:
			const Eigen::MatrixXi &B 		-	The base matrix of the protograph
			const vector<unsigned int> &P 	-	A vector containing the indices of each
												of the VNs of the protograph that are
												punctured
*/
Protograph::Protograph(const Eigen::MatrixXi &B, const vector<unsigned int> &P) : BM(B), numVNs(B.cols()), numCNs(B.rows()) {
	vector<bool> T(B.cols(), true);

	puncturedVNs = P;
	numPunctured = 0;

	for (vector<unsigned int>::const_iterator it = P.begin(); it != P.end(); ++it){
		T[*it] = false;
		numPunctured++;
	}

	symbolsToTransmit = T;

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

/* 	Return a referene to the base matrix of the protograph
	Return value:
		const Eigen::MatrixXi& 	- Reference to the protograph base matrix
*/
const Eigen::MatrixXi& Protograph::BaseMatrix() const
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
			const unsigned int numCopies 	-	The number of copies of the protograph
												to be lifted
	Return value:
		Eigen::MatrixXi 	-	The parity-check matrix of the generated LDPC code
*/
Eigen::MatrixXi Protograph::lift(const unsigned int numCopies)
{
	unsigned int N = numVNs * numCopies;
	unsigned int M = numCNs * numCopies;

	std::vector<unsigned int> PGVNDegrees;
	std::vector<unsigned int> PGCNDegrees;

	std::vector<unsigned int> VNTypeStartInds(numVNs, 0);
	std::vector<unsigned int> VNTypeEndInds(numVNs, 0);
	std::vector<unsigned int> CNTypeStartInds(numCNs, 0);
	std::vector<unsigned int> CNTypeEndInds(numCNs, 0);

	/*	VNAvailableSockets (CNAvailableSockets) stores, for each VN (CN) of the
		lifted graph, the number of sockets that have not yet been connected
	*/
	std::vector<unsigned int> VNAvailableSockets;
	std::vector<unsigned int> CNAvailableSockets;

	/* 	VNNeededConnections (CNNeededConnections) stores, for each VN (CN) type, a 
		vector listing the types of each of the CNs (VNs) the VN (CN) must be 
		connected to 
	*/
	std::vector< std::vector<unsigned int> > VNNeededConnections;
	std::vector< std::vector<unsigned int> > CNNeededConnections;

	/* 	avalailableVNIndices and avalailableCNIndices store the indices of the VNs and CNs 
		of each type which have not yet had all of their sockets filled
		Inner vector i contains the indices of all of the nodes of type i that have not yet had all of their sockets filled
		The outer vector contains one inner vector for each of the node types
	*/
	std::vector< std::vector<unsigned int> > avalailableVNIndices;
	std::vector< std::vector<unsigned int> > avalailableCNIndices;

	findPGNodeDegrees(PGVNDegrees, PGCNDegrees);

	findVNTypeConnections(PGVNDegrees, VNNeededConnections);
	findCNTypeConnections(PGCNDegrees, CNNeededConnections);

	findStartEndInds(numVNs, numCopies, VNTypeStartInds, VNTypeEndInds);
	findStartEndInds(numCNs, numCopies, CNTypeStartInds, CNTypeEndInds);

	findAvailableVNIndices(VNTypeStartInds, VNTypeEndInds, numCopies, avalailableVNIndices);
	findAvailableCNIndices(CNTypeStartInds, CNTypeEndInds, numCopies, avalailableCNIndices);

	findLiftedDegrees(VNTypeStartInds, VNTypeEndInds, PGVNDegrees, VNAvailableSockets);
	findLiftedDegrees(CNTypeStartInds, CNTypeEndInds, PGCNDegrees, CNAvailableSockets);	

	// printVector(PGVNDegrees);
	// printVector(VNTypeStartInds);
	// printVector(VNTypeEndInds);
	// printVector(PGCNDegrees);
	// printVector(CNTypeStartInds);
	// printVector(CNTypeEndInds);

	std::cout << "\n\n\n" << std::endl;

	printVector(VNAvailableSockets);
	printVector(CNAvailableSockets);
	
	Eigen::MatrixXi H(M, N);

	return H;
}

/* 	Find the degrees of the VNs and CNs of the protograph
	Arguments:
		Output:
			std::vector<unsigned int> &PGVNDegrees	-	The VN degrees
			std::vector<unsigned int> &PGCNDegrees	-	The CN degrees
*/
void Protograph::findPGNodeDegrees(std::vector<unsigned int> &PGVNDegrees, std::vector<unsigned int> &PGCNDegrees)
{
	PGVNDegrees.resize(numVNs);
	PGCNDegrees.resize(numCNs);

	for (unsigned int i = 0; i < numVNs; i++)
	{
		for (unsigned int j = 0; j < numCNs; j++)
		{
			if (BM(j, i) > 0)
			{
				PGVNDegrees[i] += BM(j, i);
				PGCNDegrees[j] += BM(j, i);
			}
		}
	}
}

/* 	Find the CN types that each VN must be connected to 
	Arguments:
		Input:
			const std::vector<unsigned int> &PGVNDegrees 	-	Vector containing the degrees
																of each of the VNs in the
																protograph
		Output:
			std::vector< std::vector<unsigned int> > &VNNeededConnections
															-	Vector containing, for each VN
																type, a vector which lists the VN types the CN must be 
																connected to
*/
void Protograph::findVNTypeConnections(const std::vector<unsigned int> &PGVNDegrees, std::vector< std::vector<unsigned int> > &VNNeededConnections)
{
	VNNeededConnections.resize(numVNs);

	/* Find the VN types that each CN must be connected to */
	for (std::vector< std::vector<unsigned int> >::iterator it1 = VNNeededConnections.begin(); it1 < VNNeededConnections.end(); ++it1)
	{
		unsigned int j = std::distance(VNNeededConnections.begin(), it1);
		unsigned int n = *(PGVNDegrees.begin() + j);
		*it1 = std::vector<unsigned int>(n);

		std::vector<unsigned int>::iterator it2 = (*it1).begin();

		for (unsigned int i = 0; i < numCNs; i++)
		{
			unsigned int temp = BM(i, j);
			while (temp > 0)
			{
				*it2++ = i;
				temp--; 
			}
		}
		printVector(*it1);
	}
}

/* 	Find the VN types that each CN must be connected to 
	Arguments:
		Input:
			const std::vector<unsigned int> &PGCNDegrees 	-	Vector containing the degrees
																of each of the CNs in the
																protograph
		Output:
			std::vector< std::vector<unsigned int> > &CNNeededConnections
															-	Vector containing, for each VN
																type, a vector which lists the CN types the VN must be 
																connected to
*/
void Protograph::findCNTypeConnections(const std::vector<unsigned int> &PGCNDegrees, std::vector< std::vector<unsigned int> > &CNNeededConnections)
{
	CNNeededConnections.resize(numCNs);

	/* Find the VN types that each CN must be connected to */
	for (std::vector< std::vector<unsigned int> >::iterator it1 = CNNeededConnections.begin(); it1 < CNNeededConnections.end(); ++it1)
	{
		unsigned int i = std::distance(CNNeededConnections.begin(), it1);
		unsigned int n = *(PGCNDegrees.begin() + i);
		*it1 = std::vector<unsigned int>(n);

		std::vector<unsigned int>::iterator it2 = (*it1).begin();

		for (unsigned int j = 0; j < numVNs; j++)
		{
			unsigned int temp = BM(i, j);
			while (temp > 0)
			{
				*it2++ = j;
				temp--; 
			}
		}
		printVector(*it1);
	}
}

/*	Generate a vector containing the degrees of each of the VNs (CNs) of the code
	created by lifting the protograph
	Arguments:
		Input:
			const std::vector<unsigned int> &TypeStartInds	-	A vector containing the indices
																of the columns (rows) of the 
																lifted parity-check matrix
																corresponding to the first 
																occurrence of each VN (CN)
																type
			const std::vector<unsigned int> &TypeEndInds	-	A vector containing the indices
																of the columns (rows) of the 
																lifted parity-check matrix
																corresponding to the final 
																occurrence of each VN (CN)
																type
			const std::vector<unsigned int> &PGDegrees 		-	A vector containing the degrees
																of each of the VNs (CNs) in the
																protograph
		Output:
			std::vector<unsigned int> &LiftedDegrees 		-	A vector containing the degrees
																of each of the VNs (CNs) in the
																lifted graph
*/
void Protograph::findLiftedDegrees(const std::vector<unsigned int> &TypeStartInds, const std::vector<unsigned int> &TypeEndInds, const std::vector<unsigned int> &PGDegrees, std::vector<unsigned int> &LiftedDegrees)
{
	unsigned int max_ind = 0;

	for (std::vector<unsigned int>::const_iterator it = TypeEndInds.begin(); it < TypeEndInds.end(); it++)
	{
		max_ind = ((*it) > max_ind) ? (*it) : max_ind;
	}

	LiftedDegrees.resize(max_ind + 1);

	std::vector<unsigned int>::iterator start_it;
	std::vector<unsigned int>::iterator end_it;

	std::vector<unsigned int>::const_iterator it1 = TypeStartInds.begin();
	std::vector<unsigned int>::const_iterator it2 = TypeEndInds.begin();

	while (it1 < TypeStartInds.end())
	{
		start_it = LiftedDegrees.begin() + *it1;
		end_it = LiftedDegrees.begin() + *it2 + 1;

		unsigned int current_degree = *(PGDegrees.begin() + std::distance(TypeStartInds.begin(), it1));

		for (std::vector<unsigned int>::iterator it3 = start_it; it3 < end_it; it3++)
		{
			(*it3) = current_degree;
		}

		it1++;
		it2++;
	}

}

/* 	Generate a list the indices of the VNs of each type that have sockets available for 		
	connection (to be used before any sockets are connected)
	Arguments:
		Input:
			const std::vector<unsigned int> &VNTypeStartInds 	-	A vector containing the 
																	indices	of the rows of 
																	the  lifted parity-check
																	matrix corresponding 
																	to the first occurrence 
																	of each VN type
			const std::vector<unsigned int> &VNTypeEndInds 		-	A vector containing the 
																	indices	of the columns of 
																	the lifted parity-check
																	matrix corresponding 
																	to the final occurrence 
																	of each VN type
			const unsigned int &numCopies						-	The number of copies of the 														protograph to be lifted
		Output:
			std::vector< std::vector<unsigned int> > &avalailableVNIndices
																-	A vector containing, for
																	each VN type, a vector of
																	indices of VNs of that
																	type. Each index is
																	repeated X times, where
																	X is the number of sockets
																	available for connection
*/
void Protograph::findAvailableVNIndices(const std::vector<unsigned int> &VNTypeStartInds, const std::vector<unsigned int> &VNTypeEndInds, const unsigned int &numCopies, std::vector< std::vector<unsigned int> > &avalailableVNIndices)
{

	avalailableVNIndices.resize(numVNs);
	std::vector<unsigned int>::const_iterator it3 = VNTypeStartInds.begin();
	unsigned int temp = 0;

	for (std::vector< std::vector<unsigned int> >::iterator it1 = avalailableVNIndices.begin(); it1 < avalailableVNIndices.end(); it1++)
	{
		(*it1).resize(numCopies);
		temp = (*it3++);

		for (std::vector<unsigned int>::iterator it2 = (*it1).begin(); it2 < (*it1).end(); it2++)
		{
			(*it2) = temp++;
		}
		std::cout << "\t\t" << std::flush;
		printVector(*it1);
	}
}

/* 	Generate a list the indices of the CNs of each type that have sockets available for 		
	connection (to be used before any sockets are connected)
	Arguments:
		Input:
			const std::vector<unsigned int> &CNTypeStartInds 	-	A vector containing the 
																	indices	of the rows of 
																	the lifted parity-check
																	matrix corresponding 
																	to the first occurrence 
																	of each CN type
			const std::vector<unsigned int> &CNTypeEndInds 		-	A vector containing the 
																	indices	of the rows of 
																	the lifted parity-check
																	matrix corresponding 
																	to the final occurrence 
																	of each CN type
			const unsigned int &numCopies						-	The number of copies of the 														protograph to be lifted
		Output:
			std::vector< std::vector<unsigned int> > &avalailableCNIndices
																-	A vector containing, for
																	each CN type, a vector of
																	indices of CNs of that
																	type. Each index is
																	repeated X times, where
																	X is the number of sockets
																	available for connection
*/
void Protograph::findAvailableCNIndices(const std::vector<unsigned int> &CNTypeStartInds, const std::vector<unsigned int> &CNTypeEndInds, const unsigned int &numCopies, std::vector< std::vector<unsigned int> > &avalailableCNIndices)
{

	avalailableCNIndices.resize(numCNs);
	std::vector<unsigned int>::const_iterator it3 = CNTypeStartInds.begin();
	unsigned int temp = 0;

	for (std::vector< std::vector<unsigned int> >::iterator it1 = avalailableCNIndices.begin(); it1 < avalailableCNIndices.end(); it1++)
	{
		(*it1).resize(numCopies);
		temp = (*it3++);

		for (std::vector<unsigned int>::iterator it2 = (*it1).begin(); it2 < (*it1).end(); it2++)
		{
			(*it2) = temp++;
		}
		std::cout << "\t\t" << std::flush;
		printVector(*it1);
	}
}

/*	Find the start and end indices of blocks of nodes of each type in the lifted parity-check 
	matrix
	Arguments:
		Input:
			const unsigned int &numNodeTypes	-	The number of types of VN (CN), i.e., the
													number of columns (rows) in the protograph
													to be lifted
			const unsigned int &numCopies		-	The number of copies of the protograph to 
													be made in the lifting process
		Output:
			std::vector<unsigned int> &TypeStartInds	-	A vector containing the indices	of 													the columns (rows) of the lifted 
															parity-check matrix corresponding 
															to the first occurrence of each VN
															(CN) type
			std::vector<unsigned int> &TypeEndInds		-	A vector containing the indices	of 													the columns (rows) of the lifted 
															parity-check matrix corresponding 
															to the final occurrence of each VN
															(CN) type
*/
void Protograph::findStartEndInds(const unsigned int &numNodeTypes, const unsigned int &numCopies, std::vector<unsigned int> &TypeStartInds, std::vector<unsigned int> &TypeEndInds)
{
	unsigned int next_start_ind = 0;

	for (unsigned int i = 0; i < numNodeTypes; i++)
	{
		TypeStartInds[i] = next_start_ind;

		next_start_ind += numCopies;
		TypeEndInds[i] = next_start_ind - 1;
	}
}