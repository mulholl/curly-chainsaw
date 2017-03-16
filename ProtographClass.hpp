#ifndef PROTOGRAPH_CLASS_HPP
#define PROTOGRAPH_CLASS_HPP

#include <vector>
#include </usr/include/eigen3/Eigen/Dense>
#include <iostream>

class Protograph{
	public:
		Protograph();
		Protograph(const Eigen::MatrixXi &);
		Protograph(const Eigen::MatrixXi &, const std::vector<bool> &);
		Protograph(const Eigen::MatrixXi &, const std::vector<unsigned int> &);
		Eigen::MatrixXi lift(const unsigned int);
	private:
		Eigen::MatrixXi BM; // The base matrix of the protograph
		unsigned int numVNs; // The number of VNs in the protograph
		unsigned int numCNs; // The number of CNs in the protograph
		std::vector<bool> symbolsToTransmit; /* Boolean value for each VN - TRUE if the 
												associated bit is transmitted, FALSE if 
												it is punctured */
		std::vector<unsigned int> puncturedVNs; // Indices of punctured VNs
		unsigned int numPunctured; // The number of punctured VNs
		double R; // The rate of the protograph
		void calcRate();
		double Rate();
		const Eigen::MatrixXi& BaseMatrix() const;
		unsigned int VNs();
		unsigned int CNs();
		template < typename T >
		static void printVector(const std::vector<T> &vec)
		{

			for (typename std::vector<T>::const_iterator i = vec.begin(); i < vec.end(); i++)
			{
				std::cout << *i << " ";
			}
			std::cout << std::endl;
		}
		void findPGNodeDegrees(std::vector<unsigned int> &, std::vector<unsigned int> &);
		void findVNTypeConnections(const std::vector<unsigned int> &, std::vector< std::vector<unsigned int> > &);
		void findCNTypeConnections(const std::vector<unsigned int> &, std::vector< std::vector<unsigned int> > &);
		void findLiftedDegrees(const std::vector<unsigned int> &, const std::vector<unsigned int> &, const std::vector<unsigned int> &, std::vector<unsigned int> &);
		void findAvailableVNIndices(const std::vector<unsigned int> &, const std::vector<unsigned int> &, const unsigned int &, std::vector< std::vector<unsigned int> > &);	
		void findAvailableCNIndices(const std::vector<unsigned int> &, const std::vector<unsigned int> &, const unsigned int &, std::vector< std::vector<unsigned int> > &);		
		void findStartEndInds(const unsigned int &, const unsigned int &, std::vector<unsigned int> &, std::vector<unsigned int> &);
};

#endif