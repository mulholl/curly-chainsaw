#ifndef PROTOGRAPH_CLASS_HPP
#define PROTOGRAPH_CLASS_HPP

#include <vector>
#include <utility>
#include </usr/include/eigen3/Eigen/Dense>
#include <iostream>
#include <random>
#include <chrono>
#include <exception>

class Protograph{
	public:
		/* Constructors */
		Protograph();
		Protograph(const Eigen::MatrixXi &);
		Protograph(const Eigen::MatrixXi &, const std::vector<bool> &);
		Protograph(const Eigen::MatrixXi &, const std::vector<unsigned int> &);

		/* Other functions */
		void lift(const unsigned int &, const unsigned int &, Eigen::MatrixXi &H);
		void lift(const unsigned int &, Eigen::MatrixXi &H);
	private:
		/* Private member variables */
		Eigen::MatrixXi BM; // The base matrix of the protograph
		unsigned int numVNs; // The number of VNs in the protograph
		unsigned int numCNs; // The number of CNs in the protograph
		std::vector<bool> symbolsToTransmit; /* Boolean value for each VN - TRUE if the 
												associated bit is transmitted, FALSE if 
												it is punctured */
		std::vector<unsigned int> puncturedVNs; // Indices of punctured VNs
		unsigned int numPunctured; // The number of punctured VNs
		double R; // The rate of the protograph
		unsigned int E; // The number of edges in the protograph

		/* Functions used by constructors to assign private member variables */
		void calcRate();
		void calcEdges();

		/* Functions that simply return private member variables */
		unsigned int VNs();
		unsigned int CNs();
		double Rate();
		unsigned int Edges();
		const Eigen::MatrixXi& BaseMatrix() const;

		/* Functions used by Protograph::lift() */
		unsigned int NNZ();
		void liftSetup(const unsigned int &, std::vector< std::vector<unsigned int> > &, std::vector<unsigned int> &, std::vector<unsigned int> &);
		void AddEdge(const unsigned int &, const unsigned int &, const unsigned int &, std::minstd_rand &, std::vector< std::vector<unsigned int> > &, std::vector<unsigned int> &, std::vector<unsigned int> &, Eigen::MatrixXi &);
		void CNsForEdgeType(const unsigned int &, const unsigned int &, std::vector<unsigned int> &);
		bool ValidateLiftedH(const Eigen::MatrixXi &);

};

class liftexc : public std::exception
{
  virtual const char* what() const throw()
  {
    return "Lifting failed!";
  }
};

/* Overloading the insertion operator for vectors of base types */
template <typename T>
static std::ostream& operator<<(std::ostream &os, const std::vector<T> &vec)
{

	if (vec.size())
	{
		os << *(vec.begin());
	}

	for (typename std::vector<T>::const_iterator i = vec.begin() + 1; i < vec.end(); i++)
	{
		os << " " << *i;
	}
	return os;
}

/* Overloading the insertion operator for vectors of pairs */
template <typename T1, typename T2>
static std::ostream& operator<<(std::ostream &os, const std::vector< std::pair<T1, T2> > &vec)
{

	if (vec.size())
	{
		os << "(" << (*(vec.begin())).first << ", " << (*(vec.begin())).second << ")";
	}

	for (typename std::vector< std::pair<T1, T2> >::const_iterator i = vec.begin() + 1; i < vec.end(); i++)
	{
		os << " (" << (*i).first << ", " << (*i).second << ")";
	}

	return os;
}

/* Overloading the insertion operator for vectors of vectors of base types */
template <typename T>
static std::ostream& operator<<(std::ostream &os, const std::vector< std::vector<T> > &vec)
{

	if (vec.size())
	{
		os << *(vec.begin());
	}

	for (typename std::vector< std::vector<T> >::const_iterator i = vec.begin() + 1; i < vec.end(); i++)
	{
		os << "\n" << *i;
	}
	return os;
}	

#endif