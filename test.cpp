#include <iostream>
#include </usr/include/eigen3/Eigen/Dense>
#include "ProtographClass.hpp"

using namespace std;

int main(void)
{
	Protograph P;

	Eigen::MatrixXi B(3, 6);
	Eigen::MatrixXi H;
	B << 1 , 2 , 2 , 1 , 1 , 0 , 2 , 1 , 2 , 1 , 1 , 0 , 1 , 1 , 0 , 0 , 0 , 1;

	cout << "B = \n" << B << endl;

	P = Protograph(B);

	for (unsigned int i = 0; i < 1; i++)
	{
		try
		{
			P.lift(1667, 1, H);
			cout << "H.rows() = " << H.rows() << "\t\tH.cols() = " << H.cols() << endl;
		}
		catch (exception &e)
		{
			std::cout << "Exception: " << e.what() << std::endl << "i = " << i << std::endl;;
			exit(0);
		}
	}
	
	return 0;
}