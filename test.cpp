#include <iostream>
#include </usr/include/eigen3/Eigen/Dense>
#include "ProtographClass.hpp"

using namespace std;

int main(void)
{
	Protograph P;

	Eigen::MatrixXi B(3, 6);
	B << 1 , 2 , 2 , 1 , 1 , 0 , 2 , 1 , 2 , 1 , 1 , 0 , 1 , 1 , 0 , 0 , 0 , 1;

	cout << "B = \n" << B << endl;

	P = Protograph(B);

	P.lift(10);

	
	return 0;
}