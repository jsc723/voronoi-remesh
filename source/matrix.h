#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
using namespace std;
class Matrix
{
public:
	Matrix(void);
	~Matrix(void);
	double mat[5][5];
	friend Matrix operator + (const Matrix&, const Matrix&);
};

#endif 

