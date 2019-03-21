#ifndef _Matrix_h
#define _Matrix_h


#include <iostream>
#include <complex>
#include"mayday.h"


using std::cout;
using std::ostream;
using std::istream;
using std::endl;
using std::complex;

typedef complex<double> cmp_dbl;

class Matrix {


private:
	cmp_dbl* mat;
	int row;
	int col;
	bool assert;
public:

friend Matrix operator+(const Matrix&, const Matrix&);
friend Matrix operator-(const Matrix&, const Matrix&);
friend Matrix operator*(const Matrix&, const Matrix&);
friend Matrix inverse(const Matrix&);

	Matrix() {
		row = col = 0;
		assert = false;

	}

	Matrix(const int col, const int row) {
		this->row = row;
		this->col = col;
		this->mat = new cmp_dbl[row*col];
		assert = true;

	}

	Matrix(const cmp_dbl* mat, const int col, const int row) {
		this->row = row;
		this->col = col;
		this->mat = new cmp_dbl[row*col];

		for(int i = 0; i < row * col; ++i)
			this->mat[i] = mat[i];

		assert = true;
 	}

	Matrix(const cmp_dbl mat1, 
		const cmp_dbl mat2,
		const cmp_dbl mat3,
		const cmp_dbl mat4) {
		this->row = 2;
		this->col = 2;
		this->mat = new cmp_dbl[4];

		this->mat[0] = mat1;
		this->mat[1] = mat2;
		this->mat[2] = mat3;
		this->mat[3] = mat4;

		assert = true;
 	}


	Matrix(const Matrix& mm) {
		this->row = mm.row;
		this->col = mm.col;
		this->mat = new cmp_dbl[row*col];

		for(int i = 0; i < row * col; ++i)
			this->mat[i] = mm.mat[i];

		assert = true;
	}



	Matrix& operator=(const Matrix& mm) {
		if(this == &mm) return *this;
		if(this->assert) 
			delete this->mat;
		this->row = mm.row;
		this->col = mm.col;
		this->mat = new cmp_dbl[row * col];
		assert = true;

		for(int i = 0; i < row * col; ++i)
			this->mat[i] = mm.mat[i];

		return *this;

	}


	Matrix& operator+=(const Matrix& mm) {
		if(this->row != mm.row || this->col != mm.col || !assert) {
			cout << "this->row != mm.row || this->col != mm.col in Matrix& operator+=(const Matrix& mm) " << endl;
				mayDayAbort("error");
		}
		for(int i = 0; i < row * col; ++i)
			this->mat[i] += mm.mat[i];
		return *this;

	}

	Matrix& operator-=(const Matrix& mm) {
		if(this->row != mm.row || this->col != mm.col || !assert) {
			cout << "this->row != mm.row || this->col != mm.col in Matrix& operator-=(const Matrix& mm) " << endl;
				mayDayAbort("error");
		}
		for(int i = 0; i < row * col; ++i)
			this->mat[i] -= mm.mat[i];
		return *this;

	}


	Matrix& operator*=(const Matrix& mm) {
		cout << "this->row != mm.row || this->col != mm.col in Matrix& operator*=(const Matrix& mm) " << endl;
		mayDayAbort("error");
		return *this;

	}

	void copy(const cmp_dbl* mat) {

		if(assert) {
			for(int i = 0; i < row * col; ++i)
				this->mat[i] = mat[i];
		} else {
			cout << " not asserted in copy " << endl;
			mayDayAbort("error");
		}
 	}

	double modul2() {
		if(this->row == 1 || this->col == 1) {
			double temp = 0;

			for(int i = 0; i < row * col; ++i)
				temp += this->mat[i].real() * this->mat[i].real() + 
				        this->mat[i].imag() * this->mat[i].imag();

			return temp;
		} else {
			cout << "You can't find modul^2 from not a Vector " << endl;
			mayDayAbort("error");

		}

 	}


	cmp_dbl det() {
		if(this->row == 2 && this->col == 2) {
			return this->mat[0] * this->mat[3] - this->mat[1] * this->mat[2];

		} else {
			cout << "You can't find determinant from not a 2x2 Matrix " << endl;
			mayDayAbort("error");

		}

 	}



	cmp_dbl getComp(const int col, const int row) {
		if(this->col >= col + 1 && this->row >= row + 1)
			return this->mat[col + row * this->col];
		else {
			cout << "wrong getComp indexes " << endl;
			mayDayAbort("error");
			return 0;
		}
	}


	void setComp(const cmp_dbl val, const int col, const int row) {
		if(this->col >= col + 1 && this->row >= row + 1)
			this->mat[col + row * this->col] = val;
		else {
			cout << "wrong setComp indexes " << endl;
			mayDayAbort("error");
		}
	}



	void write() const {
		cout << "Matrix = ";
			for(int i = 0; i < row * col; ++i) {
				cout << this->mat[i] << " ";
			}

		cout << endl;

	}

	~Matrix() {
		delete [] mat;
		mat = NULL;
		assert = false;
	}


};


Matrix operator+(const Matrix& a, const Matrix& b) {

	Matrix temp = a;
	return temp += b;

}

Matrix operator-(const Matrix& a, const Matrix& b) {

	Matrix temp = a;
	return temp -= b;

}

Matrix operator*(const Matrix& a, const Matrix& b) {
		if(a.col != b.row) {
			cout << "a.col != b.row in Matrix operator*(const Matrix& a, const Matrix& b) " << endl;
			mayDayAbort("error");
		}

		cmp_dbl* temp = new cmp_dbl[a.col * b.col];
		
		int count = 0;

		
		for(int i = 0; i < a.row; ++i) {
			for(int j = 0; j < b.col; ++j) {
				temp[count] = 0;
				for(int k = 0; k < a.col; ++k) {
					temp[count] += a.mat[k + i * a.col] * b.mat[j + k * b.col];
				}
				++count;
			}
		}


		Matrix TMP(temp, b.col, a.row);
		delete [] temp;
		return TMP;// Matrix([], col,  row)


}


Matrix inverse(const Matrix& mm) {
	if(mm.col != mm.row || mm.row != 2) {
			cout << "Only square matrix can be inverted !" << endl;
			mayDayAbort("error");
	}

	cmp_dbl det = mm.mat[0] * mm.mat[3] - mm.mat[1] * mm.mat[2];

	cmp_dbl temp[4];
	temp[0] = mm.mat[3] / det;
	temp[1] = -mm.mat[1] / det;
	temp[2] = -mm.mat[2] / det;
	temp[3] = mm.mat[0] / det;

	return Matrix(temp, 2, 2);

}




#endif


