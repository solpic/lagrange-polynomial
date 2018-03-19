/*
 * Program to first solve LaGrange polynomials given P(xi) = yi (i from 1 to n)
 * and then to solve polynomial given the case P(xi) = yi, P'(xi) = zi (i from 1 to n)
 * Includes simply polynomial math class.
 */
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class Polynomial {
	vector<double> coefficients;
	public:
	Polynomial(){}
	Polynomial(double x[], int len) {
		coefficients.resize(len);
		for(int i = 0; i<len; i++)
			coefficients[i] = x[i];
	}
	double& operator[](const int i) {
		return coefficients[i];
	}
	int size(){
		return coefficients.size();
	}
	void stretch(int n) {
		if(n>coefficients.size()) {
			coefficients.resize(n, 0);
		}
	}
	double eval(double x) {
		double result = 0;
		for(int i = 0; i<coefficients.size(); i++) {
			result += coefficients[i]*pow(x, i);
		}
		return result;
	}
	static Polynomial one() {
		double a[] = {1};
		return Polynomial(a, 1);
	}
	Polynomial operator+(Polynomial& p) {
		stretch(p.size());
		p.stretch(size());
		
		Polynomial sum;
		sum.stretch(size());
		for(int i = 0; i<size(); i++) 
			sum[i] = coefficients[i] + p[i];
		return sum;
	}
	Polynomial operator-(Polynomial& p) {
		stretch(p.size());
		p.stretch(size());
		
		Polynomial sum;
		sum.stretch(size());
		for(int i = 0; i<size(); i++) 
			sum[i] = coefficients[i] - p[i];
		return sum;
	}
	Polynomial operator*(double x) {
		Polynomial product;
		product.stretch(size());
		for(int i = 0; i<size(); i++)
			product[i] = coefficients[i]*x;
		return product;
	}
	Polynomial multiply_by_x(int n) {
		Polynomial result;
		result.stretch(n+coefficients.size());
		
		for(int i = 0; i<coefficients.size(); i++) {
			result[i+n] = coefficients[i];
		}
		
		return result;
	}
	Polynomial operator*(Polynomial& p) {
		Polynomial product;
		
		for(int i = 0; i<size(); i++) {
			Polynomial c = (p*coefficients[i]).multiply_by_x(i);
			product = product + c;
		}
		return product;
	}
	friend ostream &operator<<(ostream&, const Polynomial&);
	//Creates lagrange term of x sub i
	static Polynomial lagrange_helper(double x[], double y, int i, int n) {
		double h[] = {1};
		Polynomial result(h, 1);
		for(int j = 0; j<n; j++) {
			if(j!=i) {
				double c[] = {-x[j], 1};
				Polynomial t(c, 2);	//(x - xj)
				result = result*t;
			}
		}
		double q = result.eval(x[i]);
		result = result*(y/q);
		return result;
	}
	static Polynomial lagrange(double x[], double y[], int n) {
		Polynomial result;
		for(int i = 0; i<n; i++) {
			Polynomial t = lagrange_helper(x, y[i], i, n);
			result = result + t;
		}
		return result;
	}
	static void lagrange_test(Polynomial p, double x[], double y[], int n) {
		cout<<"Testing lagrange polynomial: "<<p<<endl;
		for(int i = 0; i<n; i++) {
			cout<<"P("<<x[i]<<") = "<<p.eval(x[i])<<", should be "<<y[i]<<endl;
		}
	}
	Polynomial derivative() {
		Polynomial result;
		result.stretch(coefficients.size()-1);
		for(int i = 1; i<coefficients.size(); i++) {
			result[i-1] = coefficients[i]*i;
		}
		return result;
	}
	static Polynomial gen_polynomial(double x[], double y[], double z[], int n) {
		Polynomial t = lagrange(x, y, n);
		// g is (x-x1)(x-x2)...(x-xn)
		Polynomial g = one();
		for(int i = 0; i<n; i++) {
			double h[] = {-x[i], 1};
			Polynomial c(h, 2);
			g = g*c;
		}
		Polynomial t_prime = t.derivative();
		Polynomial g_prime = g.derivative();
		
		//theta[i] = (zi-t'(xi))/g'(xi)
		double *theta = new double[n];
		for(int i = 0; i<n; i++) {
			theta[i] = (z[i]-t_prime.eval(x[i]))/g_prime.eval(x[i]);
		}
		
		Polynomial q = lagrange(x, theta, n);
		delete[] theta;
		Polynomial result = q*g;
		result = result + t;
		return result;
	}
#define accuracy 0.001
	static void polynomial_test(double x[], double y[], double z[], int n) {
		int errors = 0;
		Polynomial p = gen_polynomial(x, y, z, n);
		Polynomial p_prime = p.derivative();
		cout<<"Testing polynomial \""<<p<<"\":"<<endl;
		for(int i = 0; i<n; i++) {
			cout<<"x[i] = "<<x[i]<<", y[i] = "<<y[i]<<", z[i] = "<<z[i]<<endl;
			double px = p.eval(x[i]);
			cout<<"\t";
			if(fabs(px-y[i])<accuracy)
				cout<<"Correct: ";
			else{
				cout<<"Error: ";
				errors += 1;
			}
			cout<<"P(x[i]) = "<<px<<endl;
		
			double ppx = p_prime.eval(x[i]);
			cout<<"\t";
			if(fabs(ppx-z[i])<accuracy)
				cout<<"Correct: ";
			else{
				cout<<"Error: ";
				errors += 1;
			}
			cout<<"P'(x[i]) = "<<ppx<<endl;
		}
		cout<<errors<<" errors"<<endl;
	}
};

ostream &operator<<(ostream& strm, const Polynomial& obj){
	for(int i = obj.coefficients.size()-1; i>=0; i--) {
		if(i>0)
			strm<<obj.coefficients[i]<<"*x^"<<i<<" + ";
		else
			strm<<obj.coefficients[i];
	}
	if(obj.coefficients.size()==0) {
		strm<<"0";
	}
	return strm;
}

/*
 * Given xi, P(xi) and P'(xi), solve for P(x)
 */
int main() {
	double x[] = {1, 2, 3, 8, 10};
	double y[] = {1, 8, 27, -20, 100};
	double z[] = {0, 3, 2, -5, 1};
	
	Polynomial::polynomial_test(x, y, z, 5);
	return 0;
}
