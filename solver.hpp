

#pragma  once
#include <iostream>
#include <complex>
#include "solver.hpp"
using namespace std;
namespace solver{

    class RealVariable {
        double _a;      //Coefficients of quadratic equation ax^2+bx+c=0 
        double _b;
        double _c;
        bool _f = false;  //to difference between x to x^2 
         
         // Constructor work also as default Constructor
    public:
        RealVariable(const double& a= 0.0, const double& b= 1.0, const double& c= 0.0, const bool& f=false)  //b=1 for init x Coefficient
        : _a(a), _b(b), _c(c),_f(f){

        }
                   //friend operator function to have 2 sides.
        friend  const RealVariable operator*(const RealVariable &r1 , const double x);          
        friend  const RealVariable operator*(const double x , const RealVariable& r2);
        friend  const RealVariable operator*(const RealVariable& r1, const RealVariable& r2);
        friend  const RealVariable operator-(const RealVariable& r1 ,const double x);
        friend  const RealVariable operator-(const double x, const RealVariable& r1);
        friend  const RealVariable operator-(const RealVariable& r1 );
        friend  const RealVariable operator-(const RealVariable& r1 ,const RealVariable& r2);
        friend  const RealVariable operator==(const RealVariable& r1 ,const double x);
        friend  const RealVariable operator==(const double x ,const RealVariable& r1);
        friend  const RealVariable operator==(const RealVariable& r1 ,const RealVariable& r2);
        friend  const RealVariable operator^(const RealVariable& r1 ,const double x);
        friend  const RealVariable operator+(const RealVariable& r1 ,const RealVariable& r2);
        friend  const RealVariable operator+(const double x ,const RealVariable& r2);
        friend  const RealVariable operator+(const RealVariable& r1 ,const double x);
        friend  const RealVariable operator/(const RealVariable& r1 ,const double x);

       friend double solve( const RealVariable& x);     // To have access to the private fields

    };


    class ComplexVariable {
    private:
        double _a;          //Coefficients of quadratic equation ax^2+bx+c=0 
        double _b;
        double _c;
        double _i;         //imaginary part
        bool _f = false;   //to difference between x to x^2 

    public:       // Constructor work also as default Constructor
        ComplexVariable (const double& a= 0.0, const double& b= 1.0, const double& c= 0.0,    //b=1 for init x Coefficient
                         const double& i=0.0, const bool& f=false)
                : _a(a), _b(b), _c(c),_i(i),_f(f) {
        }

                      //friend operator function to have 2 sides.
        friend  const ComplexVariable operator*(double x , const ComplexVariable& r2);
        friend  const ComplexVariable operator*( const ComplexVariable& r1, const double x);
        friend  const ComplexVariable operator-(const ComplexVariable& r1, const double x);
        friend  const ComplexVariable operator-(const double x, const ComplexVariable& r1);
        friend  const ComplexVariable operator-(const ComplexVariable& r1, const ComplexVariable& r2);
        friend  const ComplexVariable operator-( const ComplexVariable& r1, const complex<double>& r2);
        friend  const ComplexVariable operator-( const complex<double>& r2, const ComplexVariable& r1);
        friend  const ComplexVariable operator==(const ComplexVariable& r1, const double x);
        friend  const ComplexVariable operator==(const ComplexVariable& r1, const ComplexVariable& r2);
        friend  const ComplexVariable operator^(const ComplexVariable& r1, const double x);
        friend  const ComplexVariable operator+(const ComplexVariable& r1, const ComplexVariable& r2);
        friend  const ComplexVariable operator+( const ComplexVariable& r1, const complex<double>& r2);
        friend  const ComplexVariable operator+( const complex<double>& r2, const ComplexVariable& r1);
        friend  const ComplexVariable operator+(const double x ,const ComplexVariable& r2);
        friend  const ComplexVariable operator+(const ComplexVariable& r1 ,const double x);
        friend  const ComplexVariable operator/(const ComplexVariable& r1 ,const double x);

        friend  complex<double> solve(const ComplexVariable& y);    // To have access to the private fields
    };

          //for using in main class
     double solve(const RealVariable& x);
     complex<double> solve(const ComplexVariable& y);

};



