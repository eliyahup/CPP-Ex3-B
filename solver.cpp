//#pragma  once

#include "solver.hpp"
#include <iostream>
#include <complex>
#include "solver.hpp"
using namespace std;
using namespace solver;



const RealVariable solver::operator*(const RealVariable &r1,const double x) {
 if(r1._f == true){  //if number*x^2
        return RealVariable (r1._a * x,r1._b , r1._c);
    }
    return RealVariable (r1._a,r1._b * x,r1._c);
}

const RealVariable solver::operator*(const double x, const RealVariable &r2) {
    if(r2._f == true){   //if number*x^2
        return RealVariable (r2._a * x,r2._b , r2._c);
    }
    return RealVariable (r2._a,r2._b * x, r2._c);
}
const RealVariable solver::operator-(const RealVariable &r1, const double x) {
    return RealVariable (r1._a, r1._b, r1._c - x);
}
const RealVariable solver::operator-(const double x ,const RealVariable &r1) {
    return RealVariable (r1._a, r1._b, r1._c - x);
}
const RealVariable solver::operator-(const RealVariable& r1 ){
    return RealVariable (r1._a, r1._b-1, r1._c);
}
const RealVariable solver::operator-(const RealVariable &r1, const RealVariable &r2) {
    return RealVariable (r1._a - r2._a,r1._b - r2._b, r1._c - r2._c);
}

const RealVariable solver::operator==(const RealVariable &r1, double x) {

    return RealVariable (r1._a, r1._b, r1._c - x);
}
const RealVariable solver::operator==(double x, const RealVariable &r1) {

    return RealVariable (r1._a, r1._b, r1._c - x);
}
const RealVariable solver::operator==(const RealVariable &r1, const RealVariable &r2) {

    return RealVariable (r1._a - r2._a,r1._b - r2._b, r1._c - r2._c);
}

const RealVariable solver::operator^(const RealVariable &r1,const double x) {
    if(x==2){   // if x^2
        return RealVariable (r1._a + 1,r1._b-1,r1._c,true);
    }else if(x==1){//if x^1
        return RealVariable (r1._a ,r1._b ,r1._c);
    }
    else{ //if x>2
        __throw_invalid_argument("the power greater than 2");
    }
}

const RealVariable solver::operator+(const RealVariable &r1, const RealVariable &r2) {
    return RealVariable (r1._a + r2._a,r1._b + r2._b, r1._c + r2._c);
}

const RealVariable solver::operator+(double x, const RealVariable &r2) {
    return RealVariable (r2._a, r2._b,  r2._c+x);
}

const RealVariable solver::operator+(const RealVariable &r1,const double x) {
    return RealVariable (r1._a ,r1._b , r1._c + x );
}

const RealVariable solver::operator/(const RealVariable &r1, const double x) {
  if(x==0){  
      __throw_invalid_argument("can't divide by 0");
  }
    return RealVariable (r1._a, r1._b / x, r1._c);
}


// I used the idea in this code https://www.programiz.com/cpp-programming/examples/quadratic-roots 

double solver::solve(const RealVariable& x){
    if(x._a==0){
        if(x._b==0){  //number == number
            __throw_invalid_argument("There is no solution");
        }
        return -x._c/x._b;
    }
    else{
        double _x1;
        double _x2;
        double discriminant = x._b*x._b - 4*x._a*x._c;
        if (discriminant > 0) {
            _x1 = (-x._b + sqrt(discriminant)) / (2*x._a);
            //_x2 = (-x._b - sqrt(discriminant)) / (2*x._a);
            return _x1;
        }
        else if (discriminant == 0) {
            _x1 = (-x._b + sqrt(discriminant)) / (2*x._a);
            return _x1;
        }
    } //if discriminant < 0
    __throw_invalid_argument("There is no real solution");
    return 0;
}

// I used the idea in this code https://www.programiz.com/cpp-programming/examples/quadratic-roots 

complex<double> solver::solve(const ComplexVariable& y){
    if(y._a==0){
        complex<double> ans;
        double realPart = -y._c/y._b;
        double imaginaryPart = -y._i/y._b;
        if(realPart==0){
            realPart=0;
        }
        if(imaginaryPart==0){
            imaginaryPart=0;
        }
        ans.real(realPart);
        ans.imag(imaginaryPart);
        return ans;
    }
    else{
        double _x1;
        double realPart;
        double imaginaryPart;
        double discriminant = y._b*y._b - 4*y._a*y._c;
        if (discriminant > 0) {
            _x1 = (-y._b + sqrt(discriminant)) / (2*y._a);
            //_x2 = (-y._b - sqrt(discriminant)) / (2*y._a);
            return _x1;
        }
        else if (discriminant == 0) {
            _x1 = (-y._b + sqrt(discriminant)) / (2*y._a);
            return _x1;
        }
        else { // if is ComplexVariable
            realPart = -y._b/(2*y._a);
            if(realPart==0){
                realPart = 0;
            }
            imaginaryPart =sqrt(-discriminant)/(2*y._a);
            complex<double> ans;
            ans.real(realPart);
            ans.imag(imaginaryPart);
            return ans;
        }
    }

}


const ComplexVariable solver::operator*(const double x, const ComplexVariable &r2) {
     if(r2._f == true){
        return ComplexVariable (r2._a * x,r2._b , r2._c);
    }
    return ComplexVariable(r2._a,r2._b * x, r2._c);
}

const ComplexVariable solver::operator*(const ComplexVariable &r1, double x) {
    if(r1._f == true){
        return ComplexVariable (r1._a * x,r1._b , r1._c);
    }
    return ComplexVariable(r1._a,r1._b * x,r1._c);
}

const ComplexVariable solver::operator-(const ComplexVariable &r1,const double x) {
    return ComplexVariable(r1._a, r1._b, r1._c - x);
}

const ComplexVariable solver::operator-(const double x ,const ComplexVariable &r1) {
    return ComplexVariable (r1._a, r1._b, r1._c - x);
}

const ComplexVariable solver::operator-(const ComplexVariable &r1, const ComplexVariable &r2) {
    return ComplexVariable (r1._a - r2._a,r1._b - r2._b, r1._c - r2._c);
}

const ComplexVariable solver::operator-( const ComplexVariable& r1, const complex<double> &r2) {
    return ComplexVariable(r1._a ,r1._b , r1._c, r1._i-r2.imag());
}
const ComplexVariable solver::operator-(const complex<double> &r2 , const ComplexVariable& r1) {
    return ComplexVariable(r1._a ,r1._b , r1._c, r1._i-r2.imag());
}

const ComplexVariable solver::operator==(const ComplexVariable &r1, const double x) {
    return ComplexVariable(r1._a, r1._b, r1._c - x);
}

const ComplexVariable solver::operator==(const ComplexVariable &r1, const ComplexVariable &r2) {
    return ComplexVariable(r1._a - r2._a,r1._b - r2._b, r1._c - r2._c,r1._i-r2._i);
}

const ComplexVariable solver::operator^(const ComplexVariable &r1, const double x) {
     if(x==2){
        return ComplexVariable (r1._a + 1,r1._b-1,r1._c,0.0,true);
    }else if(x==1){
        return ComplexVariable (r1._a ,r1._b ,r1._c);
    }
    else{
        __throw_invalid_argument("the power greater than 2");
    }
}

const ComplexVariable solver::operator+(const ComplexVariable &r1, const ComplexVariable &r2) {
    return ComplexVariable(r1._a + r2._a,r1._b + r2._b, r1._c + r2._c);
}

const ComplexVariable solver::operator+( const ComplexVariable& r1, const complex<double> &r2) {
    return ComplexVariable(r1._a ,r1._b , r1._c, r1._i+r2.imag());
}
const ComplexVariable solver::operator+(const complex<double> &r2 , const ComplexVariable& r1) {
    return ComplexVariable(r1._a ,r1._b , r1._c, r1._i+r2.imag());
}
const ComplexVariable solver::operator+(const double x, const ComplexVariable &r2) {
    return ComplexVariable(r2._a, r2._b,  r2._c+x);
}

const ComplexVariable solver::operator+(const ComplexVariable &r1, double x) {
    return ComplexVariable(r1._a ,r1._b , r1._c + x );
}

const ComplexVariable solver::operator/(const ComplexVariable &r1, const double x) {
     if(x==0){
      __throw_invalid_argument("can't divide by 0");
  }
    return ComplexVariable(r1._a, r1._b / x, r1._c);
}


