#include "Equation.hpp"
#include <vector>


Equation::Equation(double x0, double x1, double x2, double x3, double x4){
    coef_[0] = x0;
    coef_[1] = x1;
    coef_[2] = x2;
    coef_[3] = x3;
    coef_[4] = x4;
}

double Equation::function(double x){
    return (((coef_[4]*x + coef_[3])*x + coef_[2])*x + coef_[1])*x + coef_[0];
}

double Equation::derivatedFunction(double x){
    return (((4*coef_[4]*x + 3*coef_[3])*x + 2*coef_[2])*x + coef_[1]);
}

void Equation::positiveRoots(){
    int i, k = 0, v = 0, p;
    for(i = 0 ; i < 3 ; i++){
        if(coef_[i]*coef_[i+1] < 0) v++;
    }
    while(2*k <= v){
        p = v - 2*k; //Número par positivo menor que v
        k++;
    }
}