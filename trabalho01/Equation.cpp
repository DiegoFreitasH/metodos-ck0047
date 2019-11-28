#include "Equation.hpp"
#include <vector>
#include <math.h>
#include <iostream>
#include <sstream>

using namespace std;

Equation::Equation(){
    for(int i = 0 ; i < 4 ; i++) coef_[i] = 0;
}

Equation::Equation(double x0, double x1, double x2, double x3, double x4){
    coef_[0] = x0;
    coef_[1] = x1;
    coef_[2] = x2;
    coef_[3] = x3;
    coef_[4] = x4;
}

Equation::Equation(double x0, double x1, double x2, double x3, double x4, double error, int mult){
    coef_[0] = x0;
    coef_[1] = x1;
    coef_[2] = x2;
    coef_[3] = x3;
    coef_[4] = x4;
    error_ = error;
    mult_ = mult;
}

void Equation::setError(double err){
    error_ = err;
}

double Equation::getError(){
    return error_;
}

void Equation::setMult(int mult){
    mult_ = mult;
}

int Equation::getMult(){
    return mult_;
}

//Função em si
double Equation::function(double x){
    return (((coef_[4]*x + coef_[3])*x + coef_[2])*x + coef_[1])*x + coef_[0];
}

//Função derivada
double Equation::derivatedFunction(double x){
    return (((4*coef_[4]*x + 3*coef_[3])*x + 2*coef_[2])*x + coef_[1]);
}

//Função para estimar quantidades de raizes positivas
void Equation::positiveRoots(){
    int v = 0, i;
    vector<double> coeficientesNaoNulos;

    for (i = 0; i <= 4; i++) {
        if (coef_[i] != 0) coeficientesNaoNulos.push_back(coef_[i]);
    }

    for (i = 0; i < coeficientesNaoNulos.size()-1; i++) {
        if (coeficientesNaoNulos[i] * coeficientesNaoNulos[i+1] < 0) v++;
    }

    for (i = 0; i <= v; i += 2) {
        cout << v << " - p = " << i << " -> p = " << (v - i) << endl;
    }
}

//Estimar o limite máximo estimado das raizes positivas
double Equation::localizeRoots() {
    int n;
    double maximum = 0;

    for (int i = 0; i <= 4; i++) {
        if (coef_[i] != 0) n = i;
    }

    for (int i = 0; i <= n-1; i++) {
        if ((abs(coef_[i])/abs(coef_[n])) > maximum)
            maximum = (abs(coef_[i])/abs(coef_[n]));
    }

    return 1 + maximum;
}

//Isolamento, retornando o limite superior
double Equation::isolation(){
    int i;
    double positiveBound = ceil(localizeRoots());
   
    for (i = 0; i < positiveBound; i++) {
        if (function(i) * function(i+1) <= 0) {
            return (positiveBound);
        }
    }

    return -1;
}

ostream &operator<<(ostream &out, const Equation &eq) {
    std::ostringstream outs;
    outs << eq.coef_[4] << "x^4 + " << eq.coef_[3] << "x^3 + " << eq.coef_[2] << "x^2 + " << eq.coef_[1] << "x + " << eq.coef_[0];
    return out << outs.str();
}



//Printar a reação
void Equation::print() {
    int n;
	
    for (int i = 0; i <= 4; i++)
        if (coef_[i] != 0) n = i;
        
    for (int i = 4; i >= 2; i--) {
        if (i == n) cout << coef_[i] << "*X^" << i;
        else if (coef_[i] > 0) cout << " + " << coef_[i] << "*X^" << i;
        else if (coef_[i] < 0) cout << " - " << abs(coef_[i]) << "*X^" << i;
    }
   
    if (coef_[1] > 0) cout << " + " << coef_[1] << "*X";
    else if (coef_[1] < 0) cout << " - " << abs(coef_[1]) << "*X";
   
    if (coef_[0] > 0) cout << " + " << coef_[0];
    else if (coef_[0] < 0) cout << " - " << abs(coef_[0]);

    cout << endl;

}