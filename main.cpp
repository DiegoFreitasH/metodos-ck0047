#include "Equation.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

#define doubleWidth numeric_limits<double>::digits

using namespace std;
/*
int numOfRealRoots(double arr[], int len){
    int i, n = 0;
    double neg;
    for(i = 0 ; i < len -1 ; i++){
        if((i+1) % 2 != 0) neg = -arr[i]; 
        if(arr[i]*arr[i+1] < 0 || neg*arr[i+1] < 0)
            n += 1;
    }
}
*/
void desenharDivisoria(){
    cout << "+";
    for (int i = 0; i < doubleWidth; ++i) cout << "-";
    cout << "+";
    for (int i = 0; i < doubleWidth; ++i) cout << "-";
    cout << "+";
    for (int i = 0; i < doubleWidth; ++i) cout << "-";
    cout << "+" << endl;
}
void desenharCabecalho(){
    desenharDivisoria();
    /* TODO: Entender por que o x não precisa do <= */
    cout << "|";
    for (int i = 0; i < (doubleWidth-8)/2; ++i) cout << " ";
    cout << "Iteração";
    for (int i = 0; i <= (doubleWidth-8)/2; ++i) cout << " ";
    cout << "|";
    for (int i = 0; i < (doubleWidth-1)/2; ++i) cout << " ";
    cout << "x";
    for (int i = 0; i < (doubleWidth-1)/2; ++i) cout << " ";
    cout << "|";
    for (int i = 0; i < (doubleWidth-4)/2; ++i) cout << " ";
    cout << "f(x)";
    for (int i = 0; i <= (doubleWidth-4)/2; ++i) cout << " ";
    cout << "|" << endl;
    desenharDivisoria();
};

double newtonRaphson(Equation eq, double x, double p, double error){
    int k = 1;
    double xk, xk1;
    xk = x;
    if(abs(eq.function(x)) < error) return x;
    desenharCabecalho();
    while(true) {
        xk1 = xk - (p * eq.function(xk)) / eq.derivatedFunction(xk);
        cout << "|" << setw(doubleWidth) << k << "|" << setw(doubleWidth) << xk << "|" << setw(doubleWidth) << eq.function(xk) << "|" << endl;
        if (abs(eq.function(xk1)) < error || abs(xk1 - xk) < error) {
            desenharDivisoria();
            return xk1;
        }
        xk = xk1;
        k++;
    }
}

double secante(Equation eq, double x, double x1, double p, double error){
    int k = 1;
    double xk, xk1, xk2;
    xk = x;
    xk1 = x1;
    if(abs(eq.function(xk)) < error) return xk;
    if(abs(eq.function(xk1)) < error || abs(xk1 - xk) < error) return xk1;
    desenharCabecalho();
    while(true){
        xk2 = xk1 - (p*eq.function(xk1)*(xk1 - xk))/(eq.function(xk1) - eq.function(xk));
        cout << "|" << setw(doubleWidth) << k << "|" << setw(doubleWidth) << xk2 << "|" << setw(doubleWidth) << eq.function(xk2) << "|" << endl;
        if(abs(eq.function(xk2)) < error || abs(xk2 - xk1) < error) {
            desenharDivisoria();
            return xk2;
        }
        xk = xk1;
        xk1 = xk2; 
        k++;    
    }
}

double newtonPolinomios(Equation eq, double x, double error){
    int k = 1;
    double xk = x, xf, deltaX;
    desenharCabecalho();
    while(true){
        xf = eq.function(xk);
        if(abs(xf) < error) {
            desenharDivisoria();
            return xk;
        }
        deltaX = xf/(eq.derivatedFunction(xk));
        xk  = xk - deltaX;
        cout << "|" << setw(doubleWidth) << k << "|" << setw(doubleWidth) << xk << "|" << setw(doubleWidth) << eq.function(xk) << "|" << endl;
        if(abs(deltaX) < error) {
            desenharDivisoria();
            return xk;
        }
        k++;
    }
}

int main(){
    double raiz_nr, raiz_sec, raiz_np;
    Equation eq_newraph = Equation(-2, 0, 1, 0, 0);
    Equation eq_sec = Equation(3, -9, 0, 1, 0);
    Equation eq_newpoli = Equation(-2, -1, 2, 1, 0);

    cout << "Método Newton-Raphson" << endl;
    raiz_nr = newtonRaphson(eq_newraph, 1.0, 1.0, 0.000000001);

    cout << endl << "Método da Secante: " << endl;
    raiz_sec = secante(eq_sec, 0, 1, 1, 0.0005);

    cout << endl << "Método Newton para Polinomios: " << endl;
    raiz_np = newtonPolinomios(eq_newpoli, 2, 0.001);

    cout << endl << "Raizes" << endl << "-----------" << endl;
    cout << "Newton-Raphson: " << raiz_nr << endl;
    cout << "Método da Secante: " << raiz_sec << endl;
    cout << "Método Newton Polinomios: " << raiz_np << endl;

    return 0;
}