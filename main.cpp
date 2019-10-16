#include "Equation.hpp"
#include <cmath>
#include <fstream>
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
    if(x == -1){
        cout << "Sem raízes reais positivas" << endl;
        return -1;
    }
    ofstream arquivoSaida("newtonRaphson.csv");
    xk = x;
    if(abs(eq.function(xk)) < error) return xk;
    desenharCabecalho();
    arquivoSaida << "Iteração,x,f(x)" << endl;
    while(true) {
        xk1 = xk - (p * eq.function(xk) / eq.derivatedFunction(xk));
        cout << "|" << setw(doubleWidth) << k << "|" << setw(doubleWidth) << xk << "|" << setw(doubleWidth) << eq.function(xk) << "|" << endl;
        arquivoSaida << k << "," << xk1 << "," << eq.function(xk1) << endl;
        cin.get();
        if (abs(eq.function(xk1)) < error || abs(xk1 - xk) < error) {
            desenharDivisoria();
            arquivoSaida.close();
            return xk1;
        }
        xk = xk1;
        k++;
    }
}

double secante(Equation eq, double x, double x1, double p, double error){
    int k = 1;
    double xk, xk1, xk2;
    if(x == -1){
        cout << "Sem raízes reais positivas" << endl;
        return -1;
    }
    ofstream arquivoSaida("secante.csv");
    xk = x;
    xk1 = x1;
    if(abs(eq.function(xk)) < error) return xk;
    if(abs(eq.function(xk1)) < error || abs(xk1 - xk) < error) return xk1;
    desenharCabecalho();
    arquivoSaida << "Iteração,x,f(x)" << endl;
    while(true){
        xk2 = xk1 - (p*eq.function(xk1)*(xk1 - xk))/(eq.function(xk1) - eq.function(xk));
        cout << "|" << setw(doubleWidth) << k << "|" << setw(doubleWidth) << xk2 << "|" << setw(doubleWidth) << eq.function(xk2) << "|" << endl;
        arquivoSaida << k << "," << xk2 << "," << eq.function(xk2) << endl;
        if(abs(eq.function(xk2)) < error || abs(xk2 - xk1) < error) {
            desenharDivisoria();
            arquivoSaida.close();
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
    if(x == -1){
        cout << "Sem raízes reais positivas" << endl;
        return -1;
    }
    ofstream arquivoSaida("newtonPolinomial.csv");
    desenharCabecalho();
    arquivoSaida << "Iteração,x,f(x)" << endl;
    while(true){
        xf = eq.function(xk);
        if(abs(xf) < error) {
            desenharDivisoria();
            return xk;
        }
        deltaX = xf/(eq.derivatedFunction(xk));
        xk  = xk - deltaX;
        cout << "|" << setw(doubleWidth) << k << "|" << setw(doubleWidth) << xk << "|" << setw(doubleWidth) << eq.function(xk) << "|" << endl;
        arquivoSaida << k << "," << xk << "," << eq.function(xk) << endl;
        if(abs(deltaX) < error) {
            desenharDivisoria();
            arquivoSaida.close();
            return xk;
        }
        k++;
    }
}

int main(){
    int n, i, mult;
    double coef[4];
    double raiz_nr, raiz_sec, raiz_np, x0, error;

    cout << "Entre a quatidade de Equações: ";
    cin >> n;
    Equation equations[n];
    cout << "Entre a precisão das equações: ";
    cin >> error;
    for(i = 0 ; i < n ; i++){
        cout << "Entre os coeficientes separados por espaços(a0 a1 a2 a3 a4): ";
        cin >> coef[0] >> coef[1] >> coef[2] >> coef[3] >> coef[4];
        cout << "Entre a multiplicidade: ";
        cin >> mult;
        equations[i] = Equation(coef[0], coef[1], coef[2], coef[3], coef[4], error, mult);
    }

    for(i = 0 ; i < n ; i++){
        error = equations[i].getError();
        mult = equations[i].getMult();
        x0 = equations[i].isolation();
        cout << "Isolamento: "<< x0 << endl;


        cout << "Método Newton-Raphson: " << endl;
        raiz_nr = newtonRaphson(equations[i], x0, mult, error);

        cout << endl << "Método da Secante: " << endl;
        raiz_sec = secante(equations[i], x0, x0 + 1, mult, error);

        cout << endl << "Método Newton para Polinomios: " << endl;
        raiz_np = newtonPolinomios(equations[i], x0, error);

        cout << endl << "Raizes" << endl << "-----------" << endl;
        cout << "Newton-Raphson: " << raiz_nr << endl;
        cout << "Método da Secante: " << raiz_sec << endl;
        cout << "Método Newton Polinomios: " << raiz_np << endl;
        cout << "----------------------------------" << endl;
    }
    
    /*
    double raiz_nr, raiz_sec, raiz_np;
    Equation eq_newraph = Equation(-2, 0, 1, 0, 0);
    Equation eq_sec = Equation(3, -9, 0, 1, 0);
    Equation eq_newpoli = Equation(-2, -1, 2, 1, 0);
    double newraph_x0 = eq_newraph.isolation();
    double sec_x0 = eq_sec.isolation();
    double newpoli_x0 = eq_newpoli.isolation();

    cout << "Isolamento: " << endl;
    cout << "Newton-Raphson: " << newraph_x0 << endl;
    cout << "Secante: " << sec_x0 << endl;
    cout << "Newton Polinomios: "<< newpoli_x0 << endl;
    cout << "-----------------------" << endl << endl;

    cout << "Método Newton-Raphson: " << endl;
    raiz_nr = newtonRaphson(eq_newraph, newraph_x0, 1.0, 0.000000001);

    cout << endl << "Método da Secante: " << endl;
    raiz_sec = secante(eq_sec, sec_x0, sec_x0 + 1, 1, 0.0005);

    cout << endl << "Método Newton para Polinomios: " << endl;
    raiz_np = newtonPolinomios(eq_newpoli, newpoli_x0, 0.001);

    cout << endl << "Raizes" << endl << "-----------" << endl;
    cout << "Newton-Raphson: " << raiz_nr << endl;
    cout << "Método da Secante: " << raiz_sec << endl;
    cout << "Método Newton Polinomios: " << raiz_np << endl;
    cout << "Localizando as Raizes: " << endl;
    cout << eq_sec.localizeRoots() << endl; 
    */  

    return 0;
}