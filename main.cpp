#include "Equation.hpp"
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <utility>

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

pair<int, double> newtonRaphson(Equation eq, double x, double p, double error, int ordem){
    int k = 1;
    double xk, xk1;
    if(x == -1){
        return pair<int, double>(-1, -1);
    }
    ofstream arquivoSaida("iteracoesNewtonRaphson-" + to_string(ordem) + ".csv");
    xk = x;
    if(abs(eq.function(xk)) < error) return pair<int, double>(0, xk);
    arquivoSaida << "Iteração,x,f(x)" << endl;
    while(true) {
        xk1 = xk - (p * eq.function(xk) / eq.derivatedFunction(xk));
        arquivoSaida << k << "," << xk1 << "," << eq.function(xk1) << endl;
        if (abs(eq.function(xk1)) < error || abs(xk1 - xk) < error) {
            arquivoSaida.close();
            return pair<int, double>(k, xk1);
        }
        xk = xk1;
        k++;
    }
}

pair<int, double> secante(Equation eq, double x, double x1, double p, double error, int ordem){
    int k = 1;
    double xk, xk1, xk2;
    if(x == -1){
        return pair<int, double>(-1, -1);
    }
    ofstream arquivoSaida("iteracoesSecante-" + to_string(ordem) + ".csv");
    xk = x;
    xk1 = x1;
    if(abs(eq.function(xk)) < error) return pair<int, double>(0, xk);
    if(abs(eq.function(xk1)) < error || abs(xk1 - xk) < error) return pair<int, double>(0, xk1);
    arquivoSaida << "Iteração,x,f(x)" << endl;
    while(true){
        xk2 = xk1 - (p*eq.function(xk1)*(xk1 - xk))/(eq.function(xk1) - eq.function(xk));
        arquivoSaida << k << "," << xk2 << "," << eq.function(xk2) << endl;
        if(abs(eq.function(xk2)) < error || abs(xk2 - xk1) < error) {
            arquivoSaida.close();
            return pair<int, double>(k, xk2);
        }
        xk = xk1;
        xk1 = xk2; 
        k++;    
    }
}

pair<int, double> newtonPolinomios(Equation eq, double x, double error, int ordem){
    int k = 1;
    double xk = x, xf, deltaX;
    if(x == -1){
        return pair<int, double>(-1, -1);
    }
    ofstream arquivoSaida("iteracoesNewtonPolinomial-" + to_string(ordem) + ".csv");
    arquivoSaida << "Iteração,x,f(x)" << endl;
    while(true){
        xf = eq.function(xk);
        if(abs(xf) < error) {
            arquivoSaida.close();
            return pair<int, double>(k, xk);
        }
        deltaX = xf/(eq.derivatedFunction(xk));
        xk  = xk - deltaX;
        arquivoSaida << k << "," << xk << "," << eq.function(xk) << endl;
        if(abs(deltaX) < error) {
            arquivoSaida.close();
            return pair<int, double>(k, xk);
        }
        k++;
    }
}

int main(){
    int n, i, mult;
    double coef[4];
    double x0, error;
    pair<int, double> raiz_nr, raiz_sec, raiz_np;

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

    ofstream raizes ("raizes.csv");
    raizes << "Equação,Raiz Newton-Raphson,Raiz Secante,Raiz Newton Polinomial" << endl;
    ofstream comparativo("comparar.csv");
    comparativo << "Equação,x0,Iterações em Newton-Raphson,Iterações em Secante,Iteraçeõs em Newton Polinomial" << endl;
    for(i = 0 ; i < n ; i++){
        error = equations[i].getError();
        mult = equations[i].getMult();
        x0 = equations[i].isolation();

        raiz_nr = newtonRaphson(equations[i], x0, mult, error, i);
        raiz_sec = secante(equations[i], x0, x0 + 1, mult, error, i);
        raiz_np = newtonPolinomios(equations[i], x0, error, i);

        raizes << equations[i] << ",";
        if(raiz_nr.second != -1) raizes << raiz_nr.second << "," << raiz_sec.second << "," << raiz_np.second << endl;
        else raizes << "sem raizes reais positivas" << endl;
        comparativo << equations[i] << "," << x0 << ",";
        if(raiz_nr.first != -1) comparativo << raiz_nr.first << "," << raiz_sec.first << "," << raiz_np.first << endl; else comparativo << "sem raizes reais positivas a serem calculadas" << endl;
    }
    cout << "Tabelas geradas em .csv" << endl;
    
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