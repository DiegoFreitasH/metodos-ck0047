#include <iostream>
#include <vector>
#include <math.h>
#include <bits/stdc++.h>

using namespace std;

typedef vector<vector<double>> matriz;

/*
    TODO Raise Error quando pivo = 0
    TODO Refatorar para usar ponteiros ao inves de vector
    Obs: tempo de compilação está mais demorado do que eu imaginava

*/

void initPermutationVector(vector<int>& p){
    for(int i = 0 ; i < p.size() ; i++) p[i] = i;
}

void permute(vector<int>& p, matriz& matrix, const int& k, const int& r){
    int aux, i;
    double temp;
    swap(p[k], p[r]);
    for(i = 0 ; i < matrix.size() ; i++){
        swap(matrix[k][i], matrix[r][i]);
    }
}

vector<double> permute(vector<double>& D, vector<int>& p, int n){
    vector<double> out(n);
    for(int i = 0 ; i < n ; i++){
        out[i] = D[p[i]]; 
    }
    return out;
}

pair<double, int> escolher_pivo(const matriz& matrix, const int& k){
    double pv = matrix[k][k], pv_aux;
    int r = k; 
    unsigned int i;
    for( i = k+1 ; i < matrix.size() ; i++){
        pv_aux = matrix[i][k];
        if(abs(pv_aux) > pv){
            pv = pv_aux;
            r = i;
        }
    }
    return pair<double, int>(pv, r);
}

vector<double> substituicoesSucessivas(const matriz& matrix, const vector<double>& vector_b,const int& n){
    int i, j;
    double soma;
    
    vector<double> vector_out(n);
    vector_out[0] = vector_b[0]/matrix[0][0];

    for(i = 1 ; i < n ; i++){
        soma = 0;
        for(j = 0 ; j < i ; j++){
            soma += matrix[i][j] * vector_out[j];
        }
        vector_out[i] = (vector_b[i] - soma)/matrix[i][i];
    }

    return vector_out;       
}

vector<double> substituicoesRetroativas(const matriz& matrix, const vector<double>& vector_b,const int& n){
    double soma;
    int i,j;
    
    vector<double> vector_out(n);
    vector_out[n-1] = vector_b[n-1]/matrix[n-1][n-1];

    for(i = n-2 ; i > -1 ; i--){
        soma = 0;
        for(j = i+1 ; j < n ; j++){
            soma += matrix[i][j] * vector_out[j];
        }
        vector_out[i] = (vector_b[i] - soma)/matrix[i][i];
    }

    return vector_out;
}

vector<double> calcSystem(const matriz& matrix, const matriz& matrix_L, const matriz& matrix_U, const vector<double>& vector_D,const int& n){
    vector<double> vector_X, vector_Y;
    vector_Y = substituicoesSucessivas(matrix_L, vector_D, n);
    vector_X = substituicoesRetroativas(matrix_U, vector_Y, n);
    return vector_X;
}

pair<vector<double>, bool> factLU(const matriz& matrix, matriz& matrix_L, matriz& matrix_U, vector<double>& vector_D, const int& n){
    int i, j, k;
    double m;
    vector<double> vector_X;
    bool error_flag = false;
    
    pair<double, int> pivo;
    vector<int> p(n);
    vector<double> D;
    matrix_U = matrix;
    matrix_L = matriz(n, vector<double>(n));
    matrix_L[0][0] = 1;

    
    initPermutationVector(p);
    
    for(k = 0 ; k < n - 1 ; k++){
        pivo = escolher_pivo(matrix_U, k);
        
        if(pivo.first == 0) return pair<vector<double>, bool>(vector_X, true);
        
        if(pivo.second != k) permute(p, matrix_U, k, pivo.second);
        
        for(i = k + 1 ; i < n ; i++){
            m = -1 * (matrix_U[i][k] / matrix_U[k][k]);
            matrix_U[i][k] = 0;
            matrix_L[i][k] = -m;
            
            for(j = k + 1 ; j < n ; j++){
                matrix_U[i][j] = matrix_U[i][j] + (m * matrix_U[k][j]);
                matrix_L[i][j] = (i == j) ? 1 : 0; 
            }   
        }
    }

    D = permute(vector_D, p, n);
    vector_X = calcSystem(matrix, matrix_L, matrix_U, D, n);
    return pair<vector<double>, bool>(vector_X, error_flag);

}

pair<vector<double>, bool> factDoolitle(const matriz& matrix, matriz& matrix_L, matriz& matrix_U, vector<double>& vector_D, const int& n){
    int i, j, k;
    double soma;
    bool error_flag = false;

    vector<double> vector_X;
    matrix_U = matriz(n, vector<double>(n));
    matrix_L = matriz(n, vector<double>(n));

    for(i = 0 ; i < n ; i++){

        for(k = i ; k < n ; k++){
            soma = 0;
            for(j = 0 ; j < i ; j++){
                soma += (matrix_L[i][j]*matrix_U[j][k]);
            }
            matrix_U[i][k] = matrix[i][k] - soma;
            if(matrix_U[i][i] == 0) return pair<vector<double>, bool> (vector_X, true);
        }

        for(k = i ; k < n ; k++){
            if(i == k)
                matrix_L[i][i] = 1;
            else {
                soma = 0;
                for(j = 0 ; j < i ; j++ ){
                    soma += (matrix_L[k][j] * matrix_U[j][i]);
                }
                matrix_L[k][i] = (matrix[k][i] - soma) / matrix_U[i][i];
            }
        }
    }
    
    vector_X = calcSystem(matrix, matrix_L, matrix_U, vector_D, n);
    return pair<vector<double>, bool>(vector_X, error_flag);
}


bool checkResult(matriz matrix, vector<double> vector_X, vector<double> vector_D, int n){
    int i, j;
    double soma = 0;
    for(i = 0 ; i < n ; i++){
        soma = 0;
        for(j = 0 ; j < n ; j++){
            soma += vector_X[j]*matrix[i][j];
        }
        if(soma != vector_D[i]) return false;
    }
    return true;
}

void printMatrix(matriz& matrix){
    for(unsigned int i = 0 ; i < matrix.size() ; i++){
        for(int j = 0 ; j < matrix.size() ; j++){
            printf("%8.2f ", matrix[i][j]);
        }
        cout << endl;
    }
}

/*
Testes:
    LU com pivotação:
        A = {{3, 2, 4} , {1, 1, 2}, {4, 3, -2}};
        D = {1, 2, 3};
        
    Doolitle(Qualquer matriz sem pivo nulo):   
        A = {{1, -3, 2}, {-2, 8, -1}, {4, -6, 5}};
        D = {11, -15, 29};

    Entrada do Trabalho:
        A = {{20, 7, 9}, {7, 30, 8}, {9, 8, 30}}
        D = {16, 38, 38}
*/

int main(int argc, char** argv){
    // Declaração dos Dados
    int n = 3, i;
    pair<matriz, matriz> factor;
    pair<vector<double>, bool> result;
    matriz matrix_A = {{1, -3, 2}, {-2, 8, -1}, {4, -6, 5}};
    matriz matrix_L;
    matriz matrix_U;
    vector<double> vector_D = {11, -15, 29};
    vector<double> vector_X;
    bool error;
    
    // result.first -> Vetor solução
    // result.second -> boolean de erro
    result = factLU(matrix_A, matrix_L, matrix_U, vector_D, n);
    // result = factDoolitle(matrix_A, matrix_L, matrix_U, vector_D, n);
    error = result.second;
    if(error){
        cout << "Error" << endl;
        return 0;
    }
    vector_X = result.first;

    // Printa o resutado 
    for(i = 0 ; i < n ; i++) printf("c%d: %5.2f\n", i, vector_X[i]);   
    
    // Checa o resultado
    if (checkResult(matrix_A, vector_X, vector_D, n))
        cout << "Ok" << endl;
    else cout << "Error" << endl;
    
    return 0;
}