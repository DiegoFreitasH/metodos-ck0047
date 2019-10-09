#ifndef EQUATION_HPP
#define EQUATION_HPP

class Equation{
    private:
        double coef_[4];

    public:
        Equation(double x0,double x1, double x2, double x3, double x4);
        double function(double x);
        double derivatedFunction(double x);
        void positiveRoots();
        void negativeRoots();

};

#endif