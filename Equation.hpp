#ifndef EQUATION_HPP
#define EQUATION_HPP

class Equation{
    private:
        double coef_[5], error_;
        int mult_;
    public:
        Equation();
        Equation(double x0,double x1, double x2, double x3, double x4);
        Equation(double x0,double x1, double x2, double x3, double x4, double error, int mult);
        void setError(double err);
        double getError();
        void setMult(int mult);
        int getMult();
        double function(double x);
        double derivatedFunction(double x);
        void positiveRoots();
        double localizeRoots();
        double isolation();
        void print();
        //friend std::ostream& operator<<(std::ostream& out, const Equation& eq);

};

#endif