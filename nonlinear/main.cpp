#include <iostream>
#include <cmath>

using namespace std;

typedef double(*derrivative_function)(double);

double newtonMethod(double (*f)(double),
                        double (*derrivative_f)(double),
                        double epsilon,
                        double start_point,
                        int stop_method,
                        double a,
                        double b){
    double x = start_point;

    int number_of_iterations = 0;

    double new_x = x - (f(x) / derrivative_f(x));
    switch(stop_method) {
        case 1:
        while(abs(f(new_x)) > epsilon) {
            x = new_x;
            new_x = x - (f(x) / derrivative_f(x));

            number_of_iterations++;

            if(new_x < a || new_x > b){
                cout << ";-1";
                return new_x;
            }
        }
        break;
        case 2:
        while(abs(new_x - x) > epsilon) {
            x = new_x;
            new_x = x - (f(x) / derrivative_f(x));

            number_of_iterations++;

            if(new_x < a || new_x > b){
                cout << ";-1";
                return new_x;
            }
        }
    }

    cout << ";" << number_of_iterations;

    return new_x;
}

double fun(double x){
    return 30.0*x*exp(-30.0) - 30.0 * exp(-11.0*x) + 1.0 / (11*30.0);
}

double der_fun(double x) {
    double d = (30.0*exp(-30.0) + 11.0 * 30.0 * exp(-11.0*x));
    if(d == 0){
        d = 1;
    }
    return d;
}

double secantMethod(double (*f)(double), 
                    double epsilon, 
                    double start, 
                    double end, 
                    double start_point1, 
                    double start_point2,
                    int stop_method){

    double x0 = start_point1, x1 = start_point2, new_x;
    
    int number_of_iterations = 0;;

    switch(stop_method){
        case 1:
        while(abs(f(new_x)) > epsilon) {
            double a = (f(x1) - f(x0)) / (x1 - x0);
            double b = f(x1);
            new_x = x1 - b/a;
            x0 = x1;
            x1 = new_x;

            number_of_iterations++;

            if(new_x < start || new_x > end){
                cout << ";-1";
                return new_x;
            }
        }
        break;
        case 2:
        while(abs(x1-x0) > epsilon) {
            double a = (f(x1) - f(x0)) / (x1 - x0);
            double b = f(x1);
            new_x = x1 - b/a;
            x0 = x1;
            x1 = new_x;

            number_of_iterations++;

            if(new_x < start || new_x > end){
                cout << ";-1";
                return new_x;
            }
        }
        break;
    }

    cout << ";" << number_of_iterations;

    return new_x;
}


int main()
{
    //parameters
    int stop_method;        //1- f(x);  2- new_x - x
    double a;               //a
    double b;               //b
    double start_point1;
    double start_point2;    //for secant method only
    double epsilon;         //parameter for stop method

    cin >> stop_method >> epsilon;
    a = -1.25;
    b = 1.5;
    cin >> start_point1 >> start_point2;

    cout << start_point1 << ";" << start_point2;

    double x = newtonMethod(&fun, &der_fun, epsilon, start_point1, stop_method, a, b);

    cout << ";" << x;
    x = secantMethod(&fun, epsilon, a, b, start_point1, start_point2, stop_method);

    cout << ";" << x << endl;

    return 0;
}
