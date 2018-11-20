#include <iostream>
#include <cmath>
#include <boost/math/tools/polynomial.hpp>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assert.hpp>

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <string>
#include <utility>

//[polynomial_arithmetic_1
/*`and some using statements are convenient:
*/

using std::string;
using std::exception;
using std::cout;
using std::abs;
using std::pair;

using namespace boost::math;
using namespace boost::math::tools; // for polynomial
using boost::lexical_cast;
using namespace std;

using namespace boost::math;
using namespace boost::math::tools; // for polynomial
using boost::lexical_cast;

polynomial<double> calculate_polynomial(double(*f)(double), int points_number, int nodes_type, int polynomial_type);
double fun(double);
template <typename T>
string formula_format(polynomial<T> const &a);
template <typename T>
string formula_format1(polynomial<T> const &a);
double CV(polynomial<double> p, double point){
    double value = 0;
    for(int i=p.size()-1;i>=0;i--){
        value = (value) * point +  p[i];
    }

    //cout << "D: Polynomial: " << formula_format(p) << " (" << point << ") = " << value << endl;

    return value;
}
int main()
{
    int points_number;
    int nodes_type;
    int polynomial_type;
    cin >> points_number >> nodes_type >> polynomial_type;
    double(*f)(double) = &fun;

    polynomial <double> r = calculate_polynomial(f, points_number, nodes_type, polynomial_type);

    cout << "plot [-1 : 1] " << formula_format(r) << endl;
    cout << "pause -1" <<endl;
    return 0;
}

double fun(double p){
    return cos(p);
}

polynomial<double> calculate_polynomial(double(*f)(double), int points_number, int nodes_type, int polynomial_type){

    vector<pair<double, double>> table(points_number);

    //choose one path
    double point, step;
    switch(nodes_type){
        case 1:             //equal distances
        step = 2.0 / (points_number+1.0);
        point = -1.0;
        for(int i = 0;i < points_number; i++ ){
            point += step;
            table[i] = make_pair(point, f(point));
        }
        break;
        case 2:
        for(int i=points_number-1;i>=0;i--){
            point = cos(M_PI/points_number*(i+0.5));
            table[points_number-i-1] = make_pair(point, f(point));
        }
        break;
    }

    #ifdef OBD2
    cout << endl;
    for(int i = 0;i < points_number;i++) {
        cout << "(" << table[i].first << "," << table[i].second << ")" << endl;
    }
    cout << endl;
    #endif // OBD2


    //count polynomial
    polynomial<double> result;
    polynomial<double> helper;
    switch(polynomial_type){
        case 1:                 //Newton
        result = polynomial<double>({table[0].second});
        helper = polynomial<double>({1});
        for(int i=1;i<points_number;i++){
            helper *= polynomial<double>({-table[i-1].first, 1});
            double c = (table[i].second - CV(result, table[i].first))/CV(helper, table[i].first);
            result += helper * c;
        }
        break;
        case 2:
        result = polynomial<double>({0});
        for(int i = 0;i < points_number;i++){
            helper = polynomial<double>({1});
            for(int j = 0;j < points_number;j++){
                if(i!=j) {
                    helper *= polynomial<double>({-table[j].first, 1});
                    helper /= (table[i].first - table[j].first);
                }
            }
            helper *= table[i].second;
            result += helper;
            #ifdef OBD2
            cout << "D2: " << formula_format1(result) << " helper " << formula_format1(helper) << endl;
            #endif
        }
        break;
    }
    return result;
}

template <typename T>
string sign_str(T const &x)
{
  return x < 0 ? "-" : "+";
}

template <typename T>
string inner_coefficient(T const &x)
{
  string result(" " + sign_str(x) + " ");
  if (abs(x) != T(1))
      result += lexical_cast<string>(abs(x));
  return result;
}

template <typename T>
string formula_format(polynomial<T> const &a)
{
  string result;
  if (a.size() == 0)
      result += lexical_cast<string>(T(0));
  else
  {
    // First one is a special case as it may need unary negate.
    unsigned i = a.size() - 1;
    if (a[i] < 0)
        result += "-";
    if (abs(a[i]) != T(1))
        result += lexical_cast<string>(abs(a[i]));

    if (i > 0)
    {
      result += "*x";
      if (i > 1)
      {
          result += "**" + lexical_cast<string>(i);
          i--;
          for (; i != 1; i--)
              if (a[i])
                result += inner_coefficient(a[i]) + "*x**" + lexical_cast<string>(i);

          if (a[i])
            result += inner_coefficient(a[i]) + "*x";
      }
      i--;

      if (a[i])
        result += " " + sign_str(a[i]) + " " + lexical_cast<string>(abs(a[i]));
    }
  }
  return result;
} // string formula_format(polynomial<T> const &a)
template <typename T>
string formula_format1(polynomial<T> const &a)
{
  string result;
  if (a.size() == 0)
      result += lexical_cast<string>(T(0));
  else
  {
    // First one is a special case as it may need unary negate.
    unsigned i = a.size() - 1;
    if (a[i] < 0)
        result += "-";
    if (abs(a[i]) != T(1))
        result += lexical_cast<string>(abs(a[i]));

    if (i > 0)
    {
      result += "x";
      if (i > 1)
      {
          result += "^" + lexical_cast<string>(i);
          i--;
          for (; i != 1; i--)
              if (a[i])
                result += inner_coefficient(a[i]) + "x^" + lexical_cast<string>(i);

          if (a[i])
            result += inner_coefficient(a[i]) + "x";
      }
      i--;

      if (a[i])
        result += " " + sign_str(a[i]) + " " + lexical_cast<string>(abs(a[i]));
    }
  }
  return result;
} // string formula_format(polynomial<T> const &a)
