#include "l_interval.hpp"
#include "l_real.hpp"
#include "l_imath.hpp"
#include "l_rmath.hpp"
#include <iostream>
#include <math.h>
#include <limits>


using namespace cxsc;
using namespace std;



l_interval cunif(l_interval a, l_interval b);
l_interval cpareto(double alfa, l_real cons, double x1, double x2);
l_interval cexpo(l_real alpha, l_interval x);
l_interval cnormal(l_real alpha, l_real mi, double x1, double x2, int n);
l_interval cgama(double lambda, double v, double x1, double x2);

int main()
{
	cout << SetDotPrecision(17,14);

    cout << "Uniforme Exemplo 1: " << cunif(l_interval(0.0,3.0),l_interval(1.0,2.0)) << "\n";
    cout << "Uniforme Exemplo 2: " << cunif(l_interval(0.0,3.0),l_interval(0.5,4.0/7.0)) << "\n";
	cout << "Pareto Exemplo 1: " << cpareto(0.25, l_real(1.0),1.0, 2.0) << "\n";
    cout << "Pareto Exemplo 2: " << cpareto(0.50, l_real(1.0),1.0, 2.0) << "\n";
    cout << "Pareto Exemplo 3: " << cpareto(0.75, l_real(1.0),1.0, 2.0) << "\n";
    cout << "Exponencial Exemplo 1: " << cexpo(l_real(0.01),l_interval(20.0,50.0)) << "\n";
    cout << "Normal Exemplo 1: " << cnormal(l_real(1.0),l_real(0.0),0.35,1.0,50) << "\n";
    cout << "Normal Exemplo 2: " << cnormal(l_real(1.0),l_real(0.0),-1.0,1.0,50) << "\n";
    cout << "Normal Exemplo 3: " << cnormal(l_real(1.0),l_real(0.0),0.5555,0.5555,50) << "\n";
    cout << "Normal Exemplo 4: " << cnormal(l_real(1.0),l_real(0.0),-0.5,0.25,50) << "\n";
    cout << "Gama: " << cgama(0.25,5.0,20.0,40.0) << "\n";
    cout << "Gama: " << cgama(0.81,7.81,0.0,10.0) << "\n";

}

l_interval cunif(l_interval a, l_interval b) {
	l_interval x;


	x = (Sup(a&b)-Inf(a&b))/(Sup(a)-(Inf(a)));

	return x;
}

l_interval cpareto(double alfa, l_real cons, double x1, double x2)
{

    l_real up = pow(cons,l_real(alfa));

    l_real min = l_real(-up*(std::pow(x1,-alfa)));
	l_real max = l_real(-up*(std::pow(x2,-alfa)));

	return l_interval(max-min,max-min);
}

l_interval cexpo(l_real alpha, l_interval x){

    l_real min = pow(l_real(M_E),((-alpha)*Sup(x)));
    l_real max = pow(l_real(M_E),((-alpha)*Inf(x)));

    return l_interval(max-min,max-min);
}

l_interval cnormal(l_real alpha, l_real mi, double x1, double x2, int n){

    l_real numeradorMin = pow(l_real(x1) - mi,l_real(2.0));
    l_real numeradorMax = pow(l_real(x2) - mi,l_real(2.0));

    l_real denominador = l_real(2.0)*pow(alpha,l_real(2.0));


    l_real aux = numeradorMin/denominador;
    l_real aux2 = numeradorMax/denominador;


    l_real firstNormal = (l_real(1.0)/(alpha*sqrt(2*M_PI)));
    l_real secondNormalMin = pow(l_real(M_E),-aux);
    l_real secondNormalMax = pow(l_real(M_E),-aux2);

    double h = (x1 - x2)/n;
    l_real min = l_real((firstNormal*secondNormalMin) + (firstNormal*secondNormalMax));
    l_real max = l_real((firstNormal*secondNormalMin) + (firstNormal*secondNormalMax));

    double value = x2;

    //Subdivisoes
    int i=0;
    while(i<(n-1)){

        value+=h;
        //Update values
        l_real numeradorMin = pow(l_real(value) - mi,l_real(2.0));
        l_real numeradorMax = pow(l_real(value) - mi,l_real(2.0));

        aux = numeradorMin/denominador;
        aux2 = numeradorMax/denominador;

        secondNormalMin = pow(l_real(M_E),-aux);
        secondNormalMax = pow(l_real(M_E),-aux2);
        //
        if(i%2 == 0){ //Par


            min+= 4*(l_real((firstNormal*secondNormalMin)));
            max+= 4*(l_real((firstNormal*secondNormalMax)));

        }
        else{
            min+= 2*l_real((firstNormal*secondNormalMin));
            max+= 2*l_real((firstNormal*secondNormalMax));

        }
        i++;
    }

    min = -min*(l_real(h)/l_real(3.0));
    max = -max*(l_real(h)/l_real(3.0));

    return l_interval(min,max);
}

l_interval cgama(double lambda, double v, double x1, double x2) {

    double firstMin = std::pow(M_E, -lambda * x1);
    double firstMax = std::pow(M_E, -lambda * x2);

    double min = firstMin * std::pow(x1, v - 1.0);
    double max = firstMax * std::pow(x2, v - 1.0);


    double gamaMin = min + max;
    double gamaMax = min + max;


    double n = 500;
    double h = (x1 - x2) / n;

    double value = x2;



    //Primeira parte do slide para correcao final
    //Necessita simpson
    int inf = 9000000;


    double aux1Inf = std::pow(M_E, -lambda * 0);
    double aux1Max = std::pow(M_E, -lambda * inf);


    double aux2Inf = aux1Inf * std::pow(0, v - 1.0);
    double aux2Max = aux1Max * std::pow(inf, v - 1.0);

    double aux3 = aux2Inf + aux2Max;


    double n2 = 90000000;
    double tam = (-inf)/n2;
    double vAux = inf;

    for(int k=1;k<n2-1;k++){
        vAux +=tam;
        aux1Inf = std::pow(M_E, -lambda * vAux);
        aux2Inf = aux1Inf * std::pow(vAux, v - 1.0);
        if(k%2==0)
            aux3+= 4*aux2Inf;
        else
            aux3+= 2*aux2Inf;
    }
    aux3 = -(tam/3)*aux3;

    int i = 0;
    //Gama
    while (i < n - 1) {
        value += h;

        firstMin = std::pow(M_E, -lambda * value);
        firstMax = std::pow(M_E, -lambda * value);

        min = firstMin * (std::pow(value, v - 1.0));
        max = firstMax * (std::pow(value, v - 1.0));

        if (i % 2 == 0) {
            gamaMin += 4 * min;
            gamaMax += 4 * max;

        }
        else {
            gamaMin += 2 * min;
            gamaMax += 2 * max;
        }
        i++;


    }


    gamaMin = -(h / 3.0) * gamaMin;
    gamaMax = -(h / 3.0) * gamaMax;

    return l_interval(gamaMin/aux3,gamaMax/aux3);

}
//(exp(-L*a)*a**(v-1))