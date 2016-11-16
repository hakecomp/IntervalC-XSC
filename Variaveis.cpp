#include "l_interval.hpp"
#include "l_real.hpp"
#include "l_imath.hpp"
#include "l_rmath.hpp"
#include <iostream>
#include <math.h>
#include <limits>
#include<vector>
#include<time.h>


using namespace cxsc;
using namespace std;


// Utilizando Primitiva ou Simpson Convencional
l_interval cunif(l_interval a, l_interval b);
l_interval cpareto(double alfa, l_real cons, double x1, double x2);
l_interval cexpo(l_real alpha, l_interval x);
l_interval cnormal(l_real alpha, l_real mi, double x1, double x2, int n);
l_interval cgama(double lambda, double v, double x1, double x2);

//Utilizando Simpson Intervalar - CAPRANI
l_interval cgamaCaprani(double lambda, double v, double x1, double x2,int n);
l_interval cnormalCaprani(l_real alpha, l_real mi, double x1, double x2, int n);
l_interval cexpoCaprani(l_real alpha, l_interval x,int n);
l_interval cparetoCaprani(double alfa, l_real cons, double x1, double x2, int n);


//Forma Real
double cgamaReal(double lambda, double v, double x1, double x2);
double cnormalReal(double alpha, double mi, double x1, double x2);
double cparetoReal(double alfa, double cons, double x1, double x2,int n);
double cexpoReal(double alpha, double x1,double x2,int n);
double cunifReal(double x1,double x2);

//Funcao auxiliar
double firstStepGama(double lambda, double v);

//Funcao de Erros
void erroAbsoluto(double xReal, l_interval x);
void erroRelativo(double xReal, l_interval x);

int main()
{
    //Medir tempo

	cout << SetDotPrecision(8,14);

    clock_t inicio = clock();

    l_interval uni1 = cunif(l_interval(0.0,3.0),l_interval(1.0,2.0));
    cout << "Uniforme Exemplo 1: " << uni1 << "\n";
    clock_t fimUni1 = clock() - inicio;
    cout << "Tempo:" << double(fimUni1)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cunif(l_interval(0.0,3.0),l_interval(1.0,2.0))) << endl;
    double uni1Real = cunifReal(0.0,3.0);
    cout << "Erro Absoluto :";
    erroAbsoluto(uni1Real,uni1);
    cout << "Erro Relativo :";
    erroRelativo(uni1Real,uni1);
    cout << endl;

    inicio = clock();
    l_interval uni2 = cunif(l_interval(0.0,3.0),l_interval(0.5,4.0/7.0));
    cout << "Uniforme Exemplo 2: " << uni2 << "\n";
    clock_t fimUni2 = clock() - inicio;
    double uni2Real = cunifReal(0.5,4.0/7.0);
    cout << "Tempo:" << double(fimUni2)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cunif(l_interval(0.0,3.0),l_interval(0.5,4.0/7.0))) << endl;

    cout << "Erro Absoluto :";
    erroAbsoluto(uni2Real,uni2);
    cout << "Erro Relativo :";
    erroRelativo(uni2Real,uni2);
    cout << endl;

    inicio = clock();
    l_interval pareto1 = cpareto(0.25, l_real(1.0),1.0, 2.0);
    cout << "Pareto Exemplo 1: " << pareto1 << "\n";
    clock_t fimpar1 = clock() - inicio;
    double pareto1Real = cparetoReal(0.25, 1.0,1.0,2.0,50);
    cout << "Tempo:" << double(fimpar1)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cpareto(0.25, l_real(1.0),1.0, 2.0)) << endl;
    cout << "Erro Absoluto :";
    erroAbsoluto(pareto1Real,pareto1);
    cout << "Erro Relativo :";
    erroRelativo(pareto1Real,pareto1);
    cout << endl;

    inicio = clock();
    l_interval pareto2 = cparetoCaprani(0.25, l_real(1.0),1.0, 2.0,50);
    cout << "Pareto Caprani Exemplo 1: " << pareto2 << "\n";
    clock_t fimpar2 = clock() - inicio;
    cout << "Tempo:" << double(fimpar2)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cparetoCaprani(0.25, l_real(1.0),1.0, 2.0,50)) << endl;
    double pareto2Real = cparetoReal(0.25,1.0,1.0,2.0,50);
    cout << "Erro Absoluto :";
    erroAbsoluto(pareto2Real,pareto2);
    cout << "Erro Relativo :";
    erroRelativo(pareto2Real,pareto2);
    cout << endl;


    inicio = clock();
    l_interval pareto3 = cpareto(0.50, l_real(1.0),1.0, 2.0);
    cout << "Pareto Exemplo 2: " << pareto3 << "\n";
    clock_t fimpar3 = clock() - inicio;
    cout << "Tempo:" << double(fimpar3)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cpareto(0.50, l_real(1.0),1.0, 2.0)) << endl;
    cout << "Erro Absoluto :";
    double pareto3Real = cparetoReal(0.5,1.0,1.0,2.0,50);
    erroAbsoluto(pareto3Real,pareto3);
    cout << "Erro Relativo :";
    erroRelativo(pareto3Real,pareto3);
    cout << endl;


    inicio = clock();
    l_interval pareto4 = cparetoCaprani(0.50, l_real(1.0),1.0, 2.0,50);
    cout << "Pareto Caprani Exemplo 2: " << pareto4 << "\n";
    clock_t fimpar4 = clock() - inicio;
    cout << "Tempo:" << double(fimpar4)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cparetoCaprani(0.50, l_real(1.0),1.0, 2.0,50)) << endl;
    cout << "Erro Absoluto :";
    double pareto4Real = cparetoReal(0.5,1.0,1.0,2.0,50);
    erroAbsoluto(pareto4Real,pareto4);
    cout << "Erro Relativo :";
    erroRelativo(pareto4Real,pareto4);
    cout << endl;

    inicio = clock();
    l_interval pareto5 = cpareto(0.75, l_real(1.0),1.0, 2.0);
    cout << "Pareto Exemplo 3: " << pareto5 << "\n";
    clock_t fimpar5 = clock() - inicio;
    cout << "Tempo:" << double(fimpar5)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cpareto(0.75, l_real(1.0),1.0, 2.0)) << endl;
    cout << "Erro Absoluto :";
    double pareto5Real = cparetoReal(0.75,1.0,1.0,2.0,50);
    erroAbsoluto(pareto5Real,pareto5);
    cout << "Erro Relativo :";
    erroRelativo(pareto5Real,pareto5);
    cout << endl;

    inicio = clock();
    l_interval pareto6 = cparetoCaprani(0.75, l_real(1.0),1.0, 2.0,50);
    cout << "Pareto Caprani Exemplo 3: " << pareto6 << "\n";
    clock_t fimpar6 = clock() - inicio;
    cout << "Tempo:" << double(fimpar6)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cparetoCaprani(0.75, l_real(1.0),1.0, 2.0,50)) << endl;
    cout << "Erro Absoluto :";
    double pareto6Real = cparetoReal(0.75,1.0,1.0,2.0,50);
    erroAbsoluto(pareto6Real,pareto6);
    cout << "Erro Relativo :";
    erroRelativo(pareto6Real,pareto6);
    cout << endl;

    inicio = clock();
    l_interval exp6 = cexpo(l_real(0.01),l_interval(20.0,50.0));
    cout << "Exponencial Exemplo 1: " << exp6 << "\n";
    clock_t fimexp1 = clock() - inicio;
    cout << "Tempo:" << double(fimexp1)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cexpo(l_real(0.01),l_interval(20.0,50.0))) << endl;
    cout << "Erro Absoluto :";
    double exp6Real = cexpoReal(0.01,20.0,50.0,50);
    erroAbsoluto(exp6Real,exp6);
    cout << "Erro Relativo :";
    erroRelativo(exp6Real,exp6);
    cout << endl;

    inicio = clock();
    l_interval exp = cexpoCaprani(l_real(0.01),l_interval(20.0,50.0),50);
    cout << "Exponencial Caprani Exemplo 1: " << exp << "\n";
    clock_t fimexp2 = clock() - inicio;
    cout << "Tempo:" << double(fimexp2)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cexpoCaprani(l_real(0.01),l_interval(20.0,50.0),50)) << endl;
    cout << "Erro Absoluto :";
    double expReal = cexpoReal(0.01,20.0,50.0,50);
    erroAbsoluto(expReal,exp);
    cout << "Erro Relativo :";
    erroRelativo(expReal,exp);
    cout << endl;

    inicio = clock();
    l_interval normal = cnormal(l_real(1.0),l_real(0.0),0.35,1.0,50);
    cout << "Normal Exemplo 1: " << normal << "\n";
    clock_t fimn1 = clock() - inicio;
    cout << "Tempo:" << double(fimn1)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cnormal(l_real(1.0),l_real(0.0),0.35,1.0,50)) << endl;
    cout << "Erro Absoluto :";
    double normalReal = cnormalReal(1.0,0.0,0.35,1.0);
    erroAbsoluto(normalReal,normal);
    cout << "Erro Relativo :";
    erroRelativo(normalReal,normal);
    cout << endl;

    inicio = clock();
    l_interval normal2 = cnormalCaprani(l_real(1.0),l_real(0.0),0.35,1.0,50);
    cout << "Normal Caprani Exemplo 1: " << normal2 << "\n";
    clock_t fimn2 = clock() - inicio;
    cout << "Tempo:" << double(fimn2)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cnormalCaprani(l_real(1.0),l_real(0.0),0.35,1.0,50)) << endl;
    cout << "Erro Absoluto :";
    double normal2Real = cnormalReal(1.0,0.0,0.35,1.0);
    erroAbsoluto(normal2Real,normal2);
    cout << "Erro Relativo :";
    erroRelativo(normal2Real,normal2);
    cout << endl;

    inicio = clock();
    l_interval normal3 = cnormal(l_real(1.0),l_real(0.0),-1.0,1.0,50);
    cout << "Normal Exemplo 2: " << normal3 << "\n";
    clock_t fimn3 = clock() - inicio;
    cout << "Tempo:" << double(fimn3)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cnormal(l_real(1.0),l_real(0.0),-1.0,1.0,50)) << endl;
    cout << "Erro Absoluto :";
    double normal3Real = cnormalReal(1.0,0.0,-1.0,1.0);
    erroAbsoluto(normal3Real,normal3);
    cout << "Erro Relativo :";
    erroRelativo(normal3Real,normal3);
    cout << endl;

    inicio = clock();
    l_interval normal4 = cnormalCaprani(l_real(1.0),l_real(0.0),-1.0,1.0,50);
    cout << "Normal Caprani Exemplo 2: " << normal4 << "\n";
    clock_t fimn4 = clock() - inicio;
    cout << "Tempo:" << double(fimn4)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cnormalCaprani(l_real(1.0),l_real(0.0),-1.0,1.0,50)) << endl;
    cout << "Erro Absoluto :";
    double normal4Real = cnormalReal(1.0,0.0,-1.0,1.0);
    erroAbsoluto(normal4Real,normal4);
    cout << "Erro Relativo :";
    erroRelativo(normal4Real,normal4);
    cout << endl;

    inicio = clock();
    l_interval normal5 = cnormal(l_real(1.0),l_real(0.0),0.5555,0.5555,50);
    cout << "Normal Exemplo 3: " << normal5 << "\n";
    clock_t fimn5 = clock() - inicio;
    cout << "Tempo:" << double(fimn5)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cnormal(l_real(1.0),l_real(0.0),0.5555,0.5555,50)) << endl;
    cout << "Erro Absoluto :";
    double normal5Real = cnormalReal(1.0,0.0,0.5555,0.5555);
    erroAbsoluto(normal5Real,normal5);
    cout << "Erro Relativo :";
    //erroRelativo(normal5Real,normal5);
    cout << endl;


    inicio = clock();
    l_interval normal6 = cnormalCaprani(l_real(1.0),l_real(0.0),0.5555,0.5555,50);
    cout << "Normal Caprani Exemplo 3: " << normal6 << "\n";
    clock_t fimn6 = clock() - inicio;
    cout << "Tempo:" << double(fimn6)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cnormalCaprani(l_real(1.0),l_real(0.0),0.5555,0.5555,50)) << endl;
    cout << "Erro Absoluto :";
    double normal6Real = cnormalReal(1.0,0.0,0.5555,0.5555);
    erroAbsoluto(normal6Real,normal6);
    cout << "Erro Relativo :";
    //erroRelativo(normal6Real,normal6);

    cout << endl;

    inicio = clock();
    l_interval normal7 = cnormal(l_real(1.0),l_real(0.0),-0.5,0.25,50);
    cout << "Normal Exemplo 4: " << normal7 << "\n";
    clock_t fimn7 = clock() - inicio;
    cout << "Tempo:" << double(fimn7)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cnormal(l_real(1.0),l_real(0.0),-0.5,0.25,50)) << endl;
    cout << "Erro Absoluto :";
    double normal7Real = cnormalReal(1.0,0.0,-0.5,0.25);
    erroAbsoluto(normal7Real,normal7);
    cout << "Erro Relativo :";
    erroRelativo(normal7Real,normal7);
    cout << endl;

    inicio = clock();
    l_interval normal8 = cnormalCaprani(l_real(1.0),l_real(0.0),-0.5,0.25,50);
    cout << "Normal Caprani Exemplo 4: " << normal8 << "\n";
    clock_t fimn8 = clock() - inicio;
    cout << "Tempo:" << double(fimn8)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cnormalCaprani(l_real(1.0),l_real(0.0),-0.5,0.25,50)) << endl;
    cout << "Erro Absoluto :";
    double normal8Real = cnormalReal(1.0,0.0,-0.5,0.25);
    erroAbsoluto(normal8Real,normal8);
    cout << "Erro Relativo :";
    erroRelativo(normal8Real,normal8);
    cout << endl;

    inicio = clock();

    l_interval gama = cgama(0.25,5.0,20.0,40.0);
    cout << "Gama Exemplo 1 : " << gama << "\n";
    clock_t fimgama1 = clock() - inicio;
    cout << "Tempo:" << double(fimgama1)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cgama(0.25,5.0,20.0,40.0)) << endl;
    cout << "Erro Absoluto :";
    double gamaReal = cgamaReal(0.25,5.0,20.0,40.0);
    erroAbsoluto(gamaReal,gama);
    cout << "Erro Relativo :";
    erroRelativo(gamaReal,gama);
    cout << endl;

    inicio = clock();
    l_interval gama2 = cgamaCaprani(0.25,5.0,20.0,40.0,50);
    cout << "Gama Caprani Exemplo 1: " << gama2 << "\n";
    clock_t fimgama2 = clock() - inicio;
    cout << "Tempo:" << double(fimgama2)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cgamaCaprani(0.25,5.0,20.0,40.0,50)) << endl;
    cout << "Erro Absoluto :";

    erroAbsoluto(gamaReal,gama2);
    cout << "Erro Relativo :";
    erroRelativo(gamaReal,gama2);

    cout << endl;

    inicio = clock();
    l_interval gama3 = cgama(0.81,7.81,0.0,10.0);
    cout << "Gama Exemplo 2: " << gama3 << "\n";
    clock_t fimgama3 = clock() - inicio;
    cout << "Tempo:" << double(fimgama3)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cgama(0.81,7.81,0.0,10.0)) << endl;
    cout << "Erro Absoluto :";
    double gama2Real = cgamaReal(0.81,7.81,0.0,10.0);
    erroAbsoluto(gama2Real,gama3);
    cout << "Erro Relativo :";
    erroRelativo(gama2Real,gama3);
    cout << endl;

    inicio = clock();
    l_interval gama4 = cgamaCaprani(0.81,7.81,0.00001,10.0,50);
    cout << "Gama Caprani Exemplo 2: " << cgamaCaprani(0.81,7.81,0.00001,10.0,50) << "\n";
    clock_t fimgama4 = clock() - inicio;
    cout << "Tempo:" << double(fimgama4)/CLOCKS_PER_SEC<<endl;
    cout << "Diâmetro:" << diam(cgamaCaprani(0.81,7.81,0.00001,10.0,50)) << endl;
    cout << "Erro Absoluto :";
    erroAbsoluto(gama2Real,gama4);
    cout << "Erro Relativo :";
    erroRelativo(gama2Real,gama4);

    cout << endl;

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

l_interval cnormalCaprani(l_real alpha, l_real mi, double x1, double x2, int n){

    //Primeiro iniciar vetor divisores com os limites a partir do N
    double h = (x2 - x1)/n;
    double infAux = x1;

    vector<l_interval> divs;

    //Inicializa vetores com as divisoes por N
    for(int i =0;i<n;i++){
        divs.push_back(l_interval(infAux,infAux+h));
        infAux+=h;

    }

    l_interval aux;
    l_interval answer = l_interval(0.0,0.0);
    for(int i=0;i<n;i++){
        aux = l_interval(0.0,0.0);

        //Min
        l_real numeradorMin = pow(Inf(divs[i]) - mi,l_real(2.0));
        l_real denominador = l_real(2.0)*pow(alpha,l_real(2.0));
        l_real auxFormula = numeradorMin/denominador;
        l_real firstNormal = (l_real(1.0)/(alpha*sqrt(2.0*M_PI)));
        l_real secondNormalMin = pow(l_real(M_E),-auxFormula);

        l_real min = l_real(firstNormal*secondNormalMin);


        //Max
        l_real numeradorMax = pow(Sup(divs[i]) - mi,l_real(2.0));
        auxFormula = numeradorMax/denominador;
        l_real secondNormalMax = pow(l_real(M_E),-auxFormula);

        l_real max = l_real(firstNormal*secondNormalMax);

        //Calcula Metade
        l_real middleAux = (Sup(divs[i])+Inf(divs[i]))/l_real(2.0);
        l_real numeradorMid = pow(middleAux - mi,l_real(2.0));
        auxFormula = numeradorMid/denominador;
        l_real secondNormalMid = pow(l_real(M_E),-auxFormula);
        l_real mid = l_real(4.0)*l_real(firstNormal*secondNormalMid);

        l_interval formulaFinal = l_interval(min)+l_interval(max)+l_interval(mid);
        formulaFinal = formulaFinal*((Sup(divs[i])-Inf(divs[i]))/6.0);

        answer+=formulaFinal;
    }


    //while
    //Colocar metodo de Caprani, usando como inferior e superior os limites do vetor divisores
    //Recalcular diametro/2 e seguir formula de SIMPSON
    return answer;
}


l_interval cgamaCaprani(double lambda, double v, double x1, double x2,int n){
    //Primeiro iniciar vetor divisores com os limites a partir do N
    double h = (x2 - x1)/n;
    double infAux = x1;

    vector<l_interval> divs;

    //Inicializa vetores com as divisoes por N
    for(int i =0;i<n;i++){
        divs.push_back(l_interval(infAux,infAux+h));
        infAux+=h;

    }

    l_real firstStep = l_real(firstStepGama(lambda,v));
    l_interval aux;
    l_interval answer = l_interval(0.0,0.0);


    for(int i=0;i<n;i++) {
        aux = l_interval(0.0, 0.0);

        l_real firstMin = cxsc::pow(l_real(M_E), -l_real(lambda) * Inf(divs[i]));
        l_real firstMax = cxsc::pow(l_real(M_E), -l_real(lambda) * Sup(divs[i]));

        l_real middleAux = (Sup(divs[i]) + Inf(divs[i])) / l_real(2.0);
        l_real firstMid = cxsc::pow(l_real(M_E), -l_real(lambda) * middleAux);


        l_real min = firstMin * (cxsc::pow(Inf(divs[i]), l_real(v) - l_real(1.0)));
        l_real max = firstMax * (cxsc::pow(Sup(divs[i]), l_real(v) - l_real(1.0)));
        l_real mid = l_real(4.0) * (firstMid * (cxsc::pow(middleAux, v - l_real(1.0))));

        l_interval formulaFinal = l_interval(min) + l_interval(max) + l_interval(mid);
        formulaFinal = formulaFinal * ((Sup(divs[i]) - Inf(divs[i])) / 6.0);
        answer += formulaFinal;
    }

    return answer/firstStep;
}

l_interval cexpoCaprani(l_real alpha, l_interval x,int n){
    //Primeiro iniciar vetor divisores com os limites a partir do N
    l_real h = (Sup(x) - Inf(x))/l_real(n);
    l_real infAux = Inf(x);

    vector<l_interval> divs;

    //Inicializa vetores com as divisoes por N
    for(int i =0;i<n;i++){
        divs.push_back(l_interval(infAux,infAux+h));
        infAux+=h;
    }

    l_interval aux;
    l_interval answer = l_interval(0.0,0.0);
    for(int i=0;i<n;i++){
        aux = l_interval(0.0,0.0);

        l_real min = alpha*(cxsc::pow(l_real(M_E), -Inf(divs[i]) * alpha));
        l_real max = alpha*(cxsc::pow(l_real(M_E), -Sup(divs[i]) * alpha));

        l_real middleAux = (Sup(divs[i])+Inf(divs[i]))/l_real(2.0);

        l_real mid = l_real(4.0)*alpha*(cxsc::pow(l_real(M_E), (-middleAux * alpha)));

        l_interval formulaFinal = l_interval(min) + l_interval(max) + l_interval(mid);
        formulaFinal = formulaFinal * ((Sup(divs[i]) - Inf(divs[i])) / 6.0);
        answer += formulaFinal;
    }
    return answer;
}

l_interval cparetoCaprani(double alpha, l_real cons, double x1, double x2, int n) {
    //Primeiro iniciar vetor divisores com os limites a partir do N
    double h = (x2 - x1) / n;
    double infAux = x1;

    vector <l_interval> divs;

    //Inicializa vetores com as divisoes por N
    for (int i = 0; i < n; i++) {
        divs.push_back(l_interval(infAux, infAux + h));
        infAux += h;

    }

    l_interval answer = l_interval(0.0, 0.0);

    for (int i = 0; i < n; i++) {

        l_real middleAux = (Sup(divs[i])+Inf(divs[i]))/l_real(2.0);
        l_real denomiMin = (cxsc::pow(Inf(divs[i]),l_real(alpha)+l_real(1.0)));
        l_real denomiMax = (cxsc::pow(Sup(divs[i]),l_real(alpha)+l_real(1.0)));
        l_real denomiMid = (cxsc::pow(middleAux,l_real(alpha)+l_real(1.0)));


        l_real nume = l_real(alpha)*cxsc::pow(cons,l_real(alpha));


        l_real min = nume/denomiMin;

        l_real max = nume/denomiMax;

        l_real mid = l_real(4.0)*(nume/denomiMid);

        l_interval formulaFinal = l_interval(min) + l_interval(max) + l_interval(mid);
        formulaFinal = formulaFinal * ((Sup(divs[i]) - Inf(divs[i])) / 6.0);
        answer += formulaFinal;
    }
    return answer;

}



double firstStepGama(double lambda, double v){
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

    return aux3;
}
//(exp(-L*a)*a**(v-1))

double cgamaReal(double lambda, double v, double x1, double x2) {


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

    return gamaMin/aux3;

}

double cnormalReal(double alpha, double mi, double x1, double x2){


    double numeradorMin = std::pow(x1 - mi,2.0);

    double numeradorMax = std::pow(x2 - mi,2.0);

    double denominador = 2.0*std::pow(alpha,2.0);


    double aux = numeradorMin/denominador;
    double aux2 = numeradorMax/denominador;


    double firstNormal = (1.0/(alpha*sqrt(2*M_PI)));
    double secondNormalMin = std::pow(M_E,-aux);
    double secondNormalMax = std::pow(M_E,-aux2);
    int n = 500;
    double h = (x1 - x2)/n;
    double min = (firstNormal*secondNormalMin) + (firstNormal*secondNormalMax);


    double value = x2;

    //Subdivisoes
    int i=0;
    while(i<(n-1)){

        value+=h;
        //Update values
        numeradorMin = std::pow(value - mi,2.0);
        numeradorMax = std::pow(value - mi,2.0);

        aux = numeradorMin/denominador;
        aux2 = numeradorMax/denominador;

        secondNormalMin = std::pow(M_E,-aux);
        secondNormalMax = std::pow(M_E,-aux2);
        //
        if(i%2 == 0){ //Par
            min+= 4*(firstNormal*secondNormalMin);
        }
        else{
            min+= 2*(firstNormal*secondNormalMin);


        }
        i++;
    }

    min = -min*(h/3.0);


    return min;
}

double cparetoReal(double alfa, double cons, double x1, double x2,int n){

    double first = (alfa*(std::pow(cons,alfa)))/(x1*(alfa+1));
    double last = (alfa*(std::pow(cons,alfa)))/(x2*(alfa+1));

    double paretoSimp = first + last;

    double aux = x1;
    double h = (x2-x1)/n;

    for(int i =0;i < n-1 ; i++){
        aux+=h;
        if(i%2==0){
            paretoSimp+=4*((alfa*(std::pow(cons,alfa)))/(std::pow(aux,(alfa+1))));
        }
        else{
            paretoSimp+=2*((alfa*(std::pow(cons,alfa)))/(std::pow(aux,(alfa+1))));
        }

    }

    paretoSimp = (h/3)*paretoSimp;


    return paretoSimp;
}

double cexpoReal(double alpha, double x1,double x2, int n){

    double aux1 = alpha*(std::pow(M_E,-alpha*x1));
    double aux2 = alpha*(std::pow(M_E,-alpha*x2));

    double resp = aux1 + aux2;


    double h = (x2-x1)/n;
    double value = x1;

    for(int i =0;i<n-1;i++){
        value+=h;

        if(i%2==0){

            resp+= 4*(alpha*(std::pow(M_E,-alpha*value)));
        }
        else{
            resp+= 2*(alpha*(std::pow(M_E,-alpha*value)));
        }
    }

    resp = (h/3)*resp;

    return resp;

}

double cunifReal(double x1,double x2){

    return 1/(x2-x1);
}

void erroAbsoluto(double xReal, l_interval x){
    l_real abs1 = abs(xReal - mid(x));
    l_real abs2 = (Sup(x)-Inf(x))/2;
    cout << abs1 << " < " << abs2 << endl;

}
void erroRelativo(double xReal, l_interval x){
    l_real abs1 = abs(xReal - mid(x));
    l_real rel1 = abs1/xReal;
    l_real rel2;

    if(Inf(x)==l_real(0)){
       rel2 = 0.0;
    }
    else{
        rel2 = (Sup(x)-Inf(x)) / ( 2.0 * Inf(x) );}

    cout << rel1 << " <= " << rel2 << endl;

}






