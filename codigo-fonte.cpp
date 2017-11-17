#include <cmath>
#include <valarray>
#include <iostream>

typedef double type; //tipo a ser utilizado ao longo do calculo: relavante para precisao

using namespace std;

int iterBusca=0; //numero de iteracoes de bsuca (secao aurea ou Armijo)
int iterMetodo=0; //numero de iteracoes do metodo de minimizacao (gradiente, Newton ou quase-Newton)

type f(valarray<type> x){ //calcula valor da funcao objetivo f(x)=x^2+(exp(x)-y)^2 no ponto x
       return pow(x[0],2)+pow(exp(x[0])-x[1],2);
}

valarray<type> df(valarray<type> x){ //computa vetor gradiente da funcao objetivo f(x)=x^2+(exp(x)-y)^2 no ponto x
        type arr[2]={2*x[0]+2*(exp(x[0])-x[1])*exp(x[0]),2*(exp(x[0])-x[1])*(-1)};
        valarray<type> ret(arr,2);
        return ret;
}

valarray<type> ddf(valarray<type> x){ //retorna hessiana da funcao objetivo no ponto x
        type arr[4]={2+4*exp(2*x[0])-2*exp(x[0])*x[1],-2*exp(x[0]),-2*exp(x[0]),2}; //xx, xy, yx, yy
        valarray<type> ret(arr,4);
        return ret;
}

valarray<type> inverte(valarray<type> M){ //retorna inversa de M
        type arr[4]={M[3],-M[1],-M[2],M[0]};
        valarray<type> ret(arr,4);
        return ret/(M[0]*M[3]-M[1]*M[2]);
}

type secaoAurea(valarray<type> x, valarray<type> d){ //busca por secao aurea
     //COMPLETAR
}

type armijo(valarray<type> x, valarray<type> d, type gama=0.59, type n=0.4){ //busca de Armijo
       type t=1.0;
       while (f(x+t*d)>f(x)+n*t*(df(x)*d).sum()){ //outras condicoes de paradas nao sao necessarias pois e globalmente convergente
             t*=gama;
             iterBusca++;
       }
       return t;
}

valarray<type> gradiente(valarray<type> x, bool usarArmijo=true, type tol=0.00001){ //metodo do gradiente
     valarray<type> nabda=df(x);
     valarray<type> x0(0.0,2);
     while (pow(nabda[0],2)+pow(nabda[1],2)>pow(tol,2)&&(x0!=x)[0]&&(x0!=x)[1]&&iterMetodo<5000){ //condicao de parada: gradiente aprox igual a 0, duas iteracoes com mesmo otimo e numero de iteracoes maximo (5000)
           x0=x;
           if (usarArmijo)
              x+=armijo(x,-nabda)*-nabda; //chamada a busca de Armijo
           else
               x+=secaoAurea(x,-nabda)*-nabda;//chamada a busca por secao aurea
           nabda=df(x);
           iterMetodo++;
     }
     return x;
}

valarray<type> newton(valarray<type> x, bool usarArmijo=true, type tol=0.00001){ //metodo de Newton
     valarray<type> nabda=df(x);
     valarray<type> invHessiana=inverte(ddf(x));
     valarray<type> d(0.0,2);
     d[0]=-(invHessiana[0]*nabda[0]+invHessiana[1]*nabda[1]);
     d[1]=-(invHessiana[2]*nabda[0]+invHessiana[3]*nabda[1]);
     valarray<type> x0(0.0,2);
     while (pow(nabda[0],2)+pow(nabda[1],2)>pow(tol,2)&&(x0!=x)[0]&&(x0!=x)[1]&&iterMetodo<5000){ //condicao de parada: gradiente aprox igual a 0, duas iteracoes com mesmo otimo e numero de iteracoes maximo (5000)
           x0=x;
           if (usarArmijo)
              x+=armijo(x,-nabda)*d; //chamada a busca de Armijo
           else
               x+=secaoAurea(x,-nabda)*d;//chamada a busca por secao aurea
           nabda=df(x);
           invHessiana=inverte(ddf(x));
           d[0]=-(invHessiana[0]*nabda[0]+invHessiana[1]*nabda[1]);
           d[1]=-(invHessiana[2]*nabda[0]+invHessiana[3]*nabda[1]);
           iterMetodo++;
     }
     return x;
}

int main(){
    type tempX[2]={1,1}; //ponto inicial
    valarray<type> x(tempX,2);
    valarray<type> xOtimo=newton(x); //chamada do metodo (e definicao da busca)
    cout<<"iterMetodo = "<<iterMetodo<<endl;
    cout<<"iterBusca = "<<iterMetodo<<endl;
    cout<<"xOtimo = "<<xOtimo[0]<<", "<<xOtimo[1]<<endl; //valor otimo encontrado pelo metodo com a busca utilizada
    cout<<"fOtimo = "<<f(xOtimo)<<endl; //valor da funcao objetivo no ponto obtido
    valarray<type> gradXOtimo=df(xOtimo);
    cout<<"erro = "<<sqrt(pow(gradXOtimo[0],2)+pow(gradXOtimo[1],2))<<endl; //valor do gradiente no ponto computado como otimo
    getchar();
    return 0;
}
