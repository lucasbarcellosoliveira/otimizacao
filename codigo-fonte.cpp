#include <cmath>
#include <valarray>
#include <iostream>

typedef double type; //tipo a ser utilizado ao longo do calculo: relavante para precisao

using namespace std;

int iterBusca=0; //numero de iteracoes de bsuca (secao aurea ou Armijo)
int iterMetodo=0; //numero de iteracoes do metodo de minimizacao (gradiente, Newton ou quase-Newton)
bool convergiu=false; //verdadeiro se metodo foi executado e convergiu dentro do limite de iteracoes

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

valarray<type> multiplica(valarray<type> A, valarray<type> B){ //retorna o produto de A e B 2x2
        if (B.size()==2){ //multiplicacao por vetor
            valarray<type> C(0.0,2);
            C[0]=A[0]*B[0]+A[1]*B[1];
            C[1]=A[2]*B[0]+A[3]*B[1];
            return C;
        }
        else{ //multiplicacao por matriz
            valarray<type> C(0.0,4);
            C[0]=A[0]*B[0]+A[1]*B[2];
            C[1]=A[0]*B[1]+A[1]*B[3];
            C[2]=A[2]*B[0]+A[3]*B[2];
            C[3]=A[2]*B[1]+A[3]*B[3];
            return C;
        }
}

valarray<type> inverte(valarray<type> M){ //retorna inversa de M2x2
        type arr[4]={M[3],-M[1],-M[2],M[0]};
        valarray<type> ret(arr,4);
        return ret/(M[0]*M[3]-M[1]*M[2]);
}

type secaoAurea(valarray<type> x, valarray<type> d, type epsilon=0.001, type ro=1){ //busca por secao aurea
        type teta1=(3-sqrt(5))/2;
        type teta2=1-teta1;
        type a=0, s=ro, b=2*ro;
        while (f(x+b*d)<f(x+s*d)){ //obtencao do invervalo
            a=s; s=b; b*=2;
        }
        type u=a+teta1*(b-a), v=a+teta2*(b-a);
        while (b-a>epsilon){ //obtencao do tamanho do passo
            if (f(x+u*d)<f(x+v*d)){
                b=v; v=u; u=a+teta1*(b-a);
            }
            else{
                a=u; u=v; v=a+teta2*(b-a);
            }
            iterBusca++;
        }
        return (u+v)/2;
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
     valarray<type> nabla=df(x);
     valarray<type> x0(0.0,2);
     while (iterMetodo<5000){ //condicao de parada: numero de iteracoes maximo (5000)
           x0=x;
           if (usarArmijo)
               x+=armijo(x,-nabla)*-nabla; //chamada a busca de Armijo
           else
               x+=secaoAurea(x,-nabla)*-nabla;//chamada a busca por secao aurea
           nabla=df(x);
           iterMetodo++;
           if (pow(nabla[0],2)+pow(nabla[1],2)<=pow(tol,2)||pow(x0[0]-x[0],2)+pow(x0[1]-x[1],2)<=pow(tol,2)){ //condicao de parada: gradiente aprox igual a 0, duas iteracoes com mesmo otimo
               convergiu=true;
               break;
           }
     }
     return x;
}

valarray<type> newton(valarray<type> x, bool puro=false, bool usarArmijo=true, type tol=0.00001){ //metodo de Newton
     valarray<type> nabla=df(x);
     valarray<type> invHessiana=inverte(ddf(x));
     valarray<type> d(0.0,2);
     d=-multiplica(invHessiana,nabla); //multiplicacao de matrizes
     valarray<type> x0(0.0,2);
     while (iterMetodo<5000){ //condicao de parada: iteracoes maximo (5000)
           x0=x;
           if (puro) //metodo de Newton puro com tk=1
               x+=d;
           else
               if (usarArmijo)
                   x+=armijo(x,d)*d; //chamada a busca de Armijo
               else
                   x+=secaoAurea(x,d)*d;//chamada a busca por secao aurea
           nabla=df(x);
           invHessiana=inverte(ddf(x));
           d=-multiplica(invHessiana,nabla); //-hessiana^(-1)grad
           iterMetodo++;
           if (pow(nabla[0],2)+pow(nabla[1],2)<=pow(tol,2)||pow(x0[0]-x[0],2)+pow(x0[1]-x[1],2)<=pow(tol,2)){ //condicao de parada: gradiente aprox igual a 0, duas iteracoes com mesmo otimo
               convergiu=true;
               break;
           }
     }
     return x;
}

valarray<type> quasenewton(valarray<type> x, bool puro=false, bool usarArmijo=true, type tol=0.00001){ //metodo quase-Newton
     valarray<type> nabla=df(x);
     valarray<type> d(0.0,2);
     valarray<type> h(0.0,4); h[0]=1.0; h[2]=1.0; //define h=I
     d=-multiplica(h,nabla); //multiplicacao de matrizes
     valarray<type> p(0.0,4); p[0]=1.0; p[2]=1.0;  //duas diferenca de x entre iteracoes
     valarray<type> q(0.0,4); q[0]=1.0; q[2]=1.0; //duas difereca de gradientes entre iteracoes
     valarray<type> x0(0.0,2);
     while (iterMetodo<5000){ //condicao de parada: iteracoes maximo (5000)
           x0=x;
           if (usarArmijo)
               x+=armijo(x,d)*d; //chamada a busca de Armijo
           else
               x+=secaoAurea(x,d)*d;//chamada a busca por secao aurea
           p[0]=p[2]; p[1]=p[3]; p[2]=x[0]-x0[0]; p[3]=x[1]-x0[1]; //atualiza p
           q[0]=q[2]; q[1]=q[3]; q[2]=df(x)[0]-nabla[0]; q[3]=df(x)[1]-nabla[1];//atualiza q
           nabla=df(x);
           h=multiplica(p,inverte(q)); //pq^(-1)
           d=-multiplica(h,nabla); //-hgrad
           iterMetodo++;
           if (pow(nabla[0],2)+pow(nabla[1],2)<=pow(tol,2)||pow(x0[0]-x[0],2)+pow(x0[1]-x[1],2)<=pow(tol,2)){ //condicao de parada: gradiente aprox igual a 0, duas iteracoes com mesmo otimo
               convergiu=true;
               break;
           }
     }
     return x;
}

int main(){
    type tempX[2]={1,1}; //ponto inicial
    valarray<type> x(tempX,2);
    //valarray<type> xOtimo=gradiente(x, false);
    //valarray<type> xOtimo=newton(x, false, false);
    valarray<type> xOtimo=quasenewton(x, false); //chamada do metodo (e definicao da busca)
    if (convergiu)
       cout<<"Convergencia alcancada!"<<endl;
    else
        cout<<"Convergencia NAO alcancada!"<<endl;
    cout<<"iterMetodo = "<<iterMetodo<<endl;
    cout<<"iterBusca = "<<iterMetodo<<endl;
    cout<<"xOtimo = "<<xOtimo[0]<<", "<<xOtimo[1]<<endl; //valor otimo encontrado pelo metodo com a busca utilizada
    cout<<"fOtimo = "<<f(xOtimo)<<endl; //valor da funcao objetivo no ponto obtido
    valarray<type> gradXOtimo=df(xOtimo);
    cout<<"erro = "<<sqrt(pow(gradXOtimo[0],2)+pow(gradXOtimo[1],2))<<endl; //valor do gradiente no ponto computado como otimo
    getchar();
    return 0;
}
