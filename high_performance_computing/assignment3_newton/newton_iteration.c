#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<math.h>

/*
 * Checks if a point is converged with respect to its distance from roots.
 * These roots all lie on the unit circle, so we identify divergence
 * by being too far away from the unit circle (either too large or too 
 * close to 0). 
 */
int check_convergence(double complex x, const double complex* roots, int d){
// -1 keep doing newton iterations
// ix stop newton iteration convergence to a root and retrun index
// -2 stop newton iteration divergence
  int conv = -1;  
  double rex = creal(x);
  double imx = cimag(x);
  double x2 = rex*rex + imx*imx;
  if(x2 < 0.001*0.001){
    conv = -2;
    return(conv);
  }
  else if(x2 > 1e20){
    conv = -2;
    return(conv);
  }
  if(x2 < (1 - 0.001)*(1 - 0.001) || x2 > (1 + 0.001)*(1 + 0.001)){
    return conv;
  }
  // add if statment and return here for optim
  for(int ix = 0; ix < d; ix++){
    double re = creal(x-roots[ix]);
    double im = cimag(x-roots[ix]);
    if( re*re + im*im < 0.001*0.001){
      conv = ix;
      return conv;
    }
  }
  return conv;
}

/*
 * Performs the Newton step, with hardcoded algorithms depending 
 * on the polynomial degree d.
 */
double complex step_newton(double complex x_old, int d, float c){
  switch(d){
  case 1:
    return 1;
  case 2:
    return c*x_old + 1/(d*x_old);
  case 3:
    return c*x_old + 1/(d*x_old*x_old);
  case 4:
    return c*x_old + 1/(d*x_old*x_old*x_old);
  case 5:
    return c*x_old + 1/(d*x_old*x_old*x_old*x_old);
  case 6:
    return c*x_old + 1/(d*x_old*x_old*x_old*x_old*x_old);
  case 7:
    return c*x_old + 1/(d*x_old*x_old*x_old*x_old*x_old*x_old);
  case 8:
    return c*x_old + 1/(d*x_old*x_old*x_old*x_old*x_old*x_old*x_old);
  case 9:
    return c*x_old + 1/(d*x_old*x_old*x_old*x_old*x_old*x_old*x_old*x_old);
  default:
    printf( "Wrong value for the polynomial degree: d = %d \n",d);
    exit(1);
  }
}
/*
 * Perform Newton iterations for a given starting point x0.
 * Writes the converged root and the number of (possibly capped) iterations
 * to their respective arrays.
 */
void newton_iteration(char* atractor, char* convergence, double complex x0, int d, int cap, const double complex* roots){
  double complex x_old = x0;
  double complex x_new;

  double c = 1 - 1.0/d;
  int conv;
  int ix;
  for(ix = 0; ix<1000; ix++){
    conv = check_convergence(x_old, roots, d);
    if(conv != -1)
      break;
    x_new = step_newton(x_old, d, c);
    x_old = x_new;
  }

  // return the converged value and capped iteration nbr
  if(conv == -1){
    printf( "Newton did not satisfy conv condition nbr of iteration was %d\n", ix);
    printf("The latest value was x_new = %f + i%f \n", creal(x_new), cimag(x_new));
    exit(1);
  } 
  else if (conv == -2){
    *atractor = 0; // divergence always index 0
    *convergence = ix<cap ? ix : cap;
  }
  else{
    *atractor = conv + 1; // +1 since 0 is divergence
    *convergence = ix<cap ? ix : cap;
  }
}
