#include <Rcpp.h>
using namespace Rcpp;

//'@importFrom Rcpp evalCpp
//'@useDynLib cauchyTest


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//
//   +1. Chi-kwadrat (chisq.test)
//   +2. Kołomogorow (ks.test)
//   +3. Kuiper (v.test paczka truncgof) ??
//   +4. Cramer von Mises (cvm.test paczka goftest)
//   +4a. modyfikacja CvM (nie znalazłam)
//   +5. Anderson-Darling (ad.test paczka goftest)
//   +6. Watson (nie znalazłam, ale to pochodna cvm, także chyba nie będzie z nim większych problemów)
//   7. Adaptacyjny Neumana (nie znalazłam)
//   +8. ECF (nie znalazłam)
//


//============================================
// Srednia Z_i
//============================================
double mean(NumericVector tab)
{
  int n=tab.size();
  double srednia=0;
  for (int i=0; i<n; i++)
    srednia +=  tab[i]/n;
  return srednia;
}


//============================================
//   maximumu
//============================================
double max( NumericVector tab)
{
  int n = tab.size();
  double max=tab[0];
  for (int i=1;i<n;i++)
    if (tab[i] > max)
      max = tab[i];
    return max;
}


//                              RODZINA CRAMERA-VON MISESA


//'  Test Cramera-von Misesa
//'
//' @param tab niepusty, bez braków danych wektor z danymi.
//' @export
// [[Rcpp::export]]
double CM(NumericVector tab)
{
  int n = tab.size();
  double W=0;

  for (int i=0; i<n; i++)
    W= W + (tab [i] - (double)(2*i+1)/2/n)*(tab [i] - (double)(2*i+1)/2/n);
  W = W + (double) 1/12/n;
  //   printf("W = %f\n",W);
  return W;
}



//'  Test Andersona-Darlinga
//'
//' @param tab niepusty, bez braków danych wektor z danymi.
//' @export
// [[Rcpp::export]]
double AD(NumericVector tab)
{
  int n = tab.size();
  double A2=0;

  for ( int i=0; i<n; i++)
    A2 += (2*i+1)*log(tab[i]) +(2*n +1 -2*(i+1))*log(1-tab[i]);
  A2=-A2/n-n;
  //  printf("AD = %f\n", A2);
  return A2;
}


//'  Test Watsona
//'
//' @param tab niepusty, bez braków danych wektor z danymi.
//' @export
// [[Rcpp::export]]
double Watson(NumericVector tab)
{
  int n = tab.size();
  double U = CM(tab)-(mean(tab)-0.5)*(mean(tab)-0.5)*n;
  //       printf("U = %f\n",U);
  return U;
}



//' Modyfikacja testu Cramera-von Misesa
//'
//' @param tab niepusty, bez braków danych wektor z danymi.
//' @export
// [[Rcpp::export]]
double ModCM(NumericVector tab)
{
  int n = tab.size();
  double M=0;
  double pom_tab_i,pom_tab_i_minus_1;
  for ( int i=0; i<=n; i++)
  {   //----
    if (i<n)
      pom_tab_i=tab[i];
    else
      pom_tab_i=1;
    //----
    if (i==0)
      pom_tab_i_minus_1=0;
    else
      pom_tab_i_minus_1=tab[i-1];
    //----
    //(i<n?tab[i]:1)
    //(i==0?0:tab[i-1])
    M=M+(0.5*fabs(pom_tab_i - (double)i/n)*(pom_tab_i  - (double)i/n)-
      0.5*fabs(pom_tab_i_minus_1 - (double) i/n)*(pom_tab_i_minus_1 - (double)i/n));
  }
  return (double) sqrt(n)*M;
}



//                                    RODZINA KOLMOGOROVA


//'  Test Kuipera
//'
//' @param tab niepusty, bez braków danych wektor z danymi.
//' @export
// [[Rcpp::export]]
double Kuiper(NumericVector tab)
{
  int n = tab.size();
  double V=0;
  NumericVector D_plus(n);

  NumericVector D_minus(n);

  for ( int i=0; i<n; i++)
  {
    D_plus[i]=  (double)(i+1)/n - tab[i];
    D_minus[i]= tab[i] - (double)i/n;
  }
  V=max(D_plus) + max(D_minus);
  // printf ("V = %f\n",V);

  return V;
}


//' Test ECF
//'
//' @param tab niepusty, bez braków danych wektor z danymi.
//' @export
// [[Rcpp::export]]
double ECF(NumericVector tab)
{
  int n = tab.size();
  double S1=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<i; j++)
      S1 +=8./(4 + (tab[j] - tab[i])*(tab[j] - tab[i]));
  double S2=0;
  for (int i=0; i<n; i++)
    S2+=12./(9 + tab[i]*tab[i]);
  return 1 + S1/n - S2 + 0.5*n;
}




