#include <cmath>
#include <Rcpp.h>
#include<cstdlib>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double soft_thresh(double xx, double lambda) {
  double temp1;
  double temp2;
  temp1 = fabs(xx)-lambda;
  if(temp1<0){
    temp1 = 0.0;
  }
  if(xx>0){
    temp2 = 1.0;
  } else {
    temp2 = -1.0;
  }
  return(temp2*temp1);
}

double Bii(NumericMatrix S,NumericMatrix B,int i, int p){
   double sumterm, out;
   sumterm = 0;
   for(int j = 1; j<=p; j++){
     if(j!=i)
     sumterm += S(i-1,j-1)*B(i-1,j-1);
   }
  out = (-2*sumterm + sqrt(4*sumterm*sumterm + 8*S(i-1,i-1)))/(2*S(i-1,i-1));
  return(out);
}

double Bki(NumericMatrix S, NumericMatrix B, int k, int i, int p, double l){
  double sumterm, out;
  sumterm = 0;
  for(int j = 1; j<=p; j++){
    if(j!=i)
      sumterm += S(i-1,j-1)*B(k-1,j-1);
      }
  out = soft_thresh(-sumterm/S(i-1,i-1), l/(2*S(i-1,i-1)));
  return(out);
}

double qBki(NumericMatrix S,double Bki, NumericMatrix B,int k, int i,int p, double l){
  double sumterm, out1, out2;
  sumterm = 0;
  for(int j = 1; j<=p; j++){
    if(j!=i)
      sumterm += S(i-1,j-1)*B(k-1,j-1);
      }
  out1 = S(i-1,i-1)*Bki*Bki + 2*Bki*sumterm + l*fabs(Bki);
  out2 = 0;
  if(out1<out2){
    return(out1);
    } else {return(out2);}
}

double Bik(NumericMatrix S, NumericMatrix B, int i, int k, int p, double l){
  double sumterm, out;
  sumterm = 0;
  for(int j = 1; j<=p; j++){
    if(j!=k)
      sumterm += S(k-1,j-1)*B(i-1,j-1);
      }
  out = soft_thresh(-sumterm/S(k-1,k-1), l/(2*S(k-1,k-1)));
  return(out);
}

double qBik(NumericMatrix S,double Bki, NumericMatrix B, int i, int k, int p, double l){
  double sumterm, out1, out2;
  sumterm = 0;
  for(int j = 1; j<=p; j++){
    if(j!=k)
      sumterm += S(k-1,j-1)*B(i-1,j-1);
      }
  out1 = S(k-1,k-1)*Bki*Bki + 2*Bki*sumterm + l*fabs(Bki);
  out2 = 0;
  if(out1<out2){
    return(out1);
    } else {return(out2);}
}

void dfs_visit2(int **G,int u, int *color,int &found_cycle, int size){
    if(found_cycle==1){
      return;
    }
    color[u] = 1;
    for(int v = 0;v<size;v++){
        if(G[v][u]==1){
            if(color[v]==1){
              found_cycle = 1;
              return;
            }
            if(color[v]==0){
              dfs_visit2(G,v,color,found_cycle,size);
            }
        }
    }
    color[u] = 2;
}


int dfs2(int **G, int size){
  int color[size];
  for (int u = 0; u<size; u++){
    color[u] = 0;
  }
  int found_cycle = 0;
    for (int u = 0; u<size; u++){
        if(color[u]==0){
          dfs_visit2(G,u,color,found_cycle,size);
        }
        if(found_cycle){
          break;
        }
    }
    return found_cycle;
}


// [[Rcpp::export()]]
NumericMatrix ccdr_cust(NumericMatrix S, NumericMatrix B, double lambda, int maxitr, double tol, bool verbose = false){
  int p = S.ncol();
  NumericMatrix oldB;
  int **graphki = (int**) calloc (p, sizeof (int*));
  for (int i =0; i < p; i++) {
    graphki[i] = (int*) calloc (p, sizeof (int));
  }
  int **graphik = (int**) calloc (p, sizeof (int*));
  for (int i =0; i < p; i++) {
    graphik[i] = (int*) calloc (p, sizeof (int));
  }
  double Bupdateki, Bupdateik, qki, qik, min3;
  int hascycleki, hascycleik;
  oldB = clone(B);

  int r = 1;
  int converged = 0;
  while((r<maxitr)&&(converged==0)){
    if(verbose){
      Rcpp::Rcout << "ccdr at itr" << r << std::endl;
    }
    for(int i=1;i<=p;i++){
      B(i-1,i-1) = Bii(S,B,i,p);
    }
    for(int i=2;i<=p;i++){
      for(int k=1;k<=(i-1);k++){
        Bupdateki = Bki(S,B,k,i,p,lambda);
        Bupdateik = Bik(S,B,i,k,p,lambda);
        hascycleki = 0;
        hascycleik = 0;
        if(verbose){
            Rcpp::Rcout << "setting off-diagonal at i="<< i << "k = " << k <<std::endl;
            Rcpp::Rcout << Bupdateki <<std::endl;
            Rcpp::Rcout << Bupdateik <<std::endl;
        }

        for(int ii=1;ii<=p;ii++){
          for(int jj=1;jj<=p;jj++){
            if(ii!=jj){
              if(fabs(B(ii-1,jj-1))>0){
                graphki[ii-1][jj-1] = 1;
                graphik[ii-1][jj-1] = 1;
              }
            }
            else {
              graphki[ii-1][jj-1] = 0;
              graphik[ii-1][jj-1] = 0;
            }
          }
        }

        if(fabs(Bupdateki)>0){
          graphki[k-1][i-1] = 1;
          if(dfs2(graphki, p)==1){
            hascycleki = 1;
              }
        }

        if(fabs(Bupdateik)>0){
          graphik[i-1][k-1] = 1;
          if(dfs2(graphik, p)==1){
            hascycleki = 1;
              }
        }
        if(hascycleik==1){
          Bupdateik = 0.0;
                } else if(hascycleki==1){
          Bupdateki = 0.0;
                } else {
          qik = qBik(S,Bupdateik,B,i,k,p,lambda);
          qki = qBki(S,Bupdateki,B,k,i,p,lambda);
          min3 = fmin(qik,qki);
          if(min3==0){
            if(verbose){
            Rcpp::Rcout << "min3 is 0"<< min3 << std::endl;
          }
            Bupdateik = Bupdateki = 0;
          }
          else if(min3 == qik){
            if(verbose){
              Rcpp::Rcout << "min3 is qik"<< min3 << std::endl;
            }
            Bupdateik = 0;
          }
          else{
            if(verbose){
              Rcpp::Rcout << "min3 is qki"<< min3 << std::endl;
            }
            Bupdateki = 0;
          }          
          
        }
        B(i-1,k-1) = Bupdateki;
        B(k-1,i-1) = Bupdateik;
      }
    }
    if (Rcpp::max(abs(B - oldB)) < tol)
      {
        converged = 1;
      } else {
      r += 1;
    }
    oldB = clone(B);
  }
  free (graphki);
  free (graphik);
  return(B);
}


// [[Rcpp::export()]]
NumericMatrix partial_cust1(NumericMatrix S, NumericMatrix B, double lambda, int m, int maxitr, double tol,bool verbose = false){
  int p = S.ncol();
  NumericMatrix oldB;
  int **graphki = (int**) calloc (p, sizeof (int*));
  for (int i =0; i < p; i++) {
    graphki[i] = (int*) calloc (p, sizeof (int));
  }
  int **graphik = (int**) calloc (p, sizeof (int*));
  for (int i =0; i < p; i++) {
    graphik[i] = (int*) calloc (p, sizeof (int));
  }
  double Bupdateki, Bupdateik, qki, qik, min3;
  int hascycleki, hascycleik;
  oldB = clone(B);

  int r = 1;
  int converged = 0;
  while((r<maxitr)&&(converged==0)){
    if(verbose){
      Rcpp::Rcout << "partial1 at itr" << r << std::endl;
    }
    for(int i=1;i<=p;i++){
      B(i-1,i-1) = Bii(S,B,i,p);
    }
    for(int i=2;i<=m;i++){
      for(int k=1;k<=(i-1);k++){
        Bupdateki = Bki(S,B,k,i,p,lambda);
        Bupdateik = Bik(S,B,i,k,p,lambda);
        hascycleki = 0;
        hascycleik = 0;
        if(verbose){
              Rcpp::Rcout << "setting off-diagonal at i="<< i << "k = " << k <<std::endl;
              Rcpp::Rcout << Bupdateki <<std::endl;
              Rcpp::Rcout << Bupdateik <<std::endl;
          }

        for(int ii=1;ii<=p;ii++){
          for(int jj=1;jj<=p;jj++){
            if(ii!=jj){
              if(fabs(B(ii-1,jj-1))>0){
                graphki[ii-1][jj-1] = 1;
                graphik[ii-1][jj-1] = 1;
              }
            }
            else {
              graphki[ii-1][jj-1] = 0;
              graphik[ii-1][jj-1] = 0;
            }
          }
        }

        if(fabs(Bupdateki)>0){
          graphki[k-1][i-1] = 1;
          if(dfs2(graphki, p)==1){
            hascycleki = 1;
              }
        }

        if(fabs(Bupdateik)>0){
          graphik[i-1][k-1] = 1;
          if(dfs2(graphik, p)==1){
            hascycleki = 1;
              }
        }
        if(hascycleik==1){
          Bupdateik = 0.0;
                } else if(hascycleki==1){
          Bupdateki = 0.0;
                } else {
          qik = qBik(S,Bupdateik,B,i,k,p,lambda);
          qki = qBki(S,Bupdateki,B,k,i,p,lambda);
          min3 = fmin(qik,qki);
          if(min3==0){
            if(verbose){
            Rcpp::Rcout << "min3 is 0"<< min3 << std::endl;
          }
            Bupdateik = Bupdateki = 0;
          }
          else if(min3 == qik){
            if(verbose){
              Rcpp::Rcout << "min3 is qik"<< min3 << std::endl;
            }
            Bupdateik = 0;
          }
          else{
            if(verbose){
              Rcpp::Rcout << "min3 is qki"<< min3 << std::endl;
            }
            Bupdateki = 0;
          }      
        }
        B(i-1,k-1) = Bupdateki;
        B(k-1,i-1) = Bupdateik;
      }
    }


    for(int i=(m+2);i<=p;i++){
      for(int k=(m+1);k<=(i-1);k++){
        Bupdateki = Bki(S,B,k,i,p,lambda);
        Bupdateik = Bik(S,B,i,k,p,lambda);
        hascycleki = 0;
        hascycleik = 0;
        if(verbose){
              Rcpp::Rcout << "setting off-diagonal at i="<< i << "k = " << k <<std::endl;
              Rcpp::Rcout << Bupdateki <<std::endl;
              Rcpp::Rcout << Bupdateik <<std::endl;
          }

        for(int ii=1;ii<=p;ii++){
          for(int jj=1;jj<=p;jj++){
            if(ii!=jj){
              if(fabs(B(ii-1,jj-1))>0){
                graphki[ii-1][jj-1] = 1;
                graphik[ii-1][jj-1] = 1;
              }
            }
            else {
              graphki[ii-1][jj-1] = 0;
              graphik[ii-1][jj-1] = 0;
            }
          }
        }

        if(fabs(Bupdateki)>0){
          graphki[k-1][i-1] = 1;
          if(dfs2(graphki, p)==1){
            hascycleki = 1;
              }
        }

        if(fabs(Bupdateik)>0){
          graphik[i-1][k-1] = 1;
          if(dfs2(graphik, p)==1){
            hascycleki = 1;
              }
        }
        if(hascycleik==1){
          Bupdateik = 0.0;
                } else if(hascycleki==1){
          Bupdateki = 0.0;
                } else {
          qik = qBik(S,Bupdateik,B,i,k,p,lambda);
          qki = qBki(S,Bupdateki,B,k,i,p,lambda);
          min3 = fmin(qik,qki);
          if(min3==0){
            if(verbose){
            Rcpp::Rcout << "min3 is 0"<< min3 << std::endl;
          }
            Bupdateik = Bupdateki = 0;
          }
          else if(min3 == qik){
            if(verbose){
              Rcpp::Rcout << "min3 is qik"<< min3 << std::endl;
            }
            Bupdateik = 0;
          }
          else{
            if(verbose){
              Rcpp::Rcout << "min3 is qki"<< min3 << std::endl;
            }
            Bupdateki = 0;
          }      
          
        }
        B(i-1,k-1) = Bupdateki;
        B(k-1,i-1) = Bupdateik;
      }
    }

    for(int i=(m+1);i<=p;i++){
      for(int k=1;k<=m;k++){
        B(k-1,i-1) = Bik(S,B,i,k,p,lambda);
        B(i-1,k-1) = 0;
        if(verbose){
              Rcpp::Rcout << "setting off-diagonal at i="<< i << "k = " << k <<std::endl;
              Rcpp::Rcout << B(k-1,i-1) <<std::endl;
              Rcpp::Rcout << B(i-1,k-1) <<std::endl;
          }

      }
    }

    if (Rcpp::max(abs(B - oldB)) < tol)
      {
        converged = 1;
      } else {
      r += 1;
    }
    oldB = clone(B);
  }
  free (graphki);
  free (graphik);
  return(B);
}



// [[Rcpp::export()]]
NumericMatrix partial_cust2(NumericMatrix S, NumericMatrix B, double lambda, int m, int maxitr, double tol, bool verbose = false){
  int p = S.ncol();
  NumericMatrix oldB;
  int **graphki = (int**) calloc (p, sizeof (int*));
  for (int i =0; i < p; i++) {
    graphki[i] = (int*) calloc (p, sizeof (int));
  }
  int **graphik = (int**) calloc (p, sizeof (int*));
  for (int i =0; i < p; i++) {
    graphik[i] = (int*) calloc (p, sizeof (int));
  }
  double Bupdateki, Bupdateik, qki, qik, min3;
  int hascycleki, hascycleik;
  oldB = clone(B);

  int r = 1;
  int converged = 0;
  while((r<maxitr)&&(converged==0)){
    if(verbose){
      Rcpp::Rcout << "partial2 at itr" << r << std::endl;
    }
    for(int i=1;i<=p;i++){
      if(verbose){
        Rcpp::Rcout << "setting diagonal at i="<< i<<std::endl;
      }
      B(i-1,i-1) = Bii(S,B,i,p);
    }
    for(int i=2;i<=m;i++){
      for(int k=1;k<=(i-1);k++){
        Bupdateki = Bki(S,B,k,i,p,lambda);
        Bupdateik = Bik(S,B,i,k,p,lambda);
        hascycleki = 0;
        hascycleik = 0;
        if(verbose){
              Rcpp::Rcout << "setting off-diagonal at i="<< i << "k = " << k <<std::endl;
              Rcpp::Rcout << Bupdateki <<std::endl;
              Rcpp::Rcout << Bupdateik <<std::endl;
          }
          
        for(int ii=1;ii<=p;ii++){
          for(int jj=1;jj<=p;jj++){
            if(ii!=jj){
              if(fabs(B(ii-1,jj-1))>0){
                graphki[ii-1][jj-1] = 1;
                graphik[ii-1][jj-1] = 1;
              }
            }
            else {
              graphki[ii-1][jj-1] = 0;
              graphik[ii-1][jj-1] = 0;
            }
          }
        }

        if(fabs(Bupdateki)>0){
          graphki[k-1][i-1] = 1;
          if(dfs2(graphki, p)==1){
            hascycleki = 1;
              }
        }

        if(fabs(Bupdateik)>0){
          graphik[i-1][k-1] = 1;
          if(dfs2(graphik, p)==1){
            hascycleki = 1;
              }
        }
        if(hascycleik==1){
          Bupdateik = 0.0;
                } else if(hascycleki==1){
          Bupdateki = 0.0;
                } else {
          qik = qBik(S,Bupdateik,B,i,k,p,lambda);
          qki = qBki(S,Bupdateki,B,k,i,p,lambda);
          min3 = fmin(0.0,fmin(qik,qki));
          if(min3==0){
            if(verbose){
              Rcpp::Rcout << "min3 is 0"<< min3 << std::endl;
            }
            Bupdateik = Bupdateki = 0;
          }
          else if(min3 == qik){
            if(verbose){
              Rcpp::Rcout << "min3 is qik"<< min3 << std::endl;
            }
            Bupdateik = 0;
          }
          else{
            if(verbose){
              Rcpp::Rcout << "min3 is pki"<< min3 << std::endl;
            }
            Bupdateki = 0;
          }        
          
          }
        B(i-1,k-1) = Bupdateki;
        B(k-1,i-1) = Bupdateik;
      }
    }

    for(int i=(m+1);i<=p;i++){
      for(int k=1;k<=m;k++){
        B(k-1,i-1) = Bik(S,B,i,k,p,lambda);
        B(i-1,k-1) = 0;
        if(verbose){
              Rcpp::Rcout << "setting off-diagonal at i="<< i << "k = " << k <<std::endl;
              Rcpp::Rcout << B(k-1,i-1) <<std::endl;
              Rcpp::Rcout << B(i-1,k-1) <<std::endl;
          }
      }
    }


    for(int i=(m+2);i<=p;i++){
      for(int k=(m+1);k<=(i-1);k++){
        Bupdateki = Bki(S,B,k,i,p,lambda);
        Bupdateik = Bik(S,B,i,k,p,lambda);
        hascycleki = 0;
        hascycleik = 0;
        if(verbose){
              Rcpp::Rcout << "setting off-diagonal at i="<< i << "k = " << k <<std::endl;
              Rcpp::Rcout << Bupdateki <<std::endl;
              Rcpp::Rcout << Bupdateik <<std::endl;
          }

        for(int ii=1;ii<=p;ii++){
          for(int jj=1;jj<=p;jj++){
            if(ii!=jj){
              if(fabs(B(ii-1,jj-1))>0){
                graphki[ii-1][jj-1] = 1;
                graphik[ii-1][jj-1] = 1;
              }
            }
            else {
              graphki[ii-1][jj-1] = 0;
              graphik[ii-1][jj-1] = 0;
            }
          }
        }

        if(fabs(Bupdateki)>0){
          graphki[k-1][i-1] = 1;
          if(dfs2(graphki, p)==1){
            hascycleki = 1;
              }
        }

        if(fabs(Bupdateik)>0){
          graphik[i-1][k-1] = 1;
          if(dfs2(graphik, p)==1){
            hascycleki = 1;
              }
        }
        if(hascycleik==1){
          Bupdateik = 0.0;
                } else if(hascycleki==1){
          Bupdateki = 0.0;
                } else {
          qik = qBik(S,Bupdateik,B,i,k,p,lambda);
          qki = qBki(S,Bupdateki,B,k,i,p,lambda);
          min3 = fmin(0.0,fmin(qik,qki));
          if(min3==0){
            if(verbose){
              Rcpp::Rcout << "min3 is 0"<< min3 << std::endl;
            }
            Bupdateik = Bupdateki = 0;
          }
          else if(min3 == qik){
            if(verbose){
              Rcpp::Rcout << "min3 is qik"<< min3 << std::endl;
            }
            Bupdateik = 0;
          }
          else{
            if(verbose){
              Rcpp::Rcout << "min3 is pki"<< min3 << std::endl;
            }
            Bupdateki = 0;
          }
        }
        B(i-1,k-1) = Bupdateki;
        B(k-1,i-1) = Bupdateik;
      }
    }

    if (Rcpp::max(abs(B - oldB)) < tol)
      {
        converged = 1;
      } else {
      r += 1;
    }
    oldB = clone(B);
  }
  free (graphki);
  free (graphik);
  return(B);
}
