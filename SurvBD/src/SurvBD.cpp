#include <stdlib.h>
#include <R.h>
#include <stdio.h>
#include <vector>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sc_mapping(NumericVector x1, NumericVector y1,
				NumericVector x2, 
				int n1, int n2)
{
	NumericVector y2(n2);
	int j = 0;
	for(int i = 0; i<n2; i++){
		while(x2[i]>=x1[j] and j<n1){
			j += 1;
		}
		if(j == 0){
			y2[i] = 1;
		}else{
			y2[i] = y1[j-1];
		}
	}
	return y2;
}

// [[Rcpp::export]]
double sc_mapping_one(NumericVector x1, NumericVector y1,
				double x2, 
				int n1)
{
	// x1 should be sorted before use. 
	double y2;
	int j = 0;
	while(x2>=x1[j] and j<n1){
		j += 1;
	}
	if(j == 0){
		y2 = 1; 
	}else{
		y2 = y1[j-1];
	}
	return y2;
}


// [[Rcpp::export]]
double SurvBD_IPW(NumericVector time1, NumericVector time2,
                     NumericVector delta1, NumericVector delta2,
                     NumericVector sc1, NumericVector sc2,
                     NumericVector sc1AtSamp2, NumericVector sc2AtSamp1,
                     int n1, int n2)
{
  int i, j, u, v;
  double aix_x = 0.0, aix_y = 0.0, aiy_x = 0.0, aiy_y = 0.0;
  double cix_x = 0.0, cix_y = 0.0, ciy_x = 0.0, ciy_y = 0.0;
  double aix = 0.0, aiy = 0.0, cix = 0.0, ciy = 0.0;
  double A = 0.0, C = 0.0, D;
  double min_gamma_t, max_gamma_t;
  double sc1_min_gamma_t, sc1_max_gamma_t, sc2_min_gamma_t, sc2_max_gamma_t;

  // A
  for(i=0; i<n1; i++) {
    for(j=0; j<n1; j++){
      if(delta1[i] == 1 && delta1[j] == 1){
        //
        min_gamma_t = min(2*time1[i]-time1[j], time1[j]);
        max_gamma_t = max(2*time1[i]-time1[j], time1[j]);
        sc1_min_gamma_t = sc_mapping_one(time1, sc1, min_gamma_t, n1);
        sc1_max_gamma_t = sc_mapping_one(time1, sc1, max_gamma_t, n1);
        sc2_min_gamma_t = sc_mapping_one(time2, sc2, min_gamma_t, n2);
        sc2_max_gamma_t = sc_mapping_one(time2, sc2, max_gamma_t, n2);
        // Left
        // 
        for(u=0; u<n1; u++){
          if(time1[u] > min_gamma_t){
            aix_x += 1;
          }
          if(time1[u] > max_gamma_t){
            aix_y += 1;
          }
        }
        //
        if(aix_x == 0){
          aix_x = 0;
        } else{
          aix_x = aix_x / sc1_min_gamma_t;
        }
        if(aix_y == 0){
          aix_y = 0;  
        } else{
          aix_y = aix_y / sc1_max_gamma_t;
        }
        //
        aix = 1.0 / n1 * (aix_x - aix_y);
		
		// Right
        //   
        for(v=0; v<n2; v++){
          if(time2[v] > min_gamma_t){
            aiy_x += 1;
          }
          if(time2[v] > max_gamma_t){
            aiy_y += 1;
          }
        }
        //
        if(aiy_x == 0){
          aiy_x = 0;
        }else{
          aiy_x = aiy_x / sc2_min_gamma_t;
        }
        if(aiy_y == 0){
          aiy_y = 0;
        } else{
          aiy_y = aiy_y / sc2_max_gamma_t;
        }
        //
        aiy = 1.0 / n2 * (aiy_x - aiy_y);   

        //
        A += 1.0 / sc1[i] / sc1[j] * (aix - aiy) * (aix - aiy);  
        
        // reset
        aix_x = 0;
        aix_y = 0;
        aiy_x = 0;
        aiy_y = 0;
      }    
    }
  }
  A = A / n1 / n1;

  // C
  for(i=0; i<n2; i++) {
    for(j=0; j<n2; j++){
      if(delta2[i] == 1 && delta2[j] == 1){
      	//
        min_gamma_t = min(2*time2[i]-time2[j], time2[j]);
        max_gamma_t = max(2*time2[i]-time2[j], time2[j]);
        sc1_min_gamma_t = sc_mapping_one(time1, sc1, min_gamma_t, n1);
        sc1_max_gamma_t = sc_mapping_one(time1, sc1, max_gamma_t, n1);
        sc2_min_gamma_t = sc_mapping_one(time2, sc2, min_gamma_t, n2);
        sc2_max_gamma_t = sc_mapping_one(time2, sc2, max_gamma_t, n2);
      	// Left
        // 
        for(u=0; u<n1; u++){
          if(time1[u] > min_gamma_t){
            cix_x += 1;
          }
          if(time1[u] > max_gamma_t){
            cix_y += 1;
          }
        }
        //
        if(cix_x == 0){
          cix_x = 0;
        } else{
          cix_x = cix_x / sc1_min_gamma_t;
        }
        if(cix_y == 0){
          cix_y = 0;  
        } else{
          cix_y = cix_y / sc1_max_gamma_t;
        }
        //
        cix = 1.0 / n1 * (cix_x - cix_y);
        
        //
        for(v=0; v<n2; v++){
          if(time2[v] > min_gamma_t){
            ciy_x += 1;
          }
          if(time2[v] > max_gamma_t){
            ciy_y += 1;
          }
        }
        //
        if(ciy_x == 0){
          ciy_x = 0;
        }else{
          ciy_x = ciy_x / sc2_min_gamma_t;
        }
        if(ciy_y == 0){
          ciy_y = 0;
        } else{
          ciy_y = ciy_y / sc2_max_gamma_t;
        }
        //
        ciy = 1.0 / n2 * (ciy_x - ciy_y);   

        //
        C += 1.0 / sc2[i] / sc2[j] * (cix - ciy) * (cix - ciy);  
        
        // reset
        cix_x = 0;
        cix_y = 0;
        ciy_x = 0;
        ciy_y = 0;
      }    
    }
  }
  C = C / n2 / n2;

  D = A + C;
  
  return D;
}



// [[Rcpp::export]]
double SurvBD_IPW_accelerated(NumericVector time1, NumericVector time2,
                     NumericVector delta1, NumericVector delta2,
                     NumericVector sc1, NumericVector sc2,
                     NumericVector rank_X, NumericVector rank_Y, 
                     NumericVector rank_XY,
                     NumericVector rank_XX, NumericVector rank_YY, 
                     NumericVector rank_1_1, NumericVector rank_1_2,
                     NumericVector rank_2_1, NumericVector rank_2_2,
                     int n1, int n2)
{
  int i, j, ind;
  double aix_x = 0.0, aix_y = 0.0, aiy_x = 0.0, aiy_y = 0.0;
  double cix_x = 0.0, cix_y = 0.0, ciy_x = 0.0, ciy_y = 0.0;
  double aix = 0.0, aiy = 0.0, cix = 0.0, ciy = 0.0;
  double A = 0.0, C = 0.0, D;
  double sc1_min_gamma_t, sc1_max_gamma_t, sc2_min_gamma_t, sc2_max_gamma_t;

  // A
  for(i=0; i<n1; i++) {
    for(j=0; j<n1; j++){
      if(delta1[i] == 1 && delta1[j] == 1){
        //
        ind = j * n1 + i;
  
        // Left
        if(i<=j){	
          aix_x = n1 - (rank_1_1[n1+ind] - rank_XX[ind]); // sum(I(X_u > 2X_i-X_j))
          aix_y = n1 - rank_X[j]; // sum(I(X_u > X_j))
          if(rank_1_1[n1+ind] - rank_XX[ind] == 0){
            sc1_min_gamma_t = 1;
          }else{
            sc1_min_gamma_t = sc1[rank_1_1[n1+ind] - rank_XX[ind] - 1];
          }
          if(rank_X[j] == 0){
            sc1_max_gamma_t = 1;
          } else{
            sc1_max_gamma_t = sc1[rank_X[j] - 1];
          }
        } else{
        	aix_x = n1 - rank_X[j];
          aix_y = n1 - (rank_1_1[n1+ind] - rank_XX[ind]);
          if(rank_X[j] == 0){
            sc1_min_gamma_t = 1;
          } else{
            sc1_min_gamma_t = sc1[rank_X[j] - 1];
          }
          if(rank_1_1[n1+ind] - rank_XX[ind] == 0){
            sc1_max_gamma_t = 1;
          }else{
            sc1_max_gamma_t = sc1[rank_1_1[n1+ind] - rank_XX[ind] - 1];
          }
        }
        if(aix_x == 0){
          aix_x = 0;
        } else{
          aix_x = aix_x / sc1_min_gamma_t;
        }
        if(aix_y == 0){
          aix_y = 0;  
        } else{
          aix_y = aix_y / sc1_max_gamma_t;
        }
        //
        aix = 1.0 / n1 * (aix_x - aix_y);

		
		    // Right
        if(i<=j){
          aiy_x = n2 - (rank_2_1[n2+ind] - rank_XX[ind]); // sum(I(Y_v > 2X_i-X_j))
          aiy_y = n2 - (rank_XY[j] - rank_X[j]); // sum(I(Y_v > X_j))
          if(rank_2_1[n2+ind] - rank_XX[ind] == 0){
            sc2_min_gamma_t = 1;
          } else{
            sc2_min_gamma_t = sc2[rank_2_1[n2+ind] - rank_XX[ind] - 1];
          }
          if(rank_XY[j] - rank_X[j] == 0){
            sc2_max_gamma_t = 1;
          }else{
            sc2_max_gamma_t = sc2[rank_XY[j] - rank_X[j] - 1];
          }
        }else{
          aiy_x = n2 - (rank_XY[j] - rank_X[j]);
          aiy_y = n2 - (rank_2_1[n2+ind] - rank_XX[ind]);
          if(rank_XY[j] - rank_X[j] == 0){
            sc2_min_gamma_t = 1;
          }else{
            sc2_min_gamma_t = sc2[rank_XY[j] - rank_X[j] - 1];
          }
          if(rank_2_1[n2+ind] - rank_XX[ind] == 0){
            sc2_max_gamma_t = 1;
          } else{
            sc2_max_gamma_t = sc2[rank_2_1[n2+ind] - rank_XX[ind] - 1];
          }
        }
        if(aiy_x == 0){
          aiy_x = 0;
        }else{
          aiy_x = aiy_x / sc2_min_gamma_t;
        }
        if(aiy_y == 0){
          aiy_y = 0;
        } else{
          aiy_y = aiy_y / sc2_max_gamma_t;
        }
        //
        aiy = 1.0 / n2 * (aiy_x - aiy_y);   

        //
        A += 1.0 / sc1[i] / sc1[j] * (aix - aiy) * (aix - aiy);  
        
        // reset
        aix_x = 0;
        aix_y = 0;
        aiy_x = 0;
        aiy_y = 0;
      }    
    }
  }
  A = A / n1 / n1;

  // C
  for(i=0; i<n2; i++) {
    for(j=0; j<n2; j++){
      if(delta2[i] == 1 && delta2[j] == 1){

        ind = j * n2 + i;
        // Left
        if(i<=j){	
          cix_x = n1 - (rank_1_2[n1+ind] - rank_YY[ind]); // sum(I(X_u > 2Y_i-Y_j))
          cix_y = n1 - (rank_XY[n1+j] - rank_Y[j]);; // sum(I(X_u > Y_j))
          if(rank_1_2[n1+ind] - rank_YY[ind] == 0){
            sc1_min_gamma_t = 1;
          }else{
            sc1_min_gamma_t = sc1[rank_1_2[n1+ind] - rank_YY[ind] - 1];
          }
          if(rank_XY[n1+j] - rank_Y[j] == 0){
            sc1_max_gamma_t = 1;
          } else{
            sc1_max_gamma_t = sc1[rank_XY[n1+j] - rank_Y[j] - 1];
          }
        } else{
          cix_y = n1 - (rank_1_2[n1+ind] - rank_YY[ind]); // sum(I(X_u > 2Y_i-Y_j))
          cix_x = n1 - (rank_XY[n1+j] - rank_Y[j]);; // sum(I(X_u > Y_j))
          if(rank_1_2[n1+ind] - rank_YY[ind] == 0){
            sc1_max_gamma_t = 1;
          }else{
            sc1_max_gamma_t = sc1[rank_1_2[n1+ind] - rank_YY[ind] - 1];
          }
          if(rank_XY[n1+j] - rank_Y[j] == 0){
            sc1_min_gamma_t = 1;
          } else{
            sc1_min_gamma_t = sc1[rank_XY[n1+j] - rank_Y[j] - 1];
          }
        }
        if(cix_x == 0){
          cix_x = 0;
        } else{
          cix_x = cix_x / sc1_min_gamma_t;
        }
        if(cix_y == 0){
          cix_y = 0;  
        } else{
          cix_y = cix_y / sc1_max_gamma_t;
        }
        //
        cix = 1.0 / n1 * (cix_x - cix_y);

		
		// Right
        if(i<=j){
          ciy_x = n2 - (rank_2_2[n2+ind] - rank_YY[ind]); // sum(I(Y_v > 2Y_i-Y_j))
          ciy_y = n2 - rank_Y[j]; // sum(I(Y_v > Y_j))
          if(rank_2_2[n2+ind] - rank_YY[ind] == 0){
            sc2_min_gamma_t = 1;
          } else{
            sc2_min_gamma_t = sc2[rank_2_2[n2+ind] - rank_YY[ind] - 1];
          }
          if(rank_Y[j] == 0){
            sc2_max_gamma_t = 1;
          }else{
            sc2_max_gamma_t = sc2[rank_Y[j] - 1];
          }
        }else{
          ciy_y = n2 - (rank_2_2[n2+ind] - rank_YY[ind]); // sum(I(Y_v > 2Y_i-Y_j))
          ciy_x = n2 - rank_Y[j]; // sum(I(Y_v > Y_j))
          if(rank_2_2[n2+ind] - rank_YY[ind] == 0){
            sc2_max_gamma_t = 1;
          } else{
            sc2_max_gamma_t = sc2[rank_2_2[n2+ind] - rank_YY[ind] - 1];
          }
          if(rank_Y[j] == 0){
            sc2_min_gamma_t = 1;
          }else{
            sc2_min_gamma_t = sc2[rank_Y[j] - 1];
          }
        }
        if(ciy_x == 0){
          ciy_x = 0;
        }else{
          ciy_x = ciy_x / sc2_min_gamma_t;
        }
        if(ciy_y == 0){
          ciy_y = 0;
        } else{
          ciy_y = ciy_y / sc2_max_gamma_t;
        }
        ciy = 1.0 / n2 * (ciy_x - ciy_y);

        //
        C += 1.0 / sc2[i] / sc2[j] * (cix - ciy) * (cix - ciy);  
        
        // reset
        cix_x = 0;
        cix_y = 0;
        ciy_x = 0;
        ciy_y = 0;
      }    
    }
  }
  C = C / n2 / n2;

  D = A + C;
  
  return D;
}


// [[Rcpp::export]]
double SurvBD_KMW(NumericVector time1, NumericVector time2,
                     NumericVector delta1, NumericVector delta2,
                     NumericVector weight1, NumericVector weight2,
                     int n1, int n2)
{
  int i, j, u, v;
  double aix = 0.0, aiy = 0.0, cix = 0.0, ciy = 0.0;
  double A = 0.0, C = 0.0, D;
  double min_gamma_t, max_gamma_t;
  //
  for(i=0; i<n1; i++){
  	for(j=0; j<n1; j++){
		min_gamma_t = min(2*time1[i]-time1[j], time1[j]);
		max_gamma_t = max(2*time1[i]-time1[j], time1[j]);
		for(u=0; u<n1; u++){
		  if(min_gamma_t <= time1[u] && time1[u] <= max_gamma_t){
		  	aix += weight1[u];
		  }
		}
		for(v=0; v<n2; v++){
		  if(min_gamma_t <= time2[v] && time2[v] <= max_gamma_t){
		  	aiy += weight2[v];
		  }
		}
  		A += weight1[i] * weight1[j] * (aix - aiy) * (aix - aiy);
  		aix = 0;
  		aiy = 0;
  	}
  }
  //
  for(i=0; i<n2; i++){
  	for(j=0; j<n2; j++){
		min_gamma_t = min(2*time2[i]-time2[j], time2[j]);
		max_gamma_t = max(2*time2[i]-time2[j], time2[j]);
		for(u=0; u<n1; u++){
		  if(min_gamma_t <= time1[u] && time1[u] <= max_gamma_t){
		  	cix += weight1[u];
		  }
		}
		for(v=0; v<n2; v++){
		  if(min_gamma_t <= time2[v] && time2[v] <= max_gamma_t){
		  	ciy += weight2[v];
		  }
		}
  		C += weight2[i] * weight2[j] * (cix - ciy) * (cix - ciy);
  		cix = 0;
  		ciy = 0;
  	}
  }
  D = A + C;
  return D;
}
