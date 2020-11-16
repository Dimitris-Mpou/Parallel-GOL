#include <stdio.h>
#include <omp.h>

#define NUM_THREADS 4

void nextStateIn(char **cur, char **next, int n){
  int i, j, k, l, count;

  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for collapse(2)
  for(i=2; i<n; i++){   // Gia olo to grid ektos twn exwterikwn (prasinwn)
    for(j=2; j<n; j++){

      count = 0;
      for(k=i-1; k<=i+1; k++){
        for(l=j-1; l<=j+1; l++){  // Gia kathe geitona
          if(k!=i || l!=j){  // Prospernaw ton eauto tou
            if(cur[k][l] == '1'){
                count++;
            }
          }
        }
      }
      if(count < 2 || count > 3){   // Pethenei apo monaxia i uperplithismo
        next[i][j] = '0';
      }else if(count == 3){         // Epizei i genietai kainourgio
        next[i][j] = '1';
      }else if(cur[i][j] == '1'){  // Epizei
        next[i][j] = '1';
      }else{
        next[i][j] = '0';           // Paramenei nekro
      }
    }
  }
}

void nextStateOut(char **cur, char **next, int n){
  int i, j, k, l, count;

  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel
{
#pragma omp for
  for(i=1; i<=n; i++){   // Gia tin aristeri grammi
    count = 0;
    for(k=i-1; k<=i+1; k++){
      for(l=0; l<=2; l++){  // Gia kathe geitona
        if(k!=i || l!=1){  // Prospernaw ton eauto tou
          if(cur[k][l] == '1'){
              count++;
          }
        }
      }
    }
    if(count < 2 || count > 3){   // Pethenei apo monaxia i uperplithismo
      next[i][1] = '0';
    }else if(count == 3){         // Epizei i genietai kainourgio
      next[i][1] = '1';
    }else if(cur[i][1] == '1'){  // Epizei
      next[i][1] = '1';
    }else{
      next[i][1] = '0';           // Paramenei nekro
    }
  }

#pragma omp for
  for(i=1; i<=n; i++){   // Gia tin dexia stili
    count = 0;
    for(k=i-1; k<=i+1; k++){
      for(l=n-1; l<=n+1; l++){  // Gia kathe geitona
        if(k!=i || l!=n){  // Prospernaw ton eauto tou
          if(cur[k][l] == '1'){
              count++;
          }
        }
      }
    }
    if(count < 2 || count > 3){   // Pethenei apo monaxia i uperplithismo
      next[i][n] = '0';
    }else if(count == 3){         // Epizei i genietai kainourgio
      next[i][n] = '1';
    }else if(cur[i][n] == '1'){  // Epizei
      next[i][n] = '1';
    }else{
      next[i][n] = '0';           // Paramenei nekro
    }
  }

#pragma omp for
  for(j=1; j<=n; j++){    // Panw grammi
    count = 0;
    for(k=0; k<=2; k++){
      for(l=j-1; l<=j+1; l++){  // Gia kathe geitona
        if(k!=1 || l!=j){  // Prospernaw ton eauto tou
          if(cur[k][l] == '1'){
              count++;
          }
        }
      }
    }
    if(count < 2 || count > 3){   // Pethenei apo monaxia i uperplithismo
      next[1][j] = '0';
    }else if(count == 3){         // Epizei i genietai kainourgio
      next[1][j] = '1';
    }else if(cur[1][j] == '1'){  // Epizei
      next[1][j] = '1';
    }else{
      next[1][j] = '0';           // Paramenei nekro
    }
  }
#pragma omp for
  for(j=1; j<=n; j++){    // Katw grammi
    count = 0;
    for(k=n-1; k<=n+1; k++){
      for(l=j-1; l<=j+1; l++){  // Gia kathe geitona
        if(k!=n || l!=j){  // Prospernaw ton eauto tou
          if(cur[k][l] == '1'){
              count++;
          }
        }
      }
    }
    if(count < 2 || count > 3){   // Pethenei apo monaxia i uperplithismo
      next[n][j] = '0';
    }else if(count == 3){         // Epizei i genietai kainourgio
      next[n][j] = '1';
    }else if(cur[n][j] == '1'){  // Epizei
      next[n][j] = '1';
    }else{
      next[n][j] = '0';           // Paramenei nekro
    }
  }
} //  omp pararell end section
}

void findNeighbors(int rank, int sqrt_P, int P, int *up, int *down, int *left, int *right, int *up_left, int *up_right, int *down_left, int *down_right){
  int Pmod = rank%sqrt_P;
	int Pdiv = rank/sqrt_P;  //rank = div*sP+mod


	*down  = ((Pdiv + 1) % sqrt_P)*sqrt_P + Pmod;
	*right  = Pdiv*sqrt_P + (Pmod + 1)%sqrt_P;
	*down_right = ((Pdiv + 1)%sqrt_P)*sqrt_P + (Pmod + 1)%sqrt_P;
	if (Pdiv == 0) {
		*up  = P - sqrt_P + Pmod;
		*up_right = P - sqrt_P + (Pmod + 1)%sqrt_P;
		if (Pmod == 0)
			*up_left = P - 1;
		else
			*up_left = P - sqrt_P + Pmod - 1;
	}
	else {
		*up  = (Pdiv - 1)*sqrt_P + Pmod;
		*up_right = (Pdiv - 1)*sqrt_P + (Pmod + 1)%sqrt_P;
    if (Pmod == 0)
      *up_left = Pdiv*sqrt_P - 1;
    else
      *up_left = (Pdiv - 1)*sqrt_P + Pmod - 1;
	}
	if (Pmod == 0) {
		*left  = (Pdiv + 1)*sqrt_P - 1;
		*down_left =((Pdiv + 1)%sqrt_P)*sqrt_P + sqrt_P - 1;
	}
	else {
		*left  = Pdiv*sqrt_P + Pmod - 1;
		*down_left = ((Pdiv + 1)%sqrt_P)*sqrt_P + Pmod - 1;
	}
}

int sameGrid(char **cur, char **next, int blockSize){
  int i, j;

  omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for collapse(2)
  for(i=1;i<=blockSize;i++){
    for(j=1; j<=blockSize; j++){
      if(cur[i][j] != next[i][j])
        return 0;
    }
  }
  return 1;
}
