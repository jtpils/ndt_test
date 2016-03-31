#include "PSO.h"

/* メイン */
int main(void){

  Particle P;

  #ifdef TIME        // 実行時間確認
  struct timeval start, end;
  gettimeofday(&start, NULL);  // 計測開始時間
  #endif

  P = ExecPSO(-5,5);

  #ifdef TIME
  gettimeofday(&end, NULL);    // 計測終了時間
  cout << "PSOの実行時間は " << (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)*1.0E-6 << "[s]でした。" << endl;
  #endif

  for(int j=0; j<Nvariables; j++){
    printf("%f ", P->x_star[j]);
  }
  cout << "最適値 = " << P->pbest << endl;

  return 0;
}
