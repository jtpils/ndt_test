#include "PSO.h"

double GetRandom(double min,double max, int digit){
	double ten,R;

	ten = pow(10,digit-1);
	R = min*ten + (int)(rand()*((max-min)*ten+1.0)/(1.0+RAND_MAX));
	return R/ten;
}

/* 評価値計算 */
void Evaluate(Particle P){

	P->f=0.0;
	// P->f = -2*exp(-(pow((P->x[0]-2),2) + pow((P->x[1]-2),2))) - exp(-(pow((P->x[0]+0.5),2) + pow((P->x[1]+0.5),2))/3);
	P->f = -4*exp(-(pow((P->x[0]-2),2) + pow((P->x[1]-2),2))) - 2*exp(-(pow((P->x[0]+0.5),2) + pow((P->x[1]+0.5),2))/3) - 3.5*exp(-(pow((P->x[0]-4),2) + pow((P->x[1]+2),2))/2);
}

/* pbestの更新 */
void UpdateBest(Particle P){
	int j;

	for(j=0; j<Nvariables; j++){
		P->x_star[j]=P->x[j];
	}
	P->pbest=P->f;
}

/* それぞれのお魚さんの情報まとめ */
int Initialize(Particle P, int n, double min, double max){
	int i, j;
	int G;//もっともよいお魚さんの個体値番号

	G=0;
	for(i=0; i<n; i++) {
		for(j=0; j<Nvariables; j++){
			P[i].x[j] = GetRandom(min,max,6);//0から1の乱数発生
			P[i].v[j] = 0.0;//速度は0
		}
		Evaluate(&P[i]);
		UpdateBest(&P[i]);
		if(P[i].f < P[G].f){//Gよりもi番目のお魚さんの評価値がよかったらGをiに代入
			G = i;
		}
	}
	return G;
}

/* 新しいデータ構造の割り当て */
#define New(type, n, msg) (type *)NewCell(sizeof(type), n, msg)

void *NewCell(int size, int n, char *msg){
	void *NEW;
	//newの要素をn個分確保して、それがNULLだったらerrorを返す
	if((NEW=malloc(size*n))==NULL){
		fprintf(stderr, "Cannot allocate memory for %d %s\n", n, msg);
		exit(1);
	}
	return NEW;
}

/* パーティクル型のデータ構造の新しい割り当て */
Particle NewParticles(int n){
	int i;
	Particle P;

	P = New(ParticleRec, n, (char*)"particles");
	for(i=0; i<n; i++){
		P[i].x = New(double, Nvariables, (char*)"x");
		P[i].v = New(double, Nvariables, (char*)"v");
		P[i].x_star = New(double, Nvariables, (char*)"x*");
	}
	return P;
}

/* 発生させたお魚さんの位置と最もよい場所の表示 */
void Print(Particle P){
	int j;

	for(j=0; j<Nvariables; j++){
		printf("%f ", P->x_star[j]);
	}
	printf(" = %g\n", P->pbest);
}


Particle ExecPSO(double min, double max){
	int t, i, j;
	Particle P;
	int G;
	double w;

	P = NewParticles(Nparticles);
	G = Initialize(P, Nparticles, min, max);
	w=W_0;

	for(t=1; t<=T_MAX; t++){

		for(i=0; i<Nparticles; i++){

			for(j=0; j<Nvariables; j++){
				P[i].v[j] = w*P[i].v[j]
				+ c1*(double)rand()/RAND_MAX*(P[i].x_star[j]-P[i].x[j])
				+ c2*(double)rand()/RAND_MAX*(P[G].x_star[j]-P[i].x[j]);
				if(P[i].v[j] < -MAX_V){
					P[i].v[j] = -MAX_V;
				}else if(P[i].v[j] > MAX_V){
					P[i].v[j] = MAX_V;
				}
				P[i].x[j] += P[i].v[j];
			}
			Evaluate(&P[i]);
			if(P[i].f < P[i].pbest){
				if(P[i].f < P[G].pbest) G = i;
				UpdateBest(&P[i]);
			}

		}

		w -= (W_0-W_T)/T_MAX;

	}

	return P;
}


