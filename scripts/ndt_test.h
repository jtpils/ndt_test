#ifndef NDT_TEST___
#define NDT_TEST___

#define DIM 2

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

struct points
{
  int index;
  double x_pos;
  double y_pos;
};

struct gridData
{
  int xIdx;
  int yIdx;
  int size;
  points sum;
  double temp_mat[2][2]; //x * xt
};

void loadFile(std::string filename, std::vector<points>& data, std::string file_type);
void estimateBB2d(std::vector<points> data, double& min_x, double& min_y, double& max_x, double& max_y);
//void data_sorting(std::vector<points>& target_data);

//行列演算用関数（掛け算）
void multiVtM(double vec[DIM], double mat[DIM][DIM], double vec_h[DIM]);//1*2
void multiMV(double mat[DIM][DIM], double vec[DIM], double vec_v[DIM]);//2*1
void multiMV3d(double mat[DIM+1][DIM+1], double vec[DIM+1], double vec_v[DIM+1]);//3*1
double multiVtV(double vec_h[DIM], double vec_v[DIM]);

//逆行列を求める関数
bool estimateInvMat(double inv_mat[DIM][DIM], double mat[DIM][DIM]);
bool estimateInvMat3d(double mat[DIM+1][DIM+1], double inv_mat[DIM+1][DIM+1]);

//複数のポイントの重心を求める
void estimateMean(points& mean_pt, std::vector<points>& data);

//複数のポイントの共分散行列を求める
void estimateCovarianceMat(double mat[DIM][DIM], std::vector<points>& data, points mean_pt);

void saveNDT(std::vector<gridData> ndt_description, double min_x, double min_y, double max_x, double max_y);

double estimateProb(points pt, double inv_mat[DIM][DIM], points mean_pt);
double estimateScore(std::vector<points>& data, points mean_pt, double cov_mat[DIM][DIM]);

//座標変換を施す
void transformPoint(points source, double theta, double trans_x, double trans_y, points& out_put);

//ヤコビアンを返す
//idで指定した要素を返す
//thetaはradで頼む
void estimateJacov(double Jacov_vec[2], points pt, double theta, int id);

//ヤコビアンの微分を返す
//row_id, column_idを指定してやったらその要素を返す
//thetaはradで頼む
void estimateGradJacov(double GJacov_vec[2], points pt, double theta, int row_id, int column);

//評価関数の勾配を求める関数
//グリッド内の全点に適応してあとでvec_gを足してね
void estimateVecG(points pt, double transform_param[3], points mean_pt, double cov_inv_mat[2][2], double vec_g[3]);

//評価関数の勾配の勾配を求める関数
//グリッド内の全点に適応してあとでvec_gを足してね
void estimateMatH(points pt, double transform_param[3], points mean_pt, double cov_inv_mat[2][2], double mat_h[3][3]);

//行列が正定か判断する関数
//正定でなかったら rambda * I を足しとくよ
bool determineDefinite(double mat[DIM+1][DIM+1], double rambda);


void estimateTransform(std::vector<points>& target_local, std::vector<points>& input_local, double param[DIM+1], double tolerance, int max_cicle, double score);

#endif
