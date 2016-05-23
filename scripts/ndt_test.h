#ifndef NDT_TEST___
#define NDT_TEST___

#define DIM 2

#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

struct points
{
	double x_pos;
	double y_pos;
	//double z_pos;
	int xIdx;
	int yIdx;
};

void loadFile(std::string filename, std::vector<points>& data, std::string file_type);
void estimateBB2d(std::vector<points> data, double& min_x, double& min_y, double& max_x, double& max_y);
void data_sorting(std::vector<points>& target_data);
void multiVtM(double vec[DIM], double mat[DIM][DIM], double vec_h[DIM]);//1*2
void multiMV(double mat[DIM][DIM], double vec[DIM], double vec_v[DIM]);//2*1
void multiMV3d(double mat[DIM+1][DIM+1], double vec[DIM+1], double vec_v[DIM+1]);//3*1
double multiVtV(double vec_h[DIM], double vec_v[DIM]);
void estimateMean(points& mean_pt, std::vector<points>& data);
void estimateCovarianceMat(double mat[DIM][DIM], std::vector<points>& data, points mean_pt);
void estimateInvMat(double inv_mat[DIM][DIM], double mat[DIM][DIM]);
void estimateInvMat3d(double mat[DIM+1][DIM+1], double inv_mat[DIM+1][DIM+1]);
double estimateProb(points pt, double inv_mat[DIM][DIM], points mean_pt);
void transformPoint(points source, double theta, double trans_x, double trans_y, points& out_put);
double estimateScore(std::vector<points>& data, points mean_pt, double cov_mat[DIM][DIM]);
void estimateJacov(points pt, double theta, double Jacov_mat[DIM][DIM+1]);
void estimateVecG(points pt, double theta, points mean_pt, double cov_mat[DIM][DIM], double vec_g[DIM]);
void estimateMatH(points pt, double theta, points mean_pt, double cov_mat[DIM][DIM], double mat_h[DIM+1][DIM+1]);
bool determineDefinite(double mat[DIM+1][DIM+1]);
void estimateTransform(std::vector<points>& target_local, std::vector<points>& input_local, double param[DIM+1], double tolerance, int max_cicle, double score);

#endif
