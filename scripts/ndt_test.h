#ifndef NDT_TEST___
#define NDT_TEST___
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
void multiVtM(double vec[2], double mat[2][2], double vec_h[2]);//1*2
void multiMV(double mat[2][2], double vec[2], double vec_v[2]);//2*1
void multiMV3d(double mat[3][3], double vec[3], double vec_v[3]);//3*1
double multiVtV(double vec_h[2], double vec_v[2]);
void estimateMean(points& mean_pt, std::vector<points>& data);
void estimateCovarianceMat(double mat[2][2], std::vector<points>& data, points mean_pt);
void estimateInvMat(double inv_mat[2][2], double mat[2][2]);
void estimateInvMat3d(double mat[3][3], double inv_mat[3][3]);
double estimateProb(points pt, double inv_mat[2][2], points mean_pt);
void transformPoint(points source, double theta, double trans_x, double trans_y, points& out_put);
double estimateScore(std::vector<points>& data, points mean_pt, double cov_mat[2][2]);
void estimateJacov(points pt, double theta, double Jacov_mat[2][3]);
void estimateVecG(points pt, double theta, points mean_pt, double cov_mat[2][2], double vec_g[2]);
void estimateMatH(points pt, double theta, points mean_pt, double cov_mat[2][2], double mat_h[3][3]);
bool determineDefinite(double mat[3][3]);



#endif
