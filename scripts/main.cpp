#include "ndt_test.h"
#include <iostream>
#include <fstream>

using namespace std;

string filename = "../data/frame_3.pcd";
string f_type = "pcd";
double resolution = 1.0;

int main(int argc, char** argv)
{
	fstream ofile;
	std::vector<points> ayuge, transformed_ayuge;
	loadFile(filename, ayuge, f_type);
	// if(argc != 4){
	// 	cerr << "<exe> <init x> <init y> <init theta>(on build dir)" << endl;
	// 	return -1;
	// }
	double min_x, min_y, max_x, max_y;
	estimateBB2d(ayuge, min_x, min_y, max_x, max_y);

    //番号振り分け
	for(int i = 0; i < ayuge.size(); i++){
		int x_number = (int)( (ayuge.at(i).x_pos - min_x) / resolution );
		int y_number = (int)( (ayuge.at(i).y_pos - min_y) / resolution );
		ayuge.at(i).xIdx = x_number;
		ayuge.at(i).yIdx = y_number;
	}

    //データを一旦ylabelの情報を元にソートするよ
	sort(ayuge.begin(), ayuge.end(), 
		 [](const points& lhs, const points& rhs)->bool
		 { return lhs.y_pos < rhs.y_pos; });
	
	cerr << "ylabelのsortしたの？" << endl;
	int min_y_id = ayuge.at(0).yIdx;
	int x_start = 0;
	for (int i = 0; i < ayuge.size(); i++){
		if(min_y_id < ayuge.at(i).yIdx){
			int x_end = i;
			sort(&ayuge.at(x_start), &ayuge.at(x_end), 
				 [](const points& lhs, const points& rhs)->bool
				 { return lhs.x_pos < rhs.x_pos; });
			min_y_id = ayuge.at(i).yIdx;
			x_start = i;
		}
	}

	sort(&ayuge.at(x_start), &ayuge.at( ayuge.size()-1 ), 
		 [](const points& lhs, const points& rhs)->bool
		 { return lhs.x_pos < rhs.x_pos; });

	cerr << "xlabelのsortしたの？" << endl;
	// for (int i = 0; i < ayuge.size(); i++)
	// 	cerr << "[" << i << "] (" << ayuge.at(i).x_pos << "," << ayuge.at(i).y_pos << "," << ayuge.at(i).xIdx << "," << ayuge.at(i).yIdx << ")" << endl;

    //グリッド的なもののサイズ
	double x_size_f = (max_x - min_x) / resolution;
	int x_size = (int)x_size_f;
	double y_size_f = (max_y - min_y) / resolution;
	int y_size = (int)y_size_f;

    //init
	int map[y_size][x_size];
	for(int i = 0; i < y_size; i++)
		for(int j = 0; j < x_size; j++)
			map[i][j] = 0;

	for(int i = 0; i < ayuge.size(); i++){
		map[ayuge.at(i).yIdx][ayuge.at(i).xIdx]++;
	}

	int max_value = -999;
	for(int i = 0; i < y_size; i++)
		for(int j = 0; j < x_size; j++)
			if(max_value < map[i][j]) max_value = map[i][j];

	// cerr << "max value:" << max_value << endl;
	// cv::Mat cv_mat(y_size, x_size, CV_16U);
	// for(int i = 0; i < y_size; i++)
	// 	for(int j = 0; j < x_size; j++)
	// 		cv_mat.at<unsigned short>(i, j) = map[i][j];

	// cv::imshow("test", cv_mat);
	// cv::waitKey();

	// fstream output_file;
	// output_file.open("../data/2d_map.txt", ios::out);
	// for(int i = 0; i < y_size; i++)
	// 	for(int j = 0; j < x_size; j++)
	// 		output_file << i << " " << j << " " << map[i][j] << std::endl;
	// output_file.close();

// 	double init_x = atof(argv[1]);
// 	double init_y = atof(argv[2]);
// 	double init_theta = atof(argv[3]);
// 	cerr << "p (" << init_x << "," << init_y << "," << init_theta << ")" << endl; 
// 	for (int i = 0; i < ayuge.size(); i++){
// 		points trans_pt;
// 		transformPoint(ayuge.at(i), init_theta, init_x, init_y, trans_pt);
// 		transformed_ayuge.push_back(trans_pt);
// 	}

// 	points mean_pt, mean_pt2;
// 	estimateMean(mean_pt, ayuge);
// 	estimateMean(mean_pt2, transformed_ayuge);
  
// 	double mat[2][2], mat2[2][2];
// 	estimateCovarianceMat(mat, ayuge, mean_pt);
// 	estimateCovarianceMat(mat2, transformed_ayuge, mean_pt2);

// 	double inv_mat[2][2], inv_mat2[2][2];
// 	estimateInvMat(inv_mat, mat);
// 	estimateInvMat(inv_mat2, mat2);

// 	ofile.open("../data/output.txt", ios::out);
// 	for(double x = -1.0; x < 2.0; x+= 0.01)
// 		for (double y = -1.0; y < 2.0; y+= 0.01){
// 			points input_pt;
// 			input_pt.x_pos = x;
// 			input_pt.y_pos = y;
// 			double output_data = estimateProb(input_pt, inv_mat, mean_pt);
// 			ofile << x << " " << y << " " << output_data << std::endl;
// 		}
// 	ofile.close();
// 	double score_base = estimateScore(ayuge, mean_pt, mat);

// 	ofile.open("../data/output_transform.txt", ios::out);
// 	for(double x = -1.0; x < 2.0; x+= 0.01)
// 		for (double y = -2.0; y < 2.0; y+= 0.01){
// 			points input_pt;
// 			input_pt.x_pos = x;
// 			input_pt.y_pos = y;
// 			double output_data = estimateProb(input_pt, inv_mat2, mean_pt2);
// 			ofile << x << " " << y << " " << output_data << std::endl;
// 		}
// 	ofile.close();
// 	double score_trans = estimateScore(transformed_ayuge, mean_pt, mat);

// //transformed_ayugeに対するdelta_p を求める
// 	double delta_p[3];
// 	double sum_of_vec_g[2] = {0, 0};
// 	double sum_of_mat_h[3][3] ={ {0, 0, 0},
// 								 {0, 0, 0},
// 								 {0, 0, 0} };
// 	for (int index = 0; index < transformed_ayuge.size(); index++){
// 		//cerr << "[[" << index << "]](" << transformed_ayuge.at(index).x_pos << "," << transformed_ayuge.at(index).y_pos << ")" << endl;
//       //estimate g
// 		double vector_g[2];
// 		estimateVecG(transformed_ayuge.at(index), 0.0, mean_pt, mat, vector_g);
// 		sum_of_vec_g[0] += vector_g[0];
// 		sum_of_vec_g[1] += vector_g[1];
//       //estimate H
// 		double mat_h[3][3];
// 		estimateMatH(transformed_ayuge.at(index), 0.0, mean_pt, mat, mat_h);
// 		for(int i = 0; i < 3; i++)
// 			for (int j = 0; j < 3; j++)
// 				sum_of_mat_h[i][j] += mat_h[i][j];
// 	}
// 	for(int i = 0; i < 3; i++)
// 		for(int j = 0; j < 3; j++)
// 			cerr << "hessian[" << i << "][" << j<< "]=" << sum_of_mat_h[i][j] << endl;
// 	double inv_h[3][3];
// 	estimateInvMat3d(sum_of_mat_h, inv_h);
// 	cerr << inv_h[0][0] << endl;
// 	bool definite = determineDefinite(inv_h);

// 	multiMV3d(inv_h, sum_of_vec_g, delta_p);
// 	for (int i = 0; i < 3; i++)
// 		delta_p[i] = -delta_p[i];

// 	cerr << "apply ndt" << endl;
// 	double new_p[3] = {(delta_p[2]*180.0/M_PI + init_theta), (delta_p[0] + init_x), (delta_p[1] + init_y)};
//     cerr << "delta p(tx):" << delta_p[0] << endl;
// 	cerr << "delta p(ty):" << delta_p[1] << endl;
// 	cerr << "delta p(theta):" << delta_p[2]*180.0/M_PI << endl;
//     cerr << "new p (" << new_p[0] << "," << new_p[1] << "," << new_p[2] << ")" << endl;
	

// //ndt result
// 	std::vector<points> transformed_ndt;
// 	for (int i = 0; i < ayuge.size(); i++){
// 		points trans_pt;
// 		transformPoint(ayuge.at(i), (delta_p[2]*180.0/M_PI + init_theta), (delta_p[0] + init_x), (delta_p[1] + init_y), trans_pt);
// 		transformed_ndt.push_back(trans_pt);
// 	}
// 	points mean_ndt;
// 	estimateMean(mean_ndt, transformed_ndt);
// 	double cov_mat_ndt[2][2];
// 	estimateCovarianceMat(cov_mat_ndt, transformed_ndt, mean_ndt);
// 	double inv_mat_ndt[2][2];
// 	estimateInvMat(inv_mat_ndt, cov_mat_ndt);
// 	ofile.open("../data/ndt_result.txt", ios::out);
// 	for(double x = -1.0; x < 2.0; x+= 0.01)
// 		for (double y = -2.0; y < 2.0; y+= 0.01){
// 			points input_pt;
// 			input_pt.x_pos = x;
// 			input_pt.y_pos = y;
// 			double output_data = estimateProb(input_pt, inv_mat_ndt, mean_ndt);
// 			ofile << x << " " << y << " " << output_data << std::endl;
// 		}
// 	ofile.close();
// 	double result = estimateScore(transformed_ndt, mean_pt, mat);
// 	cerr << "source score:" << score_base << endl;
// 	cerr << "init score:" << score_trans << endl;
// 	cerr << "after ndt score:" << result << endl;
}



