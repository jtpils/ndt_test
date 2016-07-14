#include "ndt_test.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
	using namespace std;
	fstream ofile;
	string target_file = "../data/scan_points1.pcd";
	string input_file = "../data/scan_points2.pcd";
	
	string f_type = "pcd";

    //target data読み込み
	std::vector<points> target_data;
    loadFile(target_file, target_data, f_type);

    //ナンバリングと最大最小座標
	double min_x, min_y, max_x, max_y;
	estimateBB2d(target_data, min_x, min_y, max_x, max_y);
	for(int i = 0; i < target_data.size(); i++){
	  target_data.at(i).index = i;
	}

    //input data読み込み
	std::vector<points> input_data;
	loadFile(input_file, input_data, f_type);

    //ナンバリングと最大最小座標
	double min_x_input, min_y_input, max_x_input, max_y_input;
	estimateBB2d(input_data, min_x_input, min_y_input, max_x_input, max_y_input);
	for(int i = 0; i < input_data.size(); i++){
	  input_data.at(i).index = i;
	 }

	//参照データを保存
	std::cerr << ">>>> 参照データ保存" << std::endl;
	ofile.open("../data/target.txt", ios::out);
	for(int i = 0; i < target_data.size(); i++){
		ofile << target_data.at(i).x_pos << " " << target_data.at(i).y_pos << std::endl;
	}
	ofile.close();

    //入力データを保存
	std::cerr << ">>>> 入力データ保存" << std::endl;
	ofile.open("../data/input.txt", ios::out);
	for(int i = 0; i < input_data.size(); i++){
		ofile << input_data.at(i).x_pos << " " << input_data.at(i).y_pos << std::endl;
	}
	ofile.close();

    //グリッド的なもののサイズ
	double resolution = 2.0;
	double x_size_f = (max_x - min_x) / resolution;
	int x_size = (int)x_size_f;
	double y_size_f = (max_y - min_y) / resolution;
	int y_size = (int)y_size_f;
	std::cerr << "grid_x_size:" << x_size << std::endl;
	std::cerr << "grid_y_size:" << y_size << std::endl;

	//Normal Distributions Transform
	//target dataに対してのみ適用
	//準備
	std::vector<gridData> ndt_description;
	gridData init_data;
	init_data.xIdx = -1;
	init_data.yIdx = -1;
	init_data.size = 0;
	points dummy;
	dummy.x_pos = 0.0; 	dummy.y_pos = 0.0;
	init_data.sum = dummy;
	init_data.temp_mat[0][0] = 0.0; init_data.temp_mat[0][1] = 0.0;
	init_data.temp_mat[1][0] = 0.0; init_data.temp_mat[1][1] = 0.0;
	ndt_description.assign( (x_size + 1) * (y_size + 1), init_data);
	std::cerr << "ndt_description size:" << ndt_description.size() << std::endl;

	//計算
	for(int i = 0; i < target_data.size(); i++){
	  //その点がどのグリッドに相当するのか求める
	  points sample_pt = target_data.at(i);
	  int x_number = (int)( (sample_pt.x_pos - min_x) / resolution ); 
	  int y_number = (int)( (sample_pt.y_pos - min_y) / resolution );
	  int description_id = x_size * y_number + x_number;
	  ndt_description.at( description_id ).xIdx = x_number;
	  ndt_description.at( description_id ).yIdx = y_number;
	  ndt_description.at( description_id ).size += 1;
	  ndt_description.at( description_id ).sum.x_pos += sample_pt.x_pos;
	  ndt_description.at( description_id ).sum.y_pos += sample_pt.y_pos;
	  //temp_matに加算する行列を計算
	  double pt_square[2][2] = 
		{ {sample_pt.x_pos * sample_pt.x_pos, sample_pt.x_pos * sample_pt.y_pos}, 
		  {sample_pt.x_pos * sample_pt.y_pos, sample_pt.y_pos * sample_pt.y_pos}};
	  ndt_description.at( description_id ).temp_mat[0][0] += pt_square[0][0];
	  ndt_description.at( description_id ).temp_mat[0][1] += pt_square[0][1];
	  ndt_description.at( description_id ).temp_mat[1][0] += pt_square[1][0];
	  ndt_description.at( description_id ).temp_mat[1][1] += pt_square[1][1];
	}
	saveNDT(ndt_description, min_x, min_y, max_x, max_y);

    //local_inputをどのように動かせばいいか求めるニキ
	double opt_param[3]; //移動量
	double step; //ステップ数
	double winner_delta_p[3]; //移動量
	double score;
	double winner_score;
	int max_loop = 10;
	double Grad[3];
	double Hesse[3][3];
	int key = 0;

	std::cerr << "opt_param(" << opt_param[0] << ", " << opt_param[1] << ", " << opt_param[2] << ")" <<std::endl;

	std::cerr << "size:" << input_data.size() << std::endl;

	for(int i = 0; i < input_data.size(); i++){
	  ///点群を座標変換
	  points pt_transformed;
	  transformPoint(input_data.at(i), opt_param[2], opt_param[0], opt_param[1], pt_transformed);

	  //サンプリングした点がNDT_descriptionのどのグリッドに相当するか
	  int x_number = (int)( (pt_transformed.x_pos - min_x) / resolution ); 
	  int y_number = (int)( (pt_transformed.y_pos - min_y) / resolution );
	  int ref_id = x_size * y_number + x_number;

	  if(x_number > x_size-1 || y_number > y_size-1) {
		std::cerr << " id over flow!! size:" << ref_id << std::endl;
		continue;
	  }

	  if(ndt_description.at(ref_id).size < 3){
		std::cerr << " target grid data few" << std::endl;
		continue;
	  }

	  std::cerr << "[" << i << "] (" << input_data.at(i).x_pos << ", " << input_data.at(i).y_pos << ") size:" << ref_id << std::endl; 

	  //scoreの勾配を求める準備
	  //参照データの平均
	  points grid_mean;
	  grid_mean.x_pos = ndt_description.at(ref_id).sum.x_pos / ndt_description.at(ref_id).size;
	  grid_mean.y_pos = ndt_description.at(ref_id).sum.y_pos / ndt_description.at(ref_id).size;
	  //参照データの共分散行列の逆行列
	  double mean_pt_v[2] = {grid_mean.x_pos, grid_mean.y_pos};
	  double sum_v[2] = {ndt_description.at(ref_id).sum.x_pos, ndt_description.at(ref_id).sum.y_pos};
	  double buf_mat[2][2];
	  buf_mat[0][0] = mean_pt_v[0] * sum_v[0];
	  buf_mat[0][1] = mean_pt_v[0] * sum_v[1];
	  buf_mat[1][0] = mean_pt_v[1] * sum_v[0];
	  buf_mat[1][1] = mean_pt_v[1] * sum_v[1];
	  double cov[2][2] = 
		{ { (ndt_description.at(ref_id).temp_mat[0][0] - buf_mat[0][0]) / ndt_description.at(ref_id).size ,
			(ndt_description.at(ref_id).temp_mat[0][1] - buf_mat[0][1]) / ndt_description.at(ref_id).size },
		  { (ndt_description.at(ref_id).temp_mat[1][0] - buf_mat[1][0]) / ndt_description.at(ref_id).size ,
			(ndt_description.at(ref_id).temp_mat[1][1] - buf_mat[1][1]) / ndt_description.at(ref_id).size } };
	  double inv_cov[2][2];
	  estimateInvMat(inv_cov, cov);

	  //Gradの要素を求める
	  double grad_[3];
	  estimateVecG(input_data.at(i), opt_param, grid_mean, inv_cov, grad_);
	  //Hesseの要素を求める
	  double hesse_[3][3];
	  estimateMatH(input_data.at(i), opt_param, grid_mean, inv_cov, hesse_);

	  //Grad, Hesseをupdate
	  Grad[0] += grad_[0]; Grad[1] += grad_[1]; Grad[2] += grad_[2];
	  std::cerr << "  grad: " << Grad[0] << ", " << Grad[1] << ", " << Grad[2] << std::endl;
	  Hesse[0][0] += hesse_[0][0]; Hesse[0][1] += hesse_[0][1]; Hesse[0][2] += hesse_[0][2];
	  Hesse[1][0] += hesse_[1][0]; Hesse[1][1] += hesse_[1][1]; Hesse[1][2] += hesse_[1][2];
	  Hesse[2][0] += hesse_[2][0]; Hesse[2][1] += hesse_[2][1]; Hesse[2][2] += hesse_[2][2];
	  std::cerr << "  hesse:" << std::endl;
	  for(int a = 0; a < 3; a++){
		for(int b = 0; b < 3; b++){
		  std::cerr << Hesse[a][b] << " "; 
		}
		std::cerr<< std::endl;
	  }
	}

	//delta_pを求める
	//Hesseの逆行列
	double inv_Hesse[3][3];
	if( estimateInvMat3d(Hesse, inv_Hesse) ){
	  double delta_p[3]; //移動量の変化
	  multiMV3d(inv_Hesse, Grad, delta_p);
	  //paramのupdate
	  opt_param[0] += delta_p[0];
	  opt_param[1] += delta_p[1];
	  opt_param[2] += delta_p[2];
	}
	else{
	 
	}
	std::cerr << "update_param(" << opt_param[0] << ", " << opt_param[1] << ", " << opt_param[2] << ")" <<std::endl;

	return 0;
}

