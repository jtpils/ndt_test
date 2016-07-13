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
	double resolution = 1.0;

    //target data読み込み
	std::vector<points> target_data;
    loadFile(target_file, target_data, f_type);

    //番号振り分け
	double min_x, min_y, max_x, max_y;
	estimateBB2d(target_data, min_x, min_y, max_x, max_y);
	for(int i = 0; i < target_data.size(); i++){
		int x_number = (int)( (target_data.at(i).x_pos - min_x) / resolution );
		int y_number = (int)( (target_data.at(i).y_pos - min_y) / resolution );
		target_data.at(i).xIdx = x_number;
		target_data.at(i).yIdx = y_number;
	}

    //input data読み込み
	std::vector<points> input_data;
	loadFile(input_file, input_data, f_type);

    //番号振り分け
	double min_x_input, min_y_input, max_x_input, max_y_input;
	estimateBB2d(input_data, min_x_input, min_y_input, max_x_input, max_y_input);
	for(int i = 0; i < input_data.size(); i++){
		int x_number = (int)( (input_data.at(i).x_pos - min_x_input) / resolution );
		int y_number = (int)( (input_data.at(i).y_pos - min_y_input) / resolution );
		input_data.at(i).xIdx = x_number;
		input_data.at(i).yIdx = y_number;
	}

	//参照データを保存
	std::cerr << ">>>>>>参照データ" << std::endl;
	ofile.open("../data/target.txt", ios::out);
	for(int i = 0; i < target_data.size(); i++){
		ofile << target_data.at(i).x_pos << " " << target_data.at(i).y_pos << std::endl;
	}
	ofile.close();

    //入力データを保存
	std::cerr << ">>>>>>入力データ" << std::endl;
	ofile.open("../data/input.txt", ios::out);
	for(int i = 0; i < input_data.size(); i++){
		ofile << input_data.at(i).x_pos << " " << input_data.at(i).y_pos << std::endl;
	}
	ofile.close();

    //グリッド的なもののサイズ
	double x_size_f = (max_x - min_x) / resolution;
	int x_size = (int)x_size_f;
	double y_size_f = (max_y - min_y) / resolution;
	int y_size = (int)y_size_f;
	std::cerr << "x_size:" << x_size << std::endl;
	std::cerr << "y_size:" << y_size << std::endl;

    //local_inputをどのように動かせばいいか求めるニキ
	double delta_p[3]; //移動量
	double winner_delta_p[3]; //移動量
	double score;
	double winner_score;
	double tolerance = 0.001;
	int max_loop = 10;
	double sum_vec_g[2];
	double sum_mat_h[3][3];
	int key = 0;

    // for(int i = 0; i < y_size+1; i++){
	// 	for(int j = 0; j < x_size+1; j++){
    for(int i = 0; i < 55; i++){
		for(int j = 10; j < 20; j++){
			//i,jのIdを持つ点群をピックアップ
			std::vector<points> local_target;
			for(int n = 0; n < target_data.size(); n++)
				if( (target_data.at(n).yIdx == i) && (target_data.at(n).xIdx == j) )
					local_target.push_back(target_data.at(n));

			//i,jのIdを持つ入力点群をピックアップ
			std::vector<points> local_input;
			for(int n = 0; n < input_data.size(); n++)
				if( (input_data.at(n).yIdx == i) && (input_data.at(n).xIdx == j) )
					local_input.push_back(input_data.at(n));

            //3点以上あったらオプティマイズスタート
			if(local_input.size() >= 3 && local_target.size() >= 3){
			  cerr << "[" << i << "][" << j << "]" << endl;

			  //許容範囲に達するか、回数回すかでそいやそいや
			  estimateTransform(local_target, local_input, delta_p, tolerance, max_loop, score);
			  
			  cerr << "apply ndt matching[" << i << "][" << j << "]" << endl;
			  cerr << "delta p(tx):" << delta_p[0] << endl;
			  cerr << "delta p(ty):" << delta_p[1] << endl;
			  cerr << "delta p(theta):" << delta_p[2] << endl;
			  cerr << "score:" << score << endl;
			  
			  if(key == 0){
				winner_score = score;
				winner_delta_p[0] = delta_p[0];
				winner_delta_p[1] = delta_p[1];
				winner_delta_p[2] = delta_p[2];
				key = 1;
			  }

			  else if (score > winner_score){
				winner_score = score;
				winner_delta_p[0] = delta_p[0];
				winner_delta_p[1] = delta_p[1];
				winner_delta_p[2] = delta_p[2];
			  }
			}
		}
	}
	cerr << "result" << endl;
	cerr << "delta p(tx):" << winner_delta_p[0] << endl;
	cerr << "delta p(ty):" << winner_delta_p[1] << endl;
	cerr << "delta p(theta):" << winner_delta_p[2] << endl;
	cerr << "score:" << winner_score << endl;
	
    // //入力データtransformしたやつを表示
	// std::cerr << ">>>>>>入力データ（回転）" << std::endl;
	// ofile.open("../data/input_rot.txt", ios::out);
	// for(int i = 0; i < input_data.size(); i++){
	// 	points transformed;
	// 	transformPoint(input_data.at(i), delta_p[2], delta_p[0], delta_p[1], transformed);
	// 	//std::cerr << transformed.x_pos << " " << transformed.y_pos << std::endl;
	// 	ofile << transformed.x_pos << " " << transformed.y_pos << std::endl;
	// }
	ofile.close();
}
