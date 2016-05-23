#include "ndt_test.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
	using namespace std;
	fstream ofile;
	string target_file = "../data/data2.txt";
	string input_file = "../data/data1_rot_1.txt";
	
	string f_type = "txt";
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

	cerr << "after sorting (target data)" << std::endl;
	data_sorting(target_data);
	cerr << "after sorting (input data)" << std::endl;
	data_sorting(input_data);

    //グリッド的なもののサイズ
	double x_size_f = (max_x - min_x) / resolution;
	int x_size = (int)x_size_f;
	double y_size_f = (max_y - min_y) / resolution;
	int y_size = (int)y_size_f;

    //local_inputをどのように動かせばいいか求めるニキ
	double delta_p[3]; //移動量
	double score;
	double tolerance = 0.001;
	int max_loop = 1000;
	double sum_vec_g[2];
	double sum_mat_h[3][3];

    // for(int i = 0; i < y_size+1; i++){
	// 	for(int j = 0; j < x_size+1; j++){
	for(int i = 0; i < 1; i++){
		for(int j = 0; j < 1; j++){
			//i,jのIdを持つ点群をピックアップ
			cerr << "[" << i << "][" << j << "]" << endl;
			std::vector<points> local_target;
			for(int n = 0; n < target_data.size(); n++)
				if( (target_data.at(n).yIdx == i) && (target_data.at(n).xIdx == j) )
					local_target.push_back(target_data.at(n));

			//i,jのIdを持つ入力点群をピックアップ
			std::vector<points> local_input;
			for(int n = 0; n < input_data.size(); n++)
				if( (input_data.at(n).yIdx == i) && (input_data.at(n).xIdx == j) )
					local_input.push_back(input_data.at(n));

            //求めるお？
			cerr << "local_target.size():" << local_target.size() << endl;
			cerr << "local_input.size():" << local_input.size() << endl;

            //許容範囲に達するか、回数回すかでそいやそいや
			estimateTransform(local_target, local_input, delta_p, tolerance, max_loop, score);

			cerr << "apply ndt matching[" << i << "][" << j << "]" << endl;
			cerr << "delta p(tx):" << delta_p[0] << endl;
			cerr << "delta p(ty):" << delta_p[1] << endl;
			cerr << "delta p(theta):" << delta_p[2] << endl;
			cerr << "score:" << score << endl;

		}
	}

    //参照データを表示
	std::cerr << "****結果発表****" << std::endl;
	std::cerr << ">>>>>>参照データ" << std::endl;
	for(int i = 0; i < target_data.size(); i++){
		std::cerr << target_data.at(i).x_pos << " " << target_data.at(i).y_pos << std::endl;
	}

    //入力データを表示
	std::cerr << ">>>>>>入力データ" << std::endl;
	for(int i = 0; i < input_data.size(); i++){
		std::cerr << input_data.at(i).x_pos << " " << input_data.at(i).y_pos << std::endl;
	}

    //入力データtransformしたやつを表示
	std::cerr << ">>>>>>入力データ（回転）" << std::endl;
	for(int i = 0; i < input_data.size(); i++){
		points transformed;
		transformPoint(input_data.at(i), delta_p[2], delta_p[0], delta_p[1], transformed);
		std::cerr << transformed.x_pos << " " << transformed.y_pos << std::endl;
	}
}
