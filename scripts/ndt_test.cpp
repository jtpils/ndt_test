#include "ndt_test.h"

using namespace std;
using namespace boost;

void loadFile(string filename, std::vector<points>& ayuge, string file_type)
{
    //読み込み
	fstream file;
	string t_str;
	//file.open("input.txt", ios::in);
	file.open(filename.c_str(), ios::in);
	if (!file.is_open()) {
		cout << "data file is not exists!" << endl;
		exit(1);
	}
	if(file_type == "pcd")
	{
        //１１行読み飛ばす
		for(int i = 0; i < 11; i++)
			getline(file, t_str);
		cerr << "file type:= pcd" << endl;
	}
	if(file_type == "txt")
	{
        //何もしない
		cerr << "file type:= txt" << endl;
	}
	int count = 0;
	while(getline(file, t_str))
	{
		vector<string> ori_data;
		vector<double> data;
		//" "で区切る
		char_separator<char> seq(" ");
		tokenizer< char_separator<char> > tokens(t_str, seq);
		typedef tokenizer< char_separator<char> >::iterator Iter;
		Iter init = tokens.begin();
		//vectorに格納
		for (Iter it = tokens.begin(); it != tokens.end(); ++it){
			ori_data.push_back(*it);
		}
		points pt;
		pt.x_pos = lexical_cast<double>(ori_data[0]);
		pt.y_pos = lexical_cast<double>(ori_data[1]);
		//pt.z = lexical_cast<double>(ori_data[2]);
		//cerr << "[" << count << "] (" << pt.x_pos << "," << pt.y_pos << ") OK!" << endl;
		count++;
		if( (pt.x_pos != 0) && (pt.y_pos != 0) )
			ayuge.push_back(pt);
	}
	file.close();
	cerr << "load OK!" << endl;
}

void estimateBB2d(std::vector<points> data, double& min_x, double& min_y, double& max_x, double& max_y)
{
	min_x = 9999;
	max_x = -9999;
	min_y = 9999;
	max_y = -9999;
	for(int i = 0; i < data.size(); i++)
	{
		if(data.at(i).x_pos < min_x) min_x = data.at(i).x_pos;
		if(data.at(i).x_pos > max_x) max_x = data.at(i).x_pos;
		if(data.at(i).y_pos < min_y) min_y = data.at(i).y_pos;
		if(data.at(i).y_pos > max_y) max_y = data.at(i).y_pos;
	} 
	cerr << "(" << min_x << "," << min_y << ") (" << max_x << "," << max_y << ")" << endl;
}

void data_sorting(std::vector<points>& target_data)
{
    //データを一旦y座標の情報を元にソートするよ
	sort(target_data.begin(), target_data.end(), 
		 [](const points& lhs, const points& rhs)->bool
		 { return lhs.y_pos < rhs.y_pos; });
	
	int min_y_id = target_data.at(0).yIdx;
	int x_start = 0;
	for (int i = 0; i < target_data.size(); i++){
        //id_切り替わり時にx座標の情報をもとにソート
		if(min_y_id < target_data.at(i).yIdx){
			int x_end = i;
			sort(&target_data.at(x_start), &target_data.at(x_end), 
				 [](const points& lhs, const points& rhs)->bool
				 { return lhs.x_pos < rhs.x_pos; });
            //最低idの更新
			min_y_id = target_data.at(i).yIdx;
			x_start = i;
		}
	}
	sort(&target_data.at(x_start), &target_data.at( target_data.size()-1 ), 
		 [](const points& lhs, const points& rhs)->bool
		 { return lhs.x_pos < rhs.x_pos; });

	cerr << "sortしたの？" << endl;
	for (int i = 0; i < target_data.size(); i++)
		cerr << "[" << i << "] (" << target_data.at(i).x_pos << "," << target_data.at(i).y_pos << "," << target_data.at(i).xIdx << "," << target_data.at(i).yIdx << ")" << endl;
}


//vec * mat = 1*2
void multiVtM(double vec[2], double mat[2][2], double vec_h[2])
{
	vec_h[0] = vec[0]*mat[0][0] + vec[1]*mat[1][0];
	vec_h[1] = vec[0]*mat[0][1] + vec[1]*mat[1][1]; 
}

//mat * vec = 2*1
void multiMV(double mat[2][2], double vec[2], double vec_v[2])
{
	vec_v[0] = mat[0][0]*vec[0]+ mat[0][1]*vec[1];
	vec_v[1] = mat[1][0]*vec[0]+ mat[1][1]*vec[1];
}

//mat * vec = 3*1
void multiMV3d(double mat[3][3], double vec[3], double vec_v[3])
{
	for(int i = 0; i < 3; i++)
		vec_v[i] = mat[i][0] * vec[0] + mat[i][1]*vec[1] + mat[i][2] * vec[2];
}


//vec(1*2) * vec(2*1) = a
double multiVtV(double vec_h[2], double vec_v[2])
{
	double value = vec_h[0]*vec_v[0]+ vec_h[1]*vec_v[1];
	return value;
}

void estimateMean(points& mean_pt, std::vector<points>& data)
{
	double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
	if (data.size() == 0) 
		exit(-1);
	for(int i = 0; i < data.size(); i++){
		x_sum += data.at(i).x_pos;
		y_sum += data.at(i).y_pos;
		//z_sum += data.at(i).z;
	}
	mean_pt.x_pos = x_sum / data.size();
	mean_pt.y_pos = y_sum / data.size();
	//mean_pt.z = z_sum / data.size();
}

void estimateCovarianceMat(double mat[2][2] ,std::vector<points>& data, points mean_pt)
{
	double a00 = 0.0, a01 = 0.0;
	double a10 = 0.0, a11 = 0.0;
	for(int i = 0; i < data.size(); i++){
		a00 += pow( (data.at(i).x_pos - mean_pt.x_pos), 2);
		a11 += pow( (data.at(i).y_pos - mean_pt.y_pos), 2);
		a01 += (data.at(i).x_pos - mean_pt.x_pos) * (data.at(i).y_pos - mean_pt.y_pos);
	}
	a00 /= data.size(); a01 /= data.size(); 
	a10 = a01; a11 /= data.size(); 
	mat[0][0] = a00;
	mat[0][1] = a01;
	mat[1][0] = a10;
	mat[1][1] = a11;
}

void estimateInvMat(double inv_mat[2][2], double mat[2][2])
{
	double delta = (mat[0][0] * mat[1][1]) - (mat[0][1] * mat[1][0]);
	inv_mat[0][0] = mat[1][1] / delta;
	inv_mat[1][1] = mat[0][0] / delta;
	inv_mat[0][1] = -mat[0][1] /delta;
	inv_mat[1][0] = -mat[1][0] / delta; 
}

void estimateInvMat3d(double mat[3][3], double inv_mat[3][3])
{
	double det = mat[0][0]*mat[1][1]*mat[2][2] + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1] - mat[0][0]*mat[1][2]*mat[2][1] - mat[0][1]*mat[1][0]*mat[2][2] - mat[0][2]*mat[1][1]*mat[2][0];
	if (det == 0){
		cerr << "no inv mat" << endl;
		exit(-1);
	}
	inv_mat[0][0] = (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]) / det;
	inv_mat[0][1] = (mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2]) / det;
	inv_mat[0][2] = (mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1]) / det;
	inv_mat[1][0] = (mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2]) / det;
	inv_mat[1][1] = (mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0]) / det;
	inv_mat[1][2] = (mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2]) / det;
	inv_mat[2][0] = (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]) / det;
	inv_mat[2][1] = (mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1]) / det;
	inv_mat[2][2] = (mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]) / det;
}

double estimateProb(points pt, double inv_mat[2][2], points mean_pt)
{
	double diff_x = pt.x_pos - mean_pt.x_pos; 
	double diff_y = pt.y_pos - mean_pt.y_pos;
	double buf1 = (diff_x * inv_mat[0][0] + diff_y * inv_mat[1][0]) * diff_x;
	double buf2 = (diff_x * inv_mat[0][1] + diff_y * inv_mat[1][1]) * diff_y;
	double buf = -(buf1 + buf2) / 2; 
	double output_data = exp(buf);
	return output_data;
}

void transformPoint(points source, double theta, double trans_x, double trans_y, points& out_put)
{
	double transform_matrix[2][2];
	//transform_matrix[0][0] = cos(theta * M_PI /180.0);
	//transform_matrix[0][1] = -sin(theta * M_PI /180.0);
	transform_matrix[0][0] = cos(theta);//rad
	transform_matrix[0][1] = -sin(theta);//rad
	transform_matrix[1][0] = -transform_matrix[0][1];
	transform_matrix[1][1] = transform_matrix[0][0];
	out_put.x_pos = transform_matrix[0][0] * source.x_pos + transform_matrix[0][1] * source.y_pos + trans_x;
	out_put.y_pos = transform_matrix[1][0] * source.x_pos + transform_matrix[1][1] * source.y_pos + trans_y;
}

double estimateScore(std::vector<points>& data, points mean_pt, double cov_mat[2][2])
{
	double inv_mat[2][2];
	estimateInvMat(inv_mat, cov_mat);
	double result = 0.0;
	for (int i = 0; i < data.size(); i++){
		double prob = estimateProb(data.at(i), inv_mat, mean_pt);
		result += prob;
	}
	return result;
}

void estimateJacov(points pt, double theta, double Jacov_mat[2][3])
{
	Jacov_mat[0][0] = 1;
	Jacov_mat[0][1] = 0;
	Jacov_mat[0][2] = -pt.x_pos * sin(theta * M_PI / 180.0) - pt.y_pos * cos(theta * M_PI /180.0);
	Jacov_mat[1][0] = 0;
	Jacov_mat[1][1] = 1;
	Jacov_mat[1][2] = pt.x_pos * cos(theta * M_PI / 180.0) - pt.y_pos * sin(theta * M_PI /180.0);
}

void estimateVecG(points pt, double theta, points mean_pt, double cov_mat[2][2], double vec_g[3])
{
//contents of exp
	double vec_qit[2];
	vec_qit[0] = pt.x_pos - mean_pt.x_pos; vec_qit[1] = pt.y_pos - mean_pt.y_pos;
	double buf_vec[2];
	double inv_mat[2][2];
	estimateInvMat(inv_mat, cov_mat);
	multiVtM(vec_qit, inv_mat, buf_vec);

	double contents = -multiVtV(buf_vec, vec_qit) / 2.0;
	double jacov[2][3];
	estimateJacov(pt, theta, jacov);

	vec_g[0] = (buf_vec[0] * jacov[0][0] + buf_vec[1] * jacov[1][0]) * exp(contents);
	vec_g[1] = (buf_vec[0] * jacov[0][1] + buf_vec[1] * jacov[1][1]) * exp(contents);
	vec_g[2] = (buf_vec[0] * jacov[0][2] + buf_vec[1] * jacov[1][2]) * exp(contents);

}

void estimateMatH(points pt, double theta, points mean_pt, double cov_mat[2][2], double mat_h[3][3])
{
	double jacov[2][3];
	estimateJacov(pt, theta, jacov);
	double vec_qi[2];
	vec_qi[0] = pt.x_pos - mean_pt.x_pos; vec_qi[1] = pt.y_pos - mean_pt.y_pos;
	double buf_vec[2];
	double inv_mat[2][2];
	estimateInvMat(inv_mat, cov_mat);
	multiVtM(vec_qi, inv_mat, buf_vec);

	double contents = -multiVtV(buf_vec, vec_qi) / 2.0;
// mat_h[0][0]
// dq/dp0 = {1, 0}
// dq^2/dp0dp0 = {0, 0}
	for(int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			double pi_pn[2] = {jacov[0][i], jacov[1][i]};
			double pi_pm[2] = {jacov[0][j], jacov[1][j]};
			double pipi_pnpm[2];
			if(i == 2 && j == 2){
				pipi_pnpm[0] = -pt.x_pos * cos(theta) + pt.y_pos * sin(theta);
				pipi_pnpm[1] = -pt.x_pos * sin(theta) - pt.y_pos * cos(theta);
			}
			else{
				pipi_pnpm[0] = 0.0;
				pipi_pnpm[1] = 0.0;
			}

			// cerr << "[" << i << "][" << j << "]:" << endl 
			// 	 << "  dpi/dpn:(" << pi_pn[0] << ", " << pi_pn[1] << ")" << endl 
			// 	 << "  dpi/dpm:(" << pi_pm[0] << ", " << pi_pm[1] << ")" << endl 
			// 	 << "  d^2pi/dpndpm:(" << pipi_pnpm[0] << ", " << pipi_pnpm[1] << ")" << endl;

			double buf_vec2[2];
			multiVtM(pi_pm, inv_mat, buf_vec2);
			mat_h[i][j] = ( ( -multiVtV(buf_vec, pi_pn) ) * ( -multiVtV(buf_vec, pi_pm) )
							- multiVtV(buf_vec, pipi_pnpm) //n=m=3以外では0
							- multiVtV(buf_vec2, pi_pn) ) * exp(contents);
		}
	}
}

bool determineDefinite(double mat[3][3])
{
	double det12 = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
	double det23 = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
	double det123 = mat[0][0] * mat[1][1] * mat[2][2] 
		+ mat[0][1] * mat[1][2] * mat[2][0]
		+ mat[0][2] * mat[1][0] * mat[2][1]
		- mat[0][0] * mat[1][2] * mat[2][1]
		- mat[0][1] * mat[1][0] * mat[2][2]
		- mat[0][2] * mat[1][1] * mat[2][0]; 
	if(mat[0][0] >= 0){
		if(det12 >= 0)
			if(det23 >=0)
				if(det123 >= 0){
					cerr << "正定!!" << endl;
					return true;
				}
				else{
					cerr << "正定ではない(final stage false!)" << endl;
					return false;
				}
			else{
				cerr << "正定ではない(3rd stage false)" << endl;
				return false;
			}
		else{
			cerr << "正定ではない(2nd stage false)" << endl;
			return false;
		}
	}
	else{
		cerr << "正定ではない(1st stage false)" << endl;
		return false;
	}
}

void estimateTransform(std::vector<points>& target_local, std::vector<points>& input_local, double param[DIM+1], double tolerance, int max_cicle, double score)
{
	cerr << "local_target.size():" << target_local.size() << endl;
	cerr << "local_input.size():" << input_local.size() << endl;

	//平均を求めるお
	points mean_pt;
	estimateMean(mean_pt, target_local);
	cerr << "mean pt: " << mean_pt.x_pos << ", " << mean_pt.y_pos << std::endl;

	//共分散行列を求めるお
	double mat[DIM][DIM];
	estimateCovarianceMat(mat, target_local, mean_pt);
	cerr << "cov mat:" << endl;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			cerr << " " << mat[i][j];
		}
		cerr << endl;
	}

	//共分散行列を求めるお
	double inv_mat_test[DIM][DIM];
	estimateInvMat(inv_mat_test, mat);
	cerr << "inv cov mat:" << endl;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			cerr << " " << inv_mat_test[i][j];
		}
		cerr << endl;
	}

	double sum_vec_g[2];
	double sum_mat_h[3][3];

	if (input_local.size() <= 3){
		cerr << " entry few points cource..." << endl;
		for (int i = 0; i < DIM+1; i++)
			param[i] = 999;
		score = -999;
	}

	else {
		cerr << " entry normal cource..." << endl;
		for(int i = 0; i < max_cicle; i++){
            //gとhを求める
			for (int index = 0; index < input_local.size(); index++){		
				//estimate g
				double vector_g[2];
				estimateVecG(input_local.at(index), 0.0, mean_pt, mat, vector_g);
				//add g
				sum_vec_g[0] += vector_g[0];
				sum_vec_g[1] += vector_g[1];

				//estimate H
				double mat_h[3][3];
				estimateMatH(input_local.at(index), 0.0, mean_pt, mat, mat_h);
				
				//add H
				for(int i = 0; i < 3; i++){
					for (int j = 0; j < 3; j++)
						sum_mat_h[i][j] += mat_h[i][j];
				}
			}
			
			if(i == 0){
				cerr << "hessian:" << endl;
				for(int i = 0; i < 3; i++){
					for(int j = 0; j < 3; j++){
						cerr << " " << sum_mat_h[i][j];
					}
					cerr << endl;
				}
				determineDefinite(sum_mat_h);
			}

			//cerr << "vec_g:" << sum_vec_g[0] << " " << sum_vec_g[1] << endl;
			// cerr << "hessian:" << endl;
			// for(int i = 0; i < 3; i++){
			// 	for(int j = 0; j < 3; j++){
			// 		cerr << " " << sum_mat_h[i][j];
			// 	}
			// 	cerr << endl;
			// }
			
			if(-tolerance < sum_vec_g[0] && sum_vec_g[0] <tolerance)
				if(-tolerance < sum_vec_g[1] && sum_vec_g[1] <tolerance){
					std::cerr << "tolerance clear g = (" << sum_vec_g[0] << ", " << sum_vec_g[1] << ")" << std::endl;
					break;
				}

			double inv_h[3][3];
			estimateInvMat3d(sum_mat_h, inv_h);
			//bool definite = determineDefinite(inv_h);
			
            //変化量を求めるニキ
			multiMV3d(inv_h, sum_vec_g, param);
			param[0] = -param[0];
			param[1] = -param[1];
			param[2] = -param[2];

            //input_localを変換
			for(int i = 0; i < input_local.size(); i++)
				transformPoint(input_local.at(i), param[2], param[0], param[1], input_local.at(i));
			
			score = estimateScore(input_local, mean_pt, mat);

			if(i == max_cicle - 1)
				std::cerr << "converged!!(" << max_cicle << ")" << std::endl;
		}		
	}
}
