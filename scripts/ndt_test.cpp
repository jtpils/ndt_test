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

void estimateJacov(double Jacov_vec[2], points pt, double theta, int id)
{
  if(id == 0){
	Jacov_vec[0] = 1;
	Jacov_vec[1] = 0;
  }
  else if(id == 1){
	Jacov_vec[0] = 0;
	Jacov_vec[1] = 1;
  }
  else if(id == 2){
	Jacov_vec[0] = -pt.x_pos * sin(theta) - pt.y_pos * cos(theta);
	Jacov_vec[1] = pt.x_pos * cos(theta) - pt.y_pos * sin(theta);
  }
  else{
	std::cerr << "第四引数は0か1か2で頼む\n" << std::endl;
  }
}

void estimateGradJacov(double GJacov_vec[2], points pt, double theta, int row_id, int column)
{
  if(row_id == 2 && column == 2){
	GJacov_vec[0] = -pt.x_pos * cos(theta) + pt.y_pos * sin(theta);
	GJacov_vec[1] = -pt.x_pos * sin(theta) - pt.y_pos * cos(theta);
  }
  else{
	GJacov_vec[0] = 0.0;
	GJacov_vec[1] = 0.0;
  }
}

void estimateVecG(points pt, double transform_param[3], points mean_pt, double cov_inv_mat[2][2], double vec_g[3])
{
  //transformed point と NDTのグリッドの重心との関係を求める
  points pt_transformed;
  transformPoint(pt, transform_param[2], transform_param[0], transform_param[1], pt_transformed);
  points pt_relative;
  pt_relative.x_pos = pt_transformed.x_pos - mean_pt.x_pos;
  pt_relative.x_pos = pt_transformed.y_pos - mean_pt.y_pos;

  //jacov
  double jacov[3][2];
  estimateJacov(jacov[0], pt_relative, transform_param[2], 0);
  estimateJacov(jacov[1], pt_relative, transform_param[2], 1);
  estimateJacov(jacov[2], pt_relative, transform_param[2], 2);

  //calcurate contents of exp
  double pt_relative_v[2] = {pt_relative.x_pos, pt_relative.y_pos};
  double buf_vec[2];
  multiVtM(pt_relative_v, cov_inv_mat, buf_vec);
  double contents = - multiVtV(buf_vec, pt_relative_v) / 2.0;

  //vec_gを求める
  vec_g[0] = (buf_vec[0] * jacov[0][0] + buf_vec[1] * jacov[0][1]) * exp(contents);
  vec_g[1] = (buf_vec[0] * jacov[1][0] + buf_vec[1] * jacov[1][1]) * exp(contents);
  vec_g[2] = (buf_vec[0] * jacov[2][0] + buf_vec[1] * jacov[2][1]) * exp(contents);
}

void estimateMatH(points pt, double transform_param[3], points mean_pt, double cov_inv_mat[2][2], double mat_h[3][3])
{
  //transformed point と NDTのグリッドの重心との関係を求める
  points pt_transformed;
  transformPoint(pt, transform_param[2], transform_param[0], transform_param[1], pt_transformed);
  points pt_relative;
  pt_relative.x_pos = pt_transformed.x_pos - mean_pt.x_pos;
  pt_relative.x_pos = pt_transformed.y_pos - mean_pt.y_pos;

  //jacov
  double jacov[3][2];
  estimateJacov(jacov[0], pt_relative, transform_param[2], 0);
  estimateJacov(jacov[1], pt_relative, transform_param[2], 1);
  estimateJacov(jacov[2], pt_relative, transform_param[2], 2);

  //calcurate contents of exp
  double pt_relative_v[2] = {pt_relative.x_pos, pt_relative.y_pos};
  double buf_vec[2]; //pt_relative * cov_inv_mat
  multiVtM(pt_relative_v, cov_inv_mat, buf_vec);
  double contents = - multiVtV(buf_vec, pt_relative_v) / 2.0;

  //calcurate H
  for(int i = 0; i < 3; i++){
	for (int j = 0; j < 3; j++){
	  double jacov_i[2] = {jacov[i][0], jacov[i][1]};
	  double jacov_j[2] = {jacov[j][0], jacov[j][1]};
	  double gjacov[2];
	  estimateGradJacov(gjacov, pt_relative, transform_param[2], i, j);

	  //pt_relative * cov_inv_mat
	  double buf_vec2[2]; 
	  multiMV(cov_inv_mat, jacov_i, buf_vec2);
	  mat_h[i][j] = ( ( multiVtV(buf_vec, jacov_i) ) * ( -multiVtV(buf_vec, jacov_j) )
					  - multiVtV(buf_vec, gjacov) //n=m=3以外では0
					  - multiVtV(jacov_j, buf_vec2) ) * exp(contents);
	}
  }

  //debug用表示
  std::cerr << ">>>hesse calc..." << std::endl;
  std::cerr << "pt: (" << pt.x_pos << ", " << pt.y_pos << ")" << std::endl;
  std::cerr << "pt_trans: (" << pt_transformed.x_pos << ", " << pt_transformed.y_pos << ")" << std::endl;
  std::cerr << "pt_relative: (" << pt_relative.x_pos << ", " << pt_relative.y_pos << ")" << std::endl;
  std::cerr << "mean pt: (" << mean_pt.x_pos << ", " << mean_pt.y_pos << ")" << std::endl;
  std::cerr << "jacov[0]: (" << jacov[0][0] << ", " << jacov[0][1] << ")" << std::endl;
  std::cerr << "jacov[1]: (" << jacov[1][0] << ", " << jacov[1][1] << ")" << std::endl;
  std::cerr << "jacov[2]: (" << jacov[2][0] << ", " << jacov[2][1] << ")" << std::endl;
  std::cerr << "exp content:" << contents << std::endl;

}

bool determineDefinite(double mat[3][3], double rambda)
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
		  mat[0][0] += rambda;
		  mat[1][1] += rambda;
		  mat[2][2] += rambda;
		  // determineDefinite(mat, rambda);
		  return false;
		}
	  else{
		cerr << "正定ではない(3rd stage false)" << endl;
		mat[0][0] += rambda;
		mat[1][1] += rambda;
		mat[2][2] += rambda;
		// determineDefinite(mat, rambda);
		return false;
	  }
	else{
	  cerr << "正定ではない(2nd stage false)" << endl;
	  mat[0][0] += rambda;
	  mat[1][1] += rambda;
	  mat[2][2] += rambda;
	  // determineDefinite(mat, rambda);
	  return false;
	}
  }
  else{
	cerr << "正定ではない(1st stage false)" << endl;
	mat[0][0] += rambda;
	mat[1][1] += rambda;
	mat[2][2] += rambda;
	// determineDefinite(mat, rambda);
	return false;
  }
}

void estimateTransform(std::vector<points>& target_local, std::vector<points>& input_local, double param[DIM+1], double tolerance, int max_cicle, double score)
{
	//平均を求めるお(target)
	points mean_pt;
	estimateMean(mean_pt, target_local);
	cerr << "mean pt: " << mean_pt.x_pos << ", " << mean_pt.y_pos << std::endl;

	//共分散行列を求めるお(target)
	double mat[DIM][DIM];
	estimateCovarianceMat(mat, target_local, mean_pt);
	cerr << "cov mat:" << endl;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			cerr << " " << mat[i][j];
		}
		cerr << endl;
	}

	//共分散行列の逆行列を求めるお(target)
	double inv_cov_mat[DIM][DIM];
	estimateInvMat(inv_cov_mat, mat);
	cerr << "inv cov mat:" << endl;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			cerr << " " << inv_cov_mat[i][j];
		}
		cerr << endl;
	}

	double g_local[2];
	double h_local[3][3];

	if (input_local.size() <= 3){
	  cerr << " entry few points cource..." << endl;
	  for (int i = 0; i < DIM+1; i++)
		param[i] = 999;
	  score = -999;
	}

	else {
	  cerr << " entry normal cource..." << endl;
	  for(int i = 0; i < max_cicle; i++){
		std::cerr << i << "回目のループ" << std::endl;
 		//gとhを求める
		for (int index = 0; index < input_local.size(); index++){		
		  //estimate g
		  double v_g[2];
		  estimateVecG(input_local.at(index), param, mean_pt, inv_cov_mat, v_g);
		  //add g
		  g_local[0] += v_g[0];
		  g_local[1] += v_g[1];
		  
		  //estimate H
		  double mat_h[3][3];
		  estimateMatH(input_local.at(index), param, mean_pt, inv_cov_mat, mat_h);
		  //add H
		  for(int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++)
			  h_local[i][j] += mat_h[i][j];
		  }
		}
		
		if(i == 0){
		  // cerr << "hessian:" << endl;
		  // for(int i = 0; i < 3; i++){
		  // 	for(int j = 0; j < 3; j++){
		  // 	  cerr << " " << h_local[i][j];
		  // 	}
		  // 	cerr << endl;
		  // }
		  determineDefinite(h_local, 1.0);
		}
		
		//show g and H
		std::cerr << "vec_g:" << g_local[0] << " " << g_local[1] << std::endl;
		std::cerr << "hessian:" << std::endl;
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
			  std::cerr << " " << h_local[i][j];
			}
			std::cerr << std::endl;
		}
		
		// if(-tolerance < g_local[0] && g_local[0] <tolerance)
		//   if(-tolerance < g_local[1] && g_local[1] <tolerance){
		// 	std::cerr << "tolerance clear g = (" << g_local[0] << ", " << g_local[1] << ")" << std::endl;
		// 	break;
		//   }
		
		double inv_h[3][3];
		estimateInvMat3d(h_local, inv_h);
		//bool definite = determineDefinite(inv_h);
			
		//変化量を求めるニキ
		multiMV3d(inv_h, g_local, param);
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
