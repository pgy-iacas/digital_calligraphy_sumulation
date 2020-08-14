#include<graphics.h>
#include<conio.h>
#include<iostream>
#include<vector>
#include<time.h>
#include<math.h>
#include<iterator>
using namespace std;
class Papel
{
public:
	Papel(int a, int b, float gray_v)
	{
		this->x = a;
		this->y = b;
		this->gray_v = gray_v;
	}
	int x;         //位置
	int y;         //位置
	float gray_v;    //墨汁浓度
};
class Points_in_distance
{
public:
	Points_in_distance(int distance)
	{
		for (int i = 0; i <= distance; i++)
		{
			for (int j = 0; j <= distance; j++)
			{
				float this_dis = sqrt(i*i + j*j);
				//float this_dis = abs(i) + abs(j);
				if (this_dis <= (float)distance)
				{
					data.push_back(pair<int, int>(i, j));
				}
			}
		}
	}
	const vector<pair<int, int>> get_data(){
		return data;
	}
private:
	vector<pair<int, int>>data;
};
class Markov_matrix
{
private:
	int size, rank;
public:
	Markov_matrix(int rank, int size = 11) :size(size), rank(rank)     //rank==1   返回一次转移矩阵
	{
		for (int i = size - 1; i >= 0; i--)
		{
			vector<float>tmp(size, 0);
			matrix.push_back(tmp);
		}
		for (int i = size - 1; i>0; i--)
		{
			float b = 3 / pow(i, 3);
			for (int j = i; j >0; j--)
			{
				float tmp = b *(pow(j, 3) - pow(j - 1, 3)) / 3;
				matrix[i][j] = tmp;
			}
		}
		vector<vector <float> > copy = matrix;
		vector<vector <float> >matrix_2 = matrix;
		for (int r = 1; r<rank; r++)
		{
			for (int i = 0; i<size; i++)
			{
				for (int j = 0; j<size; j++)
				{
					float tmp = 0;
					for (int p = 0; p<size; p++)
					{
						tmp += copy[i][p] * matrix[p][j];
					}
					matrix_2[i][j] = tmp;
				}
			}
			copy = matrix_2;
		}
		matrix = matrix_2;
		matrix[0][0] = 1;
	}
	vector<vector<float>> get_matrx(){
		return matrix;
	}
	vector<vector<float>>matrix;
};
class Surface
{
private:
	void _draw_two_point(int xc, int yc, int a, int b)
	{
		auto tmp = _point_rotate(xc + a, yc + b, xc, yc);
		edge_vec.push_back(pair<int,int>(tmp.first, tmp.second));
		tmp = _point_rotate(xc - a, yc + b, xc, yc);
		edge_vec.push_back(pair<int,int>(tmp.first, tmp.second));
		//edge_vec.push_back(_point_rotate(xc + a, yc -b, xc, yc));
		//edge_vec.push_back(_point_rotate(xc - a, yc - b, xc, yc));
		for (int i = xc - a; i <= xc + a; i++)
		{
			tmp = _point_rotate(i, yc + b, xc, yc);
			inner_vec.push_back(pair<int,int>(tmp.first, tmp.second));
			//inner_vec.push_back(_point_rotate(i, yc - b, xc, yc));
		}
		return;
	}
	void _draw_eclips_and_trangle(int xc, int yc, int a, int b)
	{
		//画半椭圆
		int sqa = a*a;
		int sqb = b*b;
		int x = 0;
		int y = b;
		int d = 2 * sqb - 2 * b*sqa + sqa;
		_draw_two_point(xc, yc, x, y);
		int P_x = (int)((double)sqa / sqrt((double)(sqa + sqb)));
		while (x <= P_x)
		{
			if (d < 0)
				d += 2 * sqb*(2 * x + 3);
			else
			{
				d += 2 * sqb*(2 * x + 3) - 4 * sqa*(y - 1);
				y--;
			}
			x++;
			_draw_two_point(xc, yc, x, y);
		}
		d = sqb * (x * x + x) + sqa * (y * y - y) - sqa * sqb;
		while (y >= 0)
		{
			_draw_two_point(xc, yc, x, y);
			y--;
			if (d < 0)
			{
				x++;
				d = d - 2 * sqa*y - sqa + 2 * sqb*x + 2 * sqb;
			}
			else
				d = d - 2 * sqa*y - sqa;
		}

		//画三角形
		int height_tri = 2 * a;
		int dx, dy, n, k;
		double xinc, yinc, newx, newy;
		int x1, x2, y1, y2;
		x1 = xc - a;
		y1 = yc;
		x2 = xc;
		y2 = yc - height_tri;
		dx = x2 - x1;
		dy = y2 - y1;
		if (abs(dx) - abs(dy)>0)    //比较两参数的绝对值哪一个大，哪一个就作为步长参数（n），此参数将作为沿直线所画出点的数目
			n = abs(dx);
		else
			n = abs(dy);
		xinc = (double)dx / n;
		yinc = (double)dy / n;
		newx = x1;
		newy = y1;
		for (k = 0; k<n; k++)
		{
			int tmp = floor(newx + 0.5);
			_draw_two_point(xc, yc, xc - floor(newx + 0.5), floor(newy + 0.5) - yc);
			newx += xinc;
			newy += yinc;
		}
	}
	pair<int, int> _point_rotate(const double x, const double y, const double xc, const double yc)
	{
		//xc,yc为rotate中心,x,y为原始旋转点  负角度为逆时针转
		double new_point[2];

		new_point[0] = ((x - xc)*cos(angle) - (y - yc)*sin(angle)) + xc;
		new_point[1] = ((x - xc)*sin(angle) + (y - yc)*cos(angle)) + yc;
		return pair<int, int>((int)new_point[0], (int)new_point[1]);
	}
public:
	Surface(int x, int y, double angle , int strok_size)
	{
		this->x = x;
		this->y = y;
		this->angle = angle / 180 * 3.1415926;  //化为弧度
		a = strok_size;
		b = a / 2;
		_draw_eclips_and_trangle(x, y, a, b); //椭圆点存入向量 
	}
	inline vector<pair<int, int>>& get_edge_vec()
	{
		return edge_vec;
	}
	vector<pair<int, int>>& get_inner_vec()
	{
		return inner_vec;
	}
private:
	int x, y;
	int a ;
	int b ;
	double rotate_matrix[2][2];
	double angle;
	vector<pair<int,int>>edge_vec;
	vector<pair<int, int>>inner_vec;
};
class INK
{
private:
	inline int intense2gray(float intense)
	{
		return 255 - (int)(intense * 255);
	}
	inline float gray2intense(int grayval)
	{
		return (float)(255 - grayval) / 255;
	}
	inline void set_color(int x,int y,float tense)
	{
		int v = intense2gray(tense);
		pMem[y*length + x] = RGB(v, v, v);
	}
	inline void set_disperse_color(int xc,int yc,int x1,int y1,float tense)
	{
		int v = intense2gray(tense);
		pMem[(yc+y1)*length + xc+x1] = RGB(v, v, v);
		pMem[(yc - y1)*length + xc + x1] = RGB(v, v, v);
		pMem[(yc + y1)*length + xc - x1] = RGB(v, v, v);
		pMem[(yc - y1)*length + xc - x1] = RGB(v, v, v);
	}
	void _generate_vec( Surface tmpsurface)
	{
		edg_vec = tmpsurface.get_edge_vec();
		inner_vec = tmpsurface.get_inner_vec();
		return;
	}
	void _draw_inner_shape(DWORD* pMem, const vector<pair<int,int>> &vec)
	{
		for (auto it = vec.begin(); it != vec.end(); it++)
		{
			int x = it->first;
			int y = it->second;
			set_color(x, y, intensity);
		}
		return;
	}
	void _draw_edge_shape(DWORD* pMem, const vector<pair<int, int>> &vec)
	{
		for (auto it = vec.begin(); it != vec.end(); it++)
		{
			int x = it->first;
			int y = it->second;
			set_color(x, y, intensity);
		}
		return;
	}
	void _draw_dispersion(DWORD* pMem, const vector<pair<int, int>> &vec, int rank)
	{
		for (auto it1 = vec.begin(); it1 != vec.end(); it1++)
		{
			int ran = rand();
			if (ran%2==0)
				continue;
			for (int i = 0; i < rank; i++)
			{
				for (auto it2 = point_in_distance_set[i].begin(); it2 != point_in_distance_set[i].end(); it2++)
				{
					set_disperse_color(it1->first, it1->second, it2->first, it2->second, intensity);
				}
			}
		}
	}
public:
	INK(int L, int W, float intensity, float stoke_angle, int strok_size, int disperse_rank, 
		const vector<vector<vector<float>>>&matrix_set,
		const vector<vector<pair<int, int>>>&point_in_distance_set)
		:length(L), width(W), intensity(intensity), stroke_angle(stoke_angle),strok_size(strok_size), disperse_rank(disperse_rank)
		, matrix_set(matrix_set)
		, point_in_distance_set(point_in_distance_set)
	{
		;
	}
	const void draw_strok(const MOUSEMSG &m)
	{
		pMem = GetImageBuffer();
		Surface tmpsurface(m.x, m.y, stroke_angle, strok_size);
		_generate_vec(tmpsurface);    //将位置和像素放入edg_vec和inner_vec
		_draw_inner_shape(pMem, inner_vec);
		_draw_edge_shape(pMem, edg_vec);       //遍历vector的值并更改缓冲区
		_draw_dispersion(pMem, edg_vec, disperse_rank);
		return;
	}
	inline void change_strok_size(int size)
	{
		if (size > 0)
			strok_size += 1;
		else if (strok_size != 0)
			strok_size -= 1;
	}
private:
	int length;     //长度
	int width;      //宽度
	float intensity;  //浓度 黑色为1
	int strok_size;    // 笔触大小 
	int disperse_rank; // 扩散程度   正整数
	float stroke_angle; //笔触角度
	vector<pair<int,int>>edg_vec;
	vector<pair<int,int>>inner_vec;
	vector<vector<vector<float>>>matrix_set;
	vector<vector<pair<int, int>>>point_in_distance_set;
	DWORD* pMem;
};

int main()
{
	const int Length = 1500;
	const int Width = 800;
	const float ink_tense = 1;
	const int strok_size = 15;
	const int disperse_rank = 3;     //扩散强度 取正整数
	const float strok_angle = -55;
	vector<vector<vector<float>>>matrix_set;
	vector<vector<pair<int, int>>>point_in_distance_set;
	for (int i = 1; i <= disperse_rank; i++)           //rank==1  一次转移矩阵
	{
		Markov_matrix markov(i);
		matrix_set.push_back(markov.get_matrx());
		Points_in_distance tmp_rank(i);
		point_in_distance_set.push_back(tmp_rank.get_data());
	}
	initgraph(Length, Width);
	setbkcolor(RGB(255, 255, 255));
	cleardevice();
	//setfillcolor(RGB(0,0,0));
	//const int width = 6;
	MOUSEMSG m;		// 定义鼠标消息
	INK paper(Length, Width, ink_tense, strok_angle, strok_size, disperse_rank, matrix_set, point_in_distance_set);
	int count = 0;
	while (true)
	{
		count++;
		if (count % 100000 != 0)
			continue;
		// 获取一条鼠标消息
		m = GetMouseMsg();
		if ((m.uMsg == WM_MOUSEMOVE&&m.mkLButton == true) || m.uMsg == WM_LBUTTONDOWN)
		{
			paper.draw_strok(m);
		}
		if (m.uMsg == WM_RBUTTONDOWN)
			cleardevice();
		if (m.uMsg == WM_MBUTTONDOWN)
			break;
		if (m.wheel != 0)
			paper.change_strok_size(m.wheel);
	}
	saveimage("ink.bmp");
	// 关闭图形窗口
	closegraph();
}
