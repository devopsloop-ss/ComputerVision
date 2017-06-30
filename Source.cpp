/*The code here extracts HOG descriptor vector from an image of given size. The image pixel data has to be extracted externally
by some other program and has to be stored as a csv file. The various outputs are written into various csv files.
Parameters can be varied using given macros.
Note that if bin size is changed then suitable values have to be added to bin values array*/


#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>

#define HEIGHT 100
#define WIDTH 100
#define C_SIZE 5
#define BIN_SIZE 5
#define BLK_SIZE 2
#define blkHistSize BIN_SIZE*BLK_SIZE*BLK_SIZE

int **image, **output;

using namespace std;

struct grad
{
	float x, y;
	float mag,angle;
};

class gradient//class for calculating and storing gradients
{
public:
	int i, j;
	grad **gdt;

	gradient(int h, int w);
	~gradient()
	{
		delete gdt;
	}
};

gradient::gradient(int h, int w)
{
	gdt = new grad*[h];

	for (i = 0; i < h; i++)
	{
		gdt[i] = new grad[w];
		for (j = 0; j < w; j++)
		{
			if (i != 0 && j != 0 && i != (h - 1) && j != (w - 1))
			{
				gdt[i][j].x = (float)image[i][j - 1] - image[i][j + 1];
				gdt[i][j].y = (float)image[i + 1][j] - image[i - 1][j];
			}
			else
			{
				gdt[i][j].x = 0.0;
				gdt[i][j].y = 0.0;
			}

			if (gdt[i][j].x == 0.0)
			{
				if (gdt[i][j].y != 0.0)
					gdt[i][j].angle = 3.14159265f / 2;
				else
					gdt[i][j].angle = 5.0;
			}

			else
			gdt[i][j].angle = atan(gdt[i][j].y / gdt[i][j].x);
			gdt[i][j].mag = sqrt(gdt[i][j].x*gdt[i][j].x + gdt[i][j].y*gdt[i][j].y);
			if (gdt[i][j].angle < 0)
			{
				gdt[i][j].angle += 3.14159265f;
			}
			if(gdt[i][j].mag!=0.0)
			{
				gdt[i][j].x /= gdt[i][j].mag;
				gdt[i][j].y /= gdt[i][j].mag;
			}
			gdt[i][j].mag = (gdt[i][j].mag != 0.0)?1.0:0.0;
		}
	}
}

class cell
{
public:
	int starti, startj, i, j, imax, jmax;
	float histogram[BIN_SIZE], bin_cent[BIN_SIZE];
	void setup(int a, int b, gradient &g);
	
};

void cell::setup(int a, int b, gradient &g)
{
	starti = a;
	startj = b;
	bin_cent[0] = 0.0;
	bin_cent[1] = 0.78539816f;
	bin_cent[2] = 1.57079633f;
	bin_cent[3] = 2.35619449f;
	bin_cent[4] = 3.14159265f;
	imax = a + C_SIZE;
	jmax = b + C_SIZE;
	for (i = 0; i < BIN_SIZE; i++)
		histogram[i] = 0.0;
	for (i = starti; i < imax; i++)
	{
		for (j = startj; j < jmax; j++)
		{
			if (g.gdt[i][j].angle == 5.0)
				continue;
			int k, l;
			for (k = 4, l = 3; l > -1; k--, l--)
			{
				if (g.gdt[i][j].angle >= bin_cent[l])
					break;
			}
			float q = bin_cent[k] - g.gdt[i][j].angle, r = g.gdt[i][j].angle - bin_cent[l];
			histogram[k] += g.gdt[i][j].mag*(r / (q + r));
			histogram[l] += g.gdt[i][j].mag*(q / (q + r));
		}
	}
	cout << "Histogram : ";
	for (int k = 0; k < BIN_SIZE; k++)
	{
		cout << histogram[k] << ',';
	}
	cout << endl;
}

class block
{
public:
	float block_his[blkHistSize];
	void setup(int,int, cell**);
};

void block::setup(int a, int b, cell **c)
{
//	block_his[] = c[a][b].histogram[]+c[a+1][b].histogram[]+c[a+1][b+1].histogram[]+c[a][b+1].histogram[]
		
	for (int i = a,k=0; i <= a+(BLK_SIZE-1) ; i++)
	{
		for (int j = b; j <= b+(BLK_SIZE-1); j++)
		{
			for (int l = 0; l < BIN_SIZE; l++)
				block_his[k++] = c[i][j].histogram[l];
		}
	}
	//Normalisation
	float mag=0.0;
	for (int k = 0; k < blkHistSize; k++)
	{
		mag = mag + block_his[k]*block_his[k];
	}
	mag = sqrt(mag += 1.0);
	
	for (int k = 0; k < blkHistSize; k++)
	{
		block_his[k] /= mag;
	}

	cout << "\nBlock normalised Histogram : ";
	for (int k = 0; k < blkHistSize; k++)
	{
		cout << block_his[k] << ',';
		if ((k+1)%BIN_SIZE == 0)
			cout << "\t";
	}
	cout << endl;
}
int main()
{
	char path[200];
	ifstream filein;
	ofstream gmag, g_x, g_y, g_ang, HogV, histogram;
	cout << "Enter the path of the image : ";
	cin >> path;
	filein.open(path);
	gmag.open("gmag.csv");
	g_x.open("g_x.csv");
	g_y.open("g_y.csv");
	g_ang.open("g_ang.csv");
	HogV.open("HogV.csv");
	histogram.open("h.csv");
	string x;
	float per=0;
	image = new int*[HEIGHT];
	output = new int*[HEIGHT];
	for (int j = 0; j < HEIGHT; j++)
	{
		image[j] = new int[WIDTH];
		output[j] = new int[WIDTH];
	}
	cout << "Loading image : " << endl;
	for (int i = 0; i<HEIGHT; i++)//Reading input from image
	{

		if (!filein.good())
		{
			cout << "Error opening file\n";
			return -1;
		}
		getline(filein, x);
		stringstream str(x);

		for (int j = 0; j<WIDTH; j++)
		{
			string val;
			if (!str.good())
				break;
			getline(str, val, ',');
			stringstream trn(val);
			trn >> image[i][j];
			
		}
		int l = HEIGHT / 100;
		if ((i+1) % l == 0)
		{
			per = ((float)(i+1)/HEIGHT)*100;
			cout << "\r"<<per<<"% loaded";
			if (per == 100.00)
			{
				cout << endl;
			}
		}
	}

	for (int i = 0; i<HEIGHT; i++)//Initialising output Array
	{
		for (int j = 0; j<WIDTH; j++)
		{
			output[i][j] = 255;
		}
	}
	
	//Main Code
	gradient gr(HEIGHT,WIDTH);
	const int m = HEIGHT / C_SIZE;
	const int n = WIDTH / C_SIZE;
	cell **cl;
	cl = new cell*[m];
	
	for (int i = 0, k=0; i < m; i++,k+=C_SIZE)
	{
		cl[i] = new cell[n];
		for (int j = 0,l=0; j < n; j++, l+=C_SIZE)
		{
			cl[i][j].setup(k, l, gr);
		}
	}

	block blk[m - 1][n - 1];
	for (int i = 0; i < m-1; i++)
	{
		for (int j = 0; j < n - 1; j++)
		{
			blk[i][j].setup(i,j,cl);
		}
	}
	
	const int o = blkHistSize*(m-1)*(n-1);
	float HOGvect[o];
		
	for (int i = 0,k=0; i < m - 1; i++)
	{
		for (int j = 0; j < n - 1; j++)
		{
			for (int l = 0; l < blkHistSize; l++)
			{
				HOGvect[k++] = blk[i][j].block_his[l];
			}
		}
	}
	
	//End of main code

	//output to output files
	for (int i = 0; i<HEIGHT; i++)
	{
		for (int j = 0; j<WIDTH; j++)
		{
			//fileout << output[i][j] << ",";
			gmag << gr.gdt[i][j].mag<<",";
			g_x << gr.gdt[i][j].x << ",";
			g_y << gr.gdt[i][j].y << ",";
			g_ang << gr.gdt[i][j].angle << ",";
		}
		gmag << "\n";
		g_x << "\n";
		g_y << "\n";
		g_ang << "\n";
	}
	float h[BIN_SIZE];
	for (int i = 0; i < BIN_SIZE; i++)
		h[i] = 0.0;
	cout << endl;
	for (int i = 0,j=0; i < o; i++,j++)
	{
		HogV << HOGvect[i] << ",";
		cout << HOGvect[i] << ",";
		h[j] += HOGvect[i];
		if (j + 1 == BIN_SIZE)
		{
			j = -1;
			cout << endl;
		}
	}

	cout << "\nFinal Histogram : \n";

	for (size_t i = 0; i < BIN_SIZE; i++)
	{
		histogram << h[i]<<",";
		cout << h[i] << ",";
	}
	filein.close();
	gmag.close();
	g_x.close();
	g_y.close();
	g_ang.close();
	HogV.close();

	delete image;
	delete output;
	delete cl;

	return 0;
}
