/***********************************************************************
*  This code is part of NewPolyVest.
*
***********************************************************************/

#define PI 3.1415926536

#include <stdio.h>
#include <fstream>
#include "vol.h"
#include <time.h>

using namespace vol;

ifstream file;

double get_num(){
	char ch = ' ';
	double tmp_num = 0;
	bool isnumber = false;
	bool tmp_neg = false;
	bool decimal = false;
	int dec_len = 0;
	for (; (!file.eof()) && (isspace(ch)); file.get(ch));
	while (!file.eof()){
		if (isdigit(ch)){
			if (decimal) dec_len++;
			tmp_num = tmp_num * 10 + ch - 48;
			isnumber = true;
		}else if (ch == '-') tmp_neg = true;
		else if (ch == '.')	decimal = true;
		else if (isspace(ch)) break;
		else{
			cout << "error: unknown identifier \'" << ch << "\'." << endl;
			exit(1);
		}
		file.get(ch);
	}
	if (!isnumber){
		cout << "error: input file is incomplete, please check the size of matrix." << endl;
		exit(1);
	}
	if (tmp_neg) tmp_num = -tmp_num;
	if (dec_len > 0) tmp_num /= pow((double)10, dec_len);
	return tmp_num;
}

int factorial(int n)
{
	int fac = 1;
	for (int i = 1; i <= n; i++)
	{
		fac *= i;
	}

	return fac;
}

double gamma2(int n) /* gamma( (double) n / 2 ) */
{
	double g = 0;
	if (n % 2 == 0)
	{
		g = 1;
		for (int i = 4; i <= n; i = i + 2)
		{
			g *= ((double)i / 2 - 1);
		}
	}
	else if (n % 2 == 1)
	{
		g = sqrt(PI);

		for (int i = 3; i <= n; i = i + 2)
		{
			g *= ((double)i / 2 - 1);
		}
	}
	return g;
}

double newsin(double theta)
{
	double s = 0;
	for (int i = 1; i <= 5; i++)
	{
		s += pow(-1, i - 1) * pow(theta, 2 * i - 1) / factorial(2 * i - 1);
	}
	return s;
}

double newcos(double theta)
{
	double c = 0;
	for (int i = 0; i <= 4; i++)
	{
		c += pow(-1, i) * pow(theta, 2 * i) / factorial(2 * i);
	}
	return c;
}

rowvec sinn(int n) /* The first ten coefficients of The Taylor expansion of sin^n */
/* the i-th coefficient is the coefficient of x^(2*i+n-2) */
{
	rowvec a = zeros(1, 5);
	rowvec a1 = zeros(1, 5);
	for (int i = 1; i <= 5; i++)
	{
		a1[i - 1] = pow(-1, i - 1) / factorial(2 * i - 1);
	}
	if (n == 1)
	{
		a = a1;
	}
	else
	{
		rowvec am = sinn(n - 1);
		for (int i = 1; i <= 5; i++)
		{
			for (int j = 1; j <= 6 - i; j++)
			{
				a[i + j - 2] += am[i - 1] * a1[j - 1];
			}
		}
	}
	return a;
}

double newfun(double theta, int k) /* k>=2 */
/* an estimate of sin^(k-1)(theta)/int(sin^(k-2)(alpha),0,tehta) */
{
	double s = 0;
	if (k == 2)
	{
		s = newsin(theta) / theta;
	}
	else
	{
		rowvec b = sinn(k - 1);
		rowvec a = sinn(k - 2);
		for (int i = 1; i <= 5; i++)
		{
			a[i - 1] /= ((double)2 * i + k - 3);
		}
		rowvec c = zeros(1, 5);
		c[0] = b[0] / a[0];
		for (int i = 2; i <= 5; i++)
		{
			double sum = 0;
			for (int j = 1; j <= i - 1; j++)
			{
				sum += c[j - 1] * a[i - j];
			}
			c[i - 1] = (b[i - 1] - sum) / a[0];
		}
		for (int i = 1; i <= 5; i++)
		{
			s += c[i - 1] * pow(theta, 2 * i - 2);
		}
	}

	return s;
}

double rateofvol(double rate, int n, double theta)
{
	double p, e = 0;
	p = pow(rate, 2) - pow(newsin(theta), 2);
	e = pow(p, (double)n / 2) / pow(rate, n) / newcos(theta) / pow(newsin(theta) * sqrt(pow(rate, -2) - 1) + newcos(theta), n - 1);

	return e;
}

double solverate(double rate, int n, double epsi)/* find a proper \theta such that rateofvol(rate,n,\theta)>=epsi*/
{
	double theta, lb = 0, ub = rate, mb;
	double lv, uv, mv;
	lv = rateofvol(rate, n, lb) - epsi;
	uv = rateofvol(rate, n, ub) - epsi;
	while (1)
	{
		if (uv >= 0)
		{
			theta = ub;
			break;
		}
		else
		{
			mb = (lb + ub) / 2;
			mv = rateofvol(rate, n, mb) - epsi;
			if (lv * mv < 0)
			{
				ub = mb;
				uv = mv;
			}
			else if (lv * mv > 0)
			{
				lb = mb;
				lv = mv;
			}
			else if (mv == 0)
			{
				theta = mb;
				break;
			}
			if (abs(lv - uv) < pow(10, -8))
			{
				theta = lb;
				break;
			}
		}
	}


	return theta;
}
//-----------------------thread1: gettheta---------------------------//

polytope* p_thread;

typedef struct
{
	int samplenum;
	double rate_i;
}Rate;

pthread_mutex_t rate_mutex;

double rate = 1;

void* rate_fn(void *arg)
{
	Rate *r =(Rate*)arg;
	int count;
	for (count=1;count<=r->samplenum;count++){
		rowvec sampledata = p_thread->gaussrand();
		double distance = sampledata[0];
		double height = sampledata[1];
		if (r->rate_i > height/distance){
			r->rate_i = height/distance;
		}
	}

	pthread_mutex_lock(&rate_mutex);
	if (rate> r->rate_i){
		rate = r->rate_i;
	}
        pthread_mutex_unlock(&rate_mutex);

	return (void*)0;
}

double Gettheta(int coef,double bound, int num_thread)
{
        pthread_t* tid = (pthread_t*)malloc(num_thread * sizeof(pthread_t));
        Rate* Rate_i = (Rate*)malloc(num_thread * sizeof(Rate));

        pthread_mutex_init(&rate_mutex, NULL);

	int i;
	for (i=0; i< num_thread; i++){
		Rate_i[i] = {coef,1};
		pthread_create(&tid[i],NULL,rate_fn,&Rate_i[i]);
	}

	for(i=0; i< num_thread; i++){
		pthread_join(tid[i],NULL);
	}

	pthread_mutex_destroy(&rate_mutex);

	//free(tid);
	//free(Rate_i);

	double theta; 
	theta = solverate(rate, p_thread->n, 1 - bound);

	return theta;
}
//-----------------------thread2: NewEstimate------------------------//

double sum_total = 0;

typedef struct
{
	int samplenum;
	double theta;
}Sum;

pthread_mutex_t sum_mutex;

void* sum_fn(void *arg)
{
	Sum *r = (Sum*)arg;
	double sum_i =0 ;
	int count;
	for(count=1;count<=r->samplenum;count++){
		rowvec sampledata = p_thread->gaussrand();
		double distance = sampledata[0];
		double height = sampledata[1];
		double p = pow(height, 2) - pow(distance, 2) *
		       	pow(newsin(r->theta), 2);
		sum_i += pow(height, p_thread->n) *
		       	pow(distance, p_thread->n)
		       	/ pow(p, (double)p_thread->n / 2);
	}

	pthread_mutex_lock(&sum_mutex);
	sum_total += sum_i;
	pthread_mutex_unlock(&sum_mutex);

	return (void*)0;
}
	
double NewEstimateVol(int coef, double theta, int num_thread)
{
	pthread_t* tid = (pthread_t*)malloc(num_thread * sizeof(pthread_t));
	Sum eachsum = {coef,theta};

	pthread_mutex_init(&sum_mutex,NULL);

	int i;
	for(i=0; i<num_thread; i++){
		pthread_create(&tid[i],NULL,sum_fn,&eachsum);
	}

	for(i=0; i<num_thread; i++){
		pthread_join(tid[i],NULL);
	}

	pthread_mutex_destroy(&sum_mutex);

	//free(tid);
	
	double vol;
	vol = p_thread->determinant *sum_total * newcos(theta) * pow(PI, (double)p_thread->n / 2) / coef /num_thread / gamma2(p_thread->n + 2) / ((double)p_thread->n - 1)
	       	* newfun(theta, p_thread->n);
	return vol;
}

//-----------------------main part-----------------------------------//

int main(int argc, char *argv[]){

	cout << endl << "--------------- PolyVest ----------------" << endl;
	cout << "If you have any questions or if you have found some bugs," << endl << "please contact me at <gecj@ios.ac.cn>." << endl;
	cout << endl << "=========================================" << endl;

	if (argc <= 2 || argc >= 7){
		cout << "error: invalid arguments." << endl;
		cout << "USAGE: " << argv[0] << " <input-file>  <sample-number> <error-bound> <num-thread> [output-file]" << endl;
		return 1;
	}

	file.open(argv[1]);

	if (!file.is_open()){ cout << "Cant open file." << endl; return 1; }

	int rows, cols;
	rows = get_num();
	cols = get_num();
	
	polytope p(rows, cols);

	//cout << "Hyperplanes: " << rows << endl << "Dimensions:" << cols << endl;

	for (int i = 0; i < rows; i++){
		double t = get_num();
		p.vecb(t, i);
		for (int j = 0; j < cols; j++){
			t = get_num();
			p.matA(t, i, j);
		}
	}
	file.close();

	//p.msg_off = true;
	p.check_planes_off = true;
	cout << endl << "============= Preprocessing =============" << endl;
	p.Preprocess();
	p.Centering();
	
	cout << endl << "=============== Sampling ================" << endl;

        p_thread = &p;
	double theta = Gettheta(atoi(argv[2]), atof(argv[3]),atoi(argv[4]));

	printf("angle= %e \n", theta);

	double volume = NewEstimateVol(atoi(argv[2]), theta,atoi(argv[4]));

	cout << endl <<"Volume of polytope (estimate): " << volume <<endl;

	ofstream outfile;
	if (argc == 5){
		outfile.open("NewPolyVest.result", ios::app);
	}else{
		outfile.open(argv[4]);
	}
//	outfile << volume << endl;
	outfile.close();

	return 0;
}
