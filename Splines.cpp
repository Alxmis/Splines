#define _CRT_SECURE_NO_WARNINGS
#define EPS 0.001
#define H 1e-4
#define Jsize 2
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void read_file(FILE* fp, int size, float* array)
{
	fseek(fp, 0, SEEK_SET);
	for (int i = 0; i < size; ++i)
	{
		fscanf(fp, "%f", &array[i]);
	}
}

int count_number(FILE* fp)
{
	fseek(fp, 0, SEEK_SET);
	int counter = 0;
	while (1)
	{
		int value;
		if (fscanf(fp, "%f", &value) == 1)
			counter++;
		if (feof(fp))
			break;
	}
	return counter;
}


float* solve_slae(float** M, float* B, int n)
{
	float** coefficient = malloc(n * sizeof(float*));
	for (int j = 0; j < n; ++j)
	{
		coefficient[j] = malloc(2 * sizeof(float));
	}
	float* root = malloc(n * sizeof(float));
	float* temp_y = malloc(n * sizeof(float));

	for (int i = 0; i < n; ++i)
	{
		if (i == 0)
		{
			temp_y[i] = M[i][i]; 
			coefficient[i][0] = -M[i][1] / temp_y[i];
			coefficient[i][1] = B[i] / temp_y[i];
			continue;
		}
		else if (i == n - 1)
		{
			temp_y[i] = M[i][i] + M[i][i - 1] * coefficient[i - 1][0];
			coefficient[i][1] = (B[i] - M[i][i - 1] * coefficient[i - 1][1]) / temp_y[i];
			continue;
		}
		temp_y[i] = M[i][i] + M[i][i - 1] * coefficient[i - 1][0];
		coefficient[i][0] = -M[i][i + 1] / temp_y[i];
		coefficient[i][1] = (B[i] - M[i][i - 1] * coefficient[i - 1][1]) / temp_y[i];
	}
	
	for (int i = n - 1; i >= 0; --i)
	{
		if (i == n - 1)
		{
			root[i] = coefficient[i][1];
			continue;
		}
		root[i] = coefficient[i][0] * root[i + 1] + coefficient[i][1];
	}
	
	

	for (int j = 0; j < n; ++j)
		free(coefficient[j]);
	free(coefficient);
	free(temp_y);

	return root;
}

// Create parametric spline
float create_p_spline(float* T, float* P, float* Root_derivative, float t_value, int index)
{
	float w = (t_value - T[index]) / (T[index + 1] - T[index]);
	return (1 - w) * P[index] + w * P[index + 1] +
		((-2 * w + 3 * pow(w, 2) - pow(w, 3)) * Root_derivative[index] + (-w + pow(w, 3)) * Root_derivative[index + 1]) *
		(pow((T[index + 1] - T[index]), 2) / 6.0);
}

void create_matrixes_p_spline(float* T, float** D, float* P, float* der_B, float* Root_derivative, float* X, float* Y, char axis, int n)
{
	//Finding parametr t
	T[0] = 0.0;
	float d;
	for (int i = 1; i < n; ++i)
	{
		d = X[i] - X[i - 1];
		T[i] = T[i - 1] + d;
	}

	//Finding second derivative in n-2 points
	float d_0 = 0, d_last = 0;
	for (int i = 0; i < n - 2; ++i)
	{
		for (int j = 0; j < n - 2; ++j)
		{
			if (j == i - 1)
				D[i][j] = T[i + 1] - T[i];
			else if (j == i)
				D[i][j] = 2 * (T[i + 2] - T[i]);
			else if (j == i + 1)
				D[i][j] = T[i + 2] - T[i + 1];
			else
				D[i][j] = 0.0;
			
		}
	}

	//Finding radius-vector of each point
	for (int i = 0; i < n; ++i)
	{
		if (axis == 'x')
			P[i] = X[i];
		else if (axis == 'y')
			P[i] = Y[i];
		else
		{
			printf("Wrong arguments in create_matrixes_p_spline function.");
			exit(1);
		}
	}
	//Filling right table for finding second derivatives
	for (int i = 0; i < n - 2; ++i)
	{
		der_B[i] = 6 * ((P[i + 2] - P[i + 1]) / (T[i + 2] - T[i + 1])) -
			6 * ((P[i + 1] - P[i]) / (T[i + 1] - T[i]));
	}
	//Finding second derivative
	float* temp_array = solve_slae(D, der_B, n - 2);
	Root_derivative[0] = d_0;
	for (int i = 1; i < n - 1; ++i)
	{
		Root_derivative[i] = temp_array[i - 1];
	}
	Root_derivative[n - 1] = d_last;
	free(temp_array);



}

void x_to_parameter(float x_value, float* t_value, float* segment_value, float* T, float* P, float* Root_derivative_x, float* X, int n)
{
	//Finding value of any point of parametric spline
	for (int segment = 0; segment < n - 1; ++segment)
	{
		for (float t = x_value - 10; t < x_value + 10; t += 0.01)
		{
			if (X[segment] <= x_value && x_value <= X[segment + 1])
			{
				if (fabs(create_p_spline(T, P, Root_derivative_x, t, segment) - x_value) <= 0.001)
				{
					*t_value = t;
					*segment_value = segment;
					break;
				}
			}
			else
				break;
		}
		if (*t_value != -1000)
			break;
	}
	if (*t_value == -1000)
	{
		printf("@x_to_parameter: Too small range of 't'. Call the developer or widen the range of 't'.");
		exit(1);
	}
}

float find_dfd(float* T, float* P, float* Root_derivative, float t_value, int index)
{
	return (create_p_spline(T, P, Root_derivative, t_value + H, index) - create_p_spline(T, P, Root_derivative, t_value - H, index)) / 
		(2.0 * H);
}


int get_segment(float* t_in_x, float parameter, int n)
{
	int segment = -1000;
	for (int i = 0; i < n - 1; ++i)
	{
		if (t_in_x[i] <= parameter && parameter <= t_in_x[i + 1])
		{
			segment = i;
			break;
		}
	}
	if (segment != -1000)
		return segment;
	else
	{
		return -100000;
	}
}

int get_segment_dist(float* t_in_x, float parameter, int n)
{
	int segment = -1000;
	for (int i = 0; i < n - 1; ++i)
	{
		if (t_in_x[i] <= parameter && parameter <= t_in_x[i + 1])
		{
			segment = i;
			break;
		}
	}
	if (segment != -1000)
		return segment;
	else
	{
		printf("@get_segment_dist: Error.");
		exit(1);
	}
}



void find_reversed_Jacobian_value(float* T1, float* T2, float* P_x1, float* P_x2, float* P_y1, float* P_y2, 
	float* Root_derivative_x1, float* Root_derivative_x2, float* Root_derivative_y1, float* Root_derivative_y2, 
	int index1, int index2, float* inter_point)
{
	float J[Jsize][Jsize];
	
	//Creating identity matrix
	float I[Jsize][Jsize];
	for (int i = 0; i < Jsize; ++i)
	{
		for (int j = 0; j < Jsize; ++j)
		{
			if (i == j)
				I[i][j] = 1;
			else
				I[i][j] = 0;
		}
	}

	
	
	//Creating common Jacobian
	for (int i = 0; i < Jsize; ++i)
	{
		for (int j = 0; j < Jsize; ++j)
		{
			if (i == 0)
			{
				if (j == 0)
					J[i][j] = find_dfd(T1, P_x1, Root_derivative_x1, inter_point[0], index1);
				else if (j == 1)
					J[i][j] = -find_dfd(T2, P_x2, Root_derivative_x2, inter_point[1], index2);
			}
			else if (i == 1)
			{
				if (j == 0)
					J[i][j] = find_dfd(T1, P_y1, Root_derivative_y1, inter_point[0], index1);
				else if (j == 1)
					J[i][j] = -find_dfd(T2, P_y2, Root_derivative_y2, inter_point[1], index2);
			}
		}
	}

	//Creating reversed Jacobian
	float r, res;
	for (int k = 0; k < Jsize; ++k)
	{
		r = 1 / J[k][k];
		for (int j = 0; j < Jsize; ++j)
		{
			J[k][j] *= r;
			I[k][j] *= r;
		}
		for (int i = k + 1; i < Jsize; ++i)
		{
			res = J[i][k];
			for (int z = 0; z < Jsize; ++z)
			{
				J[i][z] = J[i][z] - J[k][z] * res;
				I[i][z] = I[i][z] - I[k][z] * res;
			}
		}
	}
	for (int k = Jsize - 1; k >= 0; --k)
	{
		for (int i = k - 1; i >= 0; --i)
		{
			res = J[i][k];
			for (int z = Jsize - 1; z >= 0; --z)
			{
				J[i][z] = J[i][z] - J[k][z] * res;
				I[i][z] = I[i][z] - I[k][z] * res;	
			}
		}
	}

	
	//Getting value
	float* temp_inter_point = malloc(Jsize * sizeof(float));
	temp_inter_point[0] = inter_point[0];
	temp_inter_point[1] = inter_point[1];
	for (int i = 0; i < Jsize; ++i)
	{
		inter_point[i] = 0;
		for (int j = 0; j < Jsize; ++j)
		{
			if (j == 0)
				inter_point[i] += I[i][j] * (create_p_spline(T1, P_x1, Root_derivative_x1, temp_inter_point[0], index1) -
					create_p_spline(T2, P_x2, Root_derivative_x2, temp_inter_point[1], index2));
			else if (j == 1)
				inter_point[i] += I[i][j] * (create_p_spline(T1, P_y1, Root_derivative_y1, temp_inter_point[0], index1) -
					create_p_spline(T2, P_y2, Root_derivative_y2, temp_inter_point[1], index2));
		}
		inter_point[i] = temp_inter_point[i] - inter_point[i];
	}

}



int main()
{
	//FILLING ARRAYS
	FILE* fp1;
	if ((fp1 = fopen("X1.txt", "rb")) == NULL)
	{
		printf("Error while opening file.\n");
		exit(1);
	}
	FILE* fp2;
	if ((fp2 = fopen("Y1.txt", "rb")) == NULL)
	{
		printf("Error while opening file.\n");
		exit(1);
	}
	FILE* fp3;
	if ((fp3 = fopen("X2.txt", "rb")) == NULL)
	{
		printf("Error while opening file.\n");
		exit(1);
	}
	FILE* fp4;
	if ((fp4 = fopen("Y2.txt", "rb")) == NULL)
	{
		printf("Error while opening file.\n");
		exit(1);
	}

	int n = count_number(fp1);
	if (n != count_number(fp2) || n != count_number(fp3) || n != count_number(fp4))
	{
		printf("Different points number");
		exit(1);
	}
	float* X1 = malloc(n * sizeof(float));
	float* Y1 = malloc(n * sizeof(float));
	float* X2 = malloc(n * sizeof(float));
	float* Y2 = malloc(n * sizeof(float));

	read_file(fp1, n, X1);
	read_file(fp2, n, Y1);
	read_file(fp3, n, X2);
	read_file(fp4, n, Y2);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	

	//REALIZING PARAMETRIZED SPLINE
	//CREATING ARRAYS AND MATRIXES
	//Finding t parameter
	float* T1 = malloc((n - 1) * sizeof(float));
	//Finding second derivative in n-2 points
	float** D1 = malloc((n - 2) * sizeof(float*));
	for (int section = 0; section < n - 2; ++section)
	{
		D1[section] = malloc((n - 2) * sizeof(float));
	}
	//Finding radius-vector_x of each point
	float* P_x1 = malloc(n * sizeof(float));
	//Finding radius-vector_y of each point
	float* P_y1 = malloc(n * sizeof(float));
	// Creating right table for finding second derivatives x
	float* der_B_x1 = malloc((n - 2) * sizeof(float));
	// Creating right table for finding second derivatives y
	float* der_B_y1 = malloc((n - 2) * sizeof(float));
	//Finding second derivative x
	float* Root_derivative_x1 = malloc(n * sizeof(float));
	//Finding second derivative y
	float* Root_derivative_y1 = malloc(n * sizeof(float));
	create_matrixes_p_spline(T1, D1, P_x1, der_B_x1, Root_derivative_x1, X1, Y1, 'x', n);
	create_matrixes_p_spline(T1, D1, P_y1, der_B_y1, Root_derivative_y1, X1, Y1, 'y', n);
	//Finding t parametr
	float* T2 = malloc((n - 1) * sizeof(float));
	//Finding second derivative in n-2 points
	float** D2 = malloc((n - 2) * sizeof(float*));
	for (int section = 0; section < n - 2; ++section)
	{
		D2[section] = malloc((n - 2) * sizeof(float));
	}
	//Finding radius-vector_x of each point
	float* P_x2 = malloc(n * sizeof(float));
	//Finding radius-vector_y of each point
	float* P_y2 = malloc(n * sizeof(float));
	// Creating right table for finding second derivatives x
	float* der_B_x2 = malloc((n - 2) * sizeof(float));
	// Creating right table for finding second derivatives y
	float* der_B_y2 = malloc((n - 2) * sizeof(float));
	//Finding second derivative x
	float* Root_derivative_x2 = malloc(n * sizeof(float));
	//Finding second derivative y
	float* Root_derivative_y2 = malloc(n * sizeof(float));
	create_matrixes_p_spline(T2, D2, P_x2, der_B_x2, Root_derivative_x2, X2, Y2, 'x', n);
	create_matrixes_p_spline(T2, D2, P_y2, der_B_y2, Root_derivative_y2, X2, Y2, 'y', n);

	
	//FINDING VALUE OF ANY POINT OF PARAMETRIZED SPLINE	
	float x_value; // !!! Point x for finding y is expected to be entered here
	printf("Enter x coordinate of THE FIRST spline to get its parameter: ");
	scanf("%f", &x_value);
	float t_value1 = -1000;
	float segment_value1 = -1;
	x_to_parameter(x_value, &t_value1, &segment_value1, T1, P_x1, Root_derivative_x1, X1, n);
	float t_value2 = -1000;
	float segment_value2 = -1;
	printf("Parameter value of spline in the point = %f\n", t_value1);
	printf("x = %f y = %f\n\n", x_value, create_p_spline(T1, P_y1, Root_derivative_y1, t_value1, segment_value1));

	

	//INTERSECTION
	printf("Intersection point: \n");
	float* t_in_x1 = malloc(n * sizeof(float));
	float* t_in_x2 = malloc(n * sizeof(float));
	float temp_segment_value = -1;
	for (int i = 0; i < n; ++i)
	{
		t_in_x1[i] = -1000;
		t_in_x2[i] = -1000;
		x_to_parameter(X1[i], &t_in_x1[i], &temp_segment_value, T1, P_x1, Root_derivative_x1, X1, n);
		x_to_parameter(X2[i], &t_in_x2[i], &temp_segment_value, T2, P_x2, Root_derivative_x2, X2, n);
	}

	//Determening zero approximation inter_point for t_value1 & t_value2
	float* inter_point = malloc(2 * sizeof(float));
	inter_point[0] = t_in_x1[2]; // Looks this way: inter_point(t_value1, t_value2)
	inter_point[1] = t_in_x2[2];
	segment_value1 = get_segment(t_in_x1, inter_point[0], n);
	segment_value2 = get_segment(t_in_x2, inter_point[1], n);

	//Usage of Newton method
	int key1 = -1;
	int key2 = -1;
	while (fabs(create_p_spline(T1, P_y1, Root_derivative_y1, inter_point[0], segment_value1) -
		create_p_spline(T2, P_y2, Root_derivative_y2, inter_point[1], segment_value2)) > EPS)
	{
		key1 = get_segment(t_in_x1, inter_point[0], n);
		key2 = get_segment(t_in_x2, inter_point[1], n);
		if (key1 == -100000 || key2 == -100000)
			break;
		segment_value1 = get_segment(t_in_x1, inter_point[0], n);
		segment_value2 = get_segment(t_in_x2, inter_point[1], n);
		find_reversed_Jacobian_value(T1, T2, P_x1, P_x2, P_y1, P_y2, Root_derivative_x1, Root_derivative_x2,
			Root_derivative_y1, Root_derivative_y2, segment_value1, segment_value2, inter_point);

	}
	if (key1 != -100000 && key2 != -100000)
		printf("Intersection found: Parameter value of first spline = %f Parameter value of second spline = %f\n", inter_point[0], inter_point[1]);
	else
		printf("@get_segment: Point is out of spline or SPLINES DO NOT HAVE INTERSECTION POINT. Or you'll want to change inter_point[0] & inter_point[1] in their creation moment.\n\n");


	
	//DISTANCE
	float mdis = 1000;
	float dis = 1000;
	float x01;
	float y01;
	float x02;
	float y02;
	float seg_val1;
	float seg_val2;
	for (float t_value1 = t_in_x1[0]; t_value1 < t_in_x1[n - 1]; t_value1 += 0.009)
	{
		seg_val1 = get_segment_dist(t_in_x1, t_value1, n);
		for (float t_value2 = t_in_x2[0]; t_value2 < t_in_x2[n - 1]; t_value2 += 0.009)
		{
			seg_val2 = get_segment_dist(t_in_x2, t_value2, n);
			x01 = create_p_spline(T1, P_x1, Root_derivative_x1, t_value1, seg_val1);
			y01 = create_p_spline(T1, P_y1, Root_derivative_y1, t_value1, seg_val1);
			x02 = create_p_spline(T2, P_x2, Root_derivative_x2, t_value2, seg_val2);
			y02 = create_p_spline(T2, P_y2, Root_derivative_y2, t_value2, seg_val2);
			dis = pow(x01 - x02, 2) + pow(y01 - y02, 2);
			if (dis < mdis)
				mdis = dis;
		}
	}
	mdis = sqrt(mdis);
	printf("Minimal distance = %f", mdis);
	




	free(X1);
	free(Y1);
	free(X2);
	free(Y2);
	return 0;
}
