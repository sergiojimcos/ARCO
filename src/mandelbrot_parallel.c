/* Conjunto Mandelbrot, M.
 * Sea c un número complejo, n un número natural y {z(n)} una sucesión de números definida por:
 * 		z(0)= 0;
 * 		z(n+1) = z(n)² + c
 * c pertenece a M ssi, para cualquier n, |z(n)| <= 2 

 * Como aproximación, tomaremos los valores de c tales que
 *  para cualquier n < MAX_ITER, |z(n)| <= 2
 
 * Compilación:
 *  gcc mandelbrot_secuencial.c -o mandelbrot_secuencial -lm -fopenmp -O1 -Wno-unused-result
 * 
 * Ejecución: ./mandelbrot_secuencial <archivo de entrada>
 * 
 * El archivo de entrada debe ajustarse al siguiente formato:
 * 	1.ª línea: <número de imágenes que queremos generar>
 * 	Una línea adicional por cada imagen que queramos generar con uno de los siguientes formatos según
 * 		la imagen que queramos sea rectangular o cuadrada, respectivamente:
 * 		1 <menor abscisa> <mayor abscisa> <menor ordenada> <mayor ordenada> <nombre archivo imagen>
 * 		2 <abscisa centro del cuadrado > <ordenada del centro> <lado del cuadrado> <nombre archivo imagen>
 * 		El nombre del archivo imagen no debe tener más de 15 caracteres ni punto.
 
 * Cada imagen será de una zona rectangular (o cuadrada) cuya abscisa menor será xmin,
 *  su abscisa mayor será xmax, su ordenada menor será ymin y su ordenada mayor será ymax.
 * Por tanto, la anchura de la zona es xmax-xmin y la altura es ymax-ymin 
 * En el archivo de imagen (que tendrá el formato ppm) el número de puntos en las direcciones
 *  horizontal y vertical será proporcional a  xmax-xmin e ymax-ymin, respectivamente. */
 
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>

#define R 0 //red colour 
#define B 1 //blue colour 
#define G 2 //green colour 
#define MAX_ITER 5000	// Máximo número de términos calculados en cada sucesión
// Cuanto mayor sea MAX_ITER más calidad de imagen y mayor tiempo de ejecución
#define RK 1000			// RK*RK número aproximado de puntos de la imagen

/* La función mandel_val determina si el número x+i·y
 * pertenece (probablemente) o no (con seguridad) al conjunto de Mandelbrot */
int mandel_val(double x, double y, int max_iter){
	int j = 1;
	double temp_real, cr, p_real = x;
	double ci, p_imag = y;

	cr = x*x; ci = y*y;	// cr+ci es el cuadrado del módulo de x + i·y
	while (((cr+ci)<=4) && (j < max_iter)){
		temp_real = cr - ci + x;
		p_imag = 2 * p_real * p_imag + y;
		p_real = temp_real;
		j++;			// z(j) = p_real + i · p_imag
		cr = p_real * p_real;	ci = p_imag * p_imag;
	}

	if((cr+ci) <= 4)	return -1; else return j;	
};

int main(int argc, char *argv[])	{
	
	int n_imag, imag, t_imag, i, j, value, valuer, valueg, valueb, opt;
	float xmin, xmax, ymin, ymax;
	double x,y,l,t;
	char archivoppm[20];
	FILE* filein;
	int tt = 2;

	// A. Abrir archivo de entrada y leer número de imágenes 	
	if (argc == 2) {
	   filein = fopen(argv[1], "r");
	   if (filein == NULL) {
		  fprintf(stderr, "No puedo abrir %s\n", argv[1]);
		  fprintf(stderr, "uso: %s <archivo de entrada> \n",argv[0]);
		  exit(0);
	   }
	   fscanf(filein,"%d",&n_imag);
	}	else {
			fprintf(stderr, "uso: %s <archivo de entrada> \n",argv[0]);
			exit(0);
		}
	
	// B. Obtener una a una las imágenes y los tiempos respectivos
	for (imag = 0; imag <n_imag; imag ++)	{

		t = omp_get_wtime();
				
		// B1. Leer tipo de imagen (1 = rectangular, 2 = cuadrada)
		fscanf(filein,"%d",&t_imag);
		if (t_imag == 1)	{
			fscanf(filein,"%f %f %f %f %s", &xmin, &xmax, &ymin, &ymax, archivoppm);
		}	else if (t_imag == 2)	{
				fscanf(filein,"%lf %lf %lf %s", &x, &y, &l, archivoppm);
				l = l/2; xmin = x-l; xmax = x+l; ymin = y-l; ymax = y+l;
			}	else {fprintf(stderr, "\n Error en formato del archivo de entrada\n");exit(0);}
			
		// B2. Calcular, en pixels, la anchura, width, y la altura, height, de la imagen
		double k = RK/sqrt((xmax-xmin)*(ymax-ymin));
		int width = k*(xmax-xmin);	// Hay conversión de double a int
		int height = k*(ymax-ymin);	// Íd.
		// width x height = RK² aproximadamente
		// k = width/(xmax-xmin) = height/(ymax-ymin);
		
		// B3. Se abre y configura el archivo de imagen ppm
		strcat(archivoppm,".ppm"); //archivoppm = archivoppm + ".ppm";
		FILE *f_imag = fopen(archivoppm, "wb");
		fprintf(f_imag, "P3\n%d %d\n255\n", width, height);
		/* Formato P3: cada pixel se representa por sus valores de rojo, verde y azul (entre 0 y 255);
		 * width x height: dimensión imagen en pixels */
		 
		// B4. Se calcula el punto (x,y) y se determina si incluir (x + i·y) en el conjunto de Mandelbrot
		// B5. Se determina el color del pixel del número (x + y·i) y se imprime en el archivo
		unsigned char value_colours[height][width][3]; //The number 3 determines number of colours that are: RED, GREEN AND BLUE
		
		#pragma omp parallel for private(x,y,i,j,value)schedule(dynamic, tt)
		for(i = 0; i < height; i++){
			for(j = 0; j < width; j++) {
				x = xmin + j/k; 
				y = ymax - i/k; 
				value = mandel_val(x, y, MAX_ITER);
				
				if(value == -1) {	// se supone que (x+i·y) pertenece a M --> pixel blanco
					value_colours[i][j][R]=255;
					value_colours[i][j][B]=255;
					value_colours[i][j][G]=255;
 
				}else
				// se sabe que (x+i·y) no pertenece a M --> color del pixel según número de iteraciones
				{	value = ((float)value/MAX_ITER) * 16777215;
					// Hacemos que value quede entre 0 y 2²⁴ y extraemos sus componentes rgb
								
					value_colours[i][j][G] = value >> 8;
					value_colours[i][j][B] = value % 256;
					value_colours[i][j][R] = value_colours[i][j][G] >> 8;
					value_colours[i][j][G] = value_colours[i][j][G] % 256;

				}
			}
		}
		#pragma omp barrier
		
		for(i = 0; i < height; i++){
			for(j = 0; j < width; j++) {
				fprintf(f_imag," %d %d %d ",value_colours[i][j][R],value_colours[i][j][G],value_colours[i][j][B]);
			}
		}
		fprintf(f_imag, "\n");	fclose(f_imag);
	
		// B6. Informamos al usuario
		printf("\n En el archivo %s se ha creado una imagen de la zona explorada", archivoppm);
		printf("\n\n Tiempo empleado %0.2f segundos\n\n",omp_get_wtime()-t);
		
	}	// fin de: for (imag = 0; imag <n; imag ++)
	fclose(filein);	
}

/* 
First of all, we have decided in order to use an unsigned char matrix because we suppose that it will
need less memory capacity for our porpuse and it could be more efficient due to we are saving characters that only indicates
the colours of each pixel. In this case, we have used a tridimensional matrix in order to clasiffy each pixel depending 
on the mandelbrot's value, the position and the colour of the pixel.

In this case, we have to take into account several things in order to assert which 
schedule will be the most optimal. In our case we have determined that is a dynamic, bacause the order
of the threads are unknown, which means that if a thread is empty, it will take another segment of code. 
However, the static schedule will spend a small time sorting the threads, and if one of them finishes, it must wait until 
the correspondace one. We have compared both types of schedule with the same tt (2-15) and the lowest value 
is the dynamic one. We made a table for comparing the average time of each iteration.

As well as, the dynamic scheduling type has higher overhead than the static scheduling type because it 
dynamically distributes the iterations during the runtime.

+--------+------+------+------+---------+			
| STATIC | IT1  | IT2  | IT3  | AVERAGE |
+--------+------+------+------+---------+
| TT=2   | 5.54s| 5.52s| 5.58s|    5.54s|
| TT=5   | 5.64s| 5.81s| 5.80s|    5.75s|
| TT=10  | 5.93s| 5.90s| 5.91s|    5.91s|
| TT=15  | 5.95s| 5.96s| 5.85s|    5.92s|		
+--------+------+------+------+---------+ 			

+---------+------+------+------+---------+
| DYNAMIC | IT1  | IT2  | IT3  | AVERAGE |
+---------+------+------+------+---------+
| TT=2    | 5.34s| 5.38s| 5.36s|    5.36s|
| TT=5    | 5.39s| 5.45s| 5.35s|    5.39s|
| TT=10   | 5.37s| 5.44s| 5.73s|    5.51s|
| TT=15   | 5.43s| 6.11s| 5.85s|    5.79s|
+---------+------+------+------+---------+

We proved different values of tt and we establish a range between 2-15 that it is the length of the
number of iterations that a thread done in cycle. The number of threads depends on the number of cores that
you have in your porcessor. This result are executed with a processor of 2 cores. Intel I7-6500u. */