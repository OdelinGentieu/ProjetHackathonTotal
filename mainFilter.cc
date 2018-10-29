// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18


#include <iostream>
#include <fstream>
#include <chrono>
#include <mpi.h>
#include "Image.h"
#include "InitMask.h"
#include "ChanVeseSchemes.h"
#include "Util.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }

	config_t c;
	parseFile(argv[1],c);

	MPI_Init(NULL,NULL);
	int me, np;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	Image* image=new Image(); //Pointeur vers l'image

  clock_t t7,t8,t9;
  auto start = chrono::high_resolution_clock::now();
  auto finish = chrono::high_resolution_clock::now();
  if (me == 0)
  {
    cout << "Debut du calcul du filtre pour l'image : " << c.imageName << endl;
    t7=clock();
    //std::cout << t7 << 'azeifgiazyegq '<<(double) t7<< endl;
    //start = chrono::high_resolution_clock::now();
  }

	image->ApplyMedianFilter(c.neighbours,c.coefFilter,c.imageName,me,np);

  if (me == 0)
  {
    cout << "Fin du calcul du filtre pour l'image : " << c.imageName << endl;
    t8=clock();
    t9=t8-t7;
    //std::cout << 'bite' << endl;
    cout << t9*1e-6 << "sec" << endl;
    //finish = chrono::high_resolution_clock::now();
    //double t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
    //cout << "temps pour s'executer :" << t*0.000001 << "secondes"<<endl;
  }

  //image->WriteImage();

	return 0;
}
