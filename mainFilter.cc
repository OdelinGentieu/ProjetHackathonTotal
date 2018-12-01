// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18


#include <iostream>
#include <fstream>
#include <chrono>
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

  int me=0, np=1;

  Image* image=new Image(); //Pointeur vers l'image

  cout << "Debut du calcul du filtre pour l'image : " << c.imageName << endl;

  image->ApplyMedianFilter(c.neighbours,c.coefFilter,c.imageName,me,np);

  
  cout << "Fin du calcul du filtre pour l'image : " << c.imageName << endl;

  //image->WriteImage();

	return 0;
}
