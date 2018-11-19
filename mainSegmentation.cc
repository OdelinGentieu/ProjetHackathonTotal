// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18


#include <iostream>
#include <fstream>
#include <Dense>
#include <chrono>
#include "Image.h"
#include "InitMask.h"
#include "ChanVeseSchemes.h"
#include "Util.h"
#include <cmath>
#include <mpi.h>
using namespace std;


int main(int argc, char** argv)
{
  // Lecture du fichier de données
  if (argc < 2)
  {
    cout << "Il faut indiquer le nom du fichier de données." << endl;
    abort();
  }
  config_t c;
  parseFile(argv[1],c);

  // Construction des différents noms de fichiers de sortie
  std::string extension(c.imageName);
  extension = extension.substr(extension.find_last_of(".") + 1);
  std::string namewithoutextension(c.imageName);
  namewithoutextension.erase(namewithoutextension.rfind('.'));
  c.imageName = namewithoutextension+"_filtered.tiff";
  std::string imagemaskdistance((namewithoutextension+"_filtered_distance_mask.vtk").c_str());

  // Créer un dossier
  system("mkdir -p ./Results");
  // et supprimer les anciens fichiers
  system("rm -f ./Results/*.vtk");

  // Lecture de l'image
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Lecture de l'image " << c.imageName << std::endl;
  Image* image = new Image();
  image->ReadImage(c.imageName);
  saveVTKFile(image->GetImage(), "Results/image.vtk");
  std::cout << "-------------------------------------------------" << std::endl;

  // Initialisation ou lecture du masque redistancié
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Construction ou lecture du masque redistancié" << std::endl;
  InitMask* initMask = new InitMask();
  int rows(image->GetImage().rows()), cols(image->GetImage().cols());
  initMask->BuildMaskAndRedistancing(rows, cols, imagemaskdistance);
  std::cout << "-------------------------------------------------" << std::endl;

  // Chan Vese method
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Initialisation de la méthode de Chan Vese" << std::endl;
  ChanVeseSchemes* chanVese=new ChanVeseSchemes(image);
  field phi(initMask->GetRedistMask());
  saveVTKFile(phi, "Results/sol_0.vtk");
  std::cout << "-------------------------------------------------" << std::endl;

  std::cout << "Iteration -- 1" << std::endl;





  /////////////////// !! BOUCLE: SCHEMA + COMMUNICATION + CL !! ////////////////////
  //Initialisation des termes de corrections commun a tous les proc/IdBlock
  double C1,C2;
  (C1,C2)= chanVese->Correction(phi);


  // Schema
  //#pragma acc parallel loop reduction(max:diff)
  field newphi = chanVese->ExplicitScheme(phi,c.dt,c.mu,c.nu,c.l1,c.l2,C1,C2);

  // Condition Neumann au bords + nouvelle correction TODO buffer for boundary
  //#pragma acc update self(newphi)

  newphi.block(0,0,1,newphi.cols()) = newphi.block(1,0,1,newphi.cols()); // CL haut
  newphi.block(newphi.rows(),0,1,newphi.cols()) = newphi.block(newphi.rows() -1,0,1,newphi.cols()); // CL Bas
  newphi.block(0,0,newphi.rows(),1) = newphi.block(0,1,newphi.rows(),1); // CL gauche
  newphi.block(0,newphi.cols(),newphi.rows(),1) = newphi.block(0,newphi.cols()-1,newphi.rows(),1); // CL droite

  //#pragma acc update device(newphi)
  ////////////////////////////// !! FIN BOUCLE !! //////////////////////////////


  double diff = (((newphi>=0).cast<double>()-0.5)*2. - ((phi>=0).cast<double>()-0.5)*2.).matrix().norm()
                /(phi.rows()*phi.cols());
  phi = newphi;

  int i(2);
  //////////////////// !! PREPARATION GPU/CPU !! ////////////////////
  #pragma acc data copy(phi) create(newphi)
  //////////////////// !! FIN PREPARATION !! ////////////////////////

  #pragma acc parallel loop reduction(max:diff)
  while ( (diff > 5e-6) && (i < 100) )
  {
    // Schema
    newphi = chanVese->ExplicitScheme(phi,c.dt,c.mu,c.nu,c.l1,c.l2,C1,C2);

    // Condition Neumann au bords + nouvelle correction TODO buffer for boundary
    #pragma acc update self(newphi)

    newphi.block(0,0,1,newphi.cols()) = newphi.block(1,0,1,newphi.cols()); // CL haut
    newphi.block(newphi.rows(),0,1,newphi.cols()) = newphi.block(newphi.rows() -1,0,1,newphi.cols()); // CL Bas
    newphi.block(0,0,newphi.rows(),1) = newphi.block(0,1,newphi.rows(),1); // CL gauche
    newphi.block(0,newphi.cols(),newphi.rows(),1) = newphi.block(0,newphi.cols()-1,newphi.rows(),1); // CL droite

    double diff = (((newphi>=0).cast<double>()-0.5)*2. - ((phi>=0).cast<double>()-0.5)*2.).matrix().norm()
                  /(phi.rows()*phi.cols());

    if (i%10 == 0)
    {
      newphi = ((newphi>=0).cast<double>()-0.5)*2;
      LevelSet lv(newphi);
      cout << "Redistanciation... " << endl;
      lv.redistancing(10);
      saveVTKFile(newphi, ("Results/sol_" + to_string(i)+ ".vtk").c_str());
      std::cout << "Evolution of phi : " << diff << std::endl;
    }
    phi = newphi;
    i++;

    #pragma acc update device(phi)
  }


  newphi = ((newphi>=0).cast<double>()-0.5)*2;
  saveVTKFile(newphi, "Results/lastsol.vtk");
  image->WriteImage(chanVese->AbsGradPhi(newphi), (namewithoutextension + "_filtered_with_contour." + extension).c_str());

  cout << "Fin de la segmentation pour l'image." << endl;

  return 0;
}
