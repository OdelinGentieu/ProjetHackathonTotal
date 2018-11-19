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

  //////////////////// !! PREPARATION GPU/CPU !! ////////////////////////
  MPI_Status status;
  MPI_Init(&argc,&argv);
  int me, Np;
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  int chargex = chanVese->charge(me,Np,phi.rows()); //TODO
  std::vector<int> charge_min(Np);
  std::vector<int> Point_min(Np);


  // Communication Proc me/ Idblock.X --> Proc 0 [-------,chargex,------] a  [2*me:2*me +1] place
  if (me==0){
    charge_min[0]=chargex;
    for (size_t i = 1; i < Np; i++) {
      MPI_Recv(&charge_min[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
    };
  }
  else{
    MPI_Send(&chargex, 1, MPI_INT, 0, me, MPI_COMM_WORLD);
  };

  if (me==0){
    Point_min[0]=0;
    for (size_t i = 1; i < Np; i++) {
      Point_min[i]=Point_min[i-1]+ charge_min[i-1];
    }
  }

  // Communication Proc 0 -> Proc me/IdBlock.X [------,x_min(me),------] --> FACULTATIF ??
  double xmin;
  if (me==0){
    xmin=0;
    for (size_t i = 1; i < Np; i++) {
      MPI_Send(&Point_min[i], 1, MPI_INT, i, i+1000, MPI_COMM_WORLD);
    };
  }
  else{
    MPI_Recv(&xmin, 1, MPI_INT, 0, me+1000, MPI_COMM_WORLD, &status);
  };
  //////////////////// !! FIN PREPARATION !! ////////////////////////



  /////////////////// !! BOUCLE: SCHEMA + COMMUNICATION + CL !! ////////////////////
  //Initialisation des termes de corrections commun a tous les proc/IdBlock
  double C1,C2;
  if (me==0){
    (C1,C2)= chanVese->Correction(phi);
    for (size_t i = 1; i < Np; i++) {
      MPI_Send(&(C1,C2), 2, MPI_FLOAT, i, i+2000, MPI_COMM_WORLD);
    };
  }
  else{
    MPI_Recv(&(C1,C2), 2, MPI_FLOAT, 0, me+2000, MPI_COMM_WORLD, &status);
  }

  //Initialisation phi_pti = petit "block" de phi calcul par chaque Block/proc
  field phi_pti(chargex+2, phi.cols());
  phi_pti = phi.block(xmin,0,chargex+2,phi.cols());

  // Schema
  std::cout << "Explicit scheme" << std::endl;
  phi_pti = chanVese->ExplicitScheme(phi_pti,c.dt,c.mu,c.nu,c.l1,c.l2,C1,C2);

  // Recréation du phi au temps T+1 + Communication proc me/IdBlock.X --> Proc 0
  //newphi.block(xmin,0,chargex+2,phi.cols()) =phi_pti;
  if (me==0){
    field newphi(phi.rows(),phi.cols());
    float *Buf_phi;
    for (size_t i = 1; i < Np; i++) {
      Buf_phi = (float*) malloc((charge_min[i]+2)*phi.cols()*sizeof(float));
      MPI_Recv(Buf_phi  , (charge_min[i]+2)*phi.cols(), MPI_FLOAT, i, i+3000, MPI_COMM_WORLD, &status);
      newphi.block(Point_min[i],0,charge_min[i]+2,phi.cols()).data()=*Buf_phi;
      free(Buf_phi);
    };
  }
  else{
    MPI_Send(&phi_pti,(chargex+2)*phi.cols() , MPI_FLOAT, 0, me+3000, MPI_COMM_WORLD);
  }

  // Condition Neumann au bords + nouvelle correction
  if (me==0){
    newphi.block(0,0,1,newphi.cols()) = newphi.block(1,0,1,newphi.cols()); // CL haut
    newphi.block(newphi.rows(),0,1,newphi.cols()) = newphi.block(newphi.rows() -1,0,1,newphi.cols()); // CL Bas
    newphi.block(0,0,newphi.rows(),1) = newphi.block(0,1,newphi.rows(),1); // CL gauche
    newphi.block(0,newphi.cols(),newphi.rows(),1) = newphi.block(0,newphi.cols()-1,newphi.rows(),1); // CL droite
    (C1,C2)= chanVese->Correction(newphi);
    for (size_t i = 1; i < Np; i++) {
      MPI_Send(&(C1,C2), 2, MPI_FLOAT, i, i+4000, MPI_COMM_WORLD);
    };
  }
  else{
    MPI_Recv(&(C1,C2), 2, MPI_FLOAT, 0, me+4000, MPI_COMM_WORLD, &status);
  }
  ////////////////////////////// !! FIN BOUCLE !! //////////////////////////////


  double diff = (((newphi>=0).cast<double>()-0.5)*2. - ((phi>=0).cast<double>()-0.5)*2.).matrix().norm()
                /(phi.rows()*phi.cols());

  phi = newphi;
  int i(2);
  while ( (diff > 5e-6) && (i < 100) )
  {
    if (i%10 == 0) { std::cout << "Iteration -- " << i << std::endl;}

    /////////////////// !! BOUCLE: SCHEMA + COMMUNICATION + CL !! ////////////////////
    //Initialisation phi_pti = petit "block" de phi calcul par chaque Block/proc
    phi_pti = phi.block(xmin,0,chargex+2,phi.cols());

    // Schema
    phi_pti = chanVese->ExplicitScheme(phi_pti,c.dt,c.mu,c.nu,c.l1,c.l2,C1,C2);

    // Recréation du phi au temps T+1 + Communication proc me/IdBlock.X --> Proc 0
    newphi.block(xmin,0,chargex+2,phi.cols()) = phi_pti;

    // Condition Neumann au bords + nouvelle correction
    newphi.block(0,0,1,newphi.cols()) = newphi.block(1,0,1,newphi.cols()); // CL haut
    newphi.block(newphi.rows(),0,1,newphi.cols()) = newphi.block(newphi.rows() -1,0,1,newphi.cols()); // CL Bas
    newphi.block(0,0,newphi.rows(),1) = newphi.block(0,1,newphi.rows(),1); // CL gauche
    newphi.block(0,newphi.cols(),newphi.rows(),1) = newphi.block(0,newphi.cols()-1,newphi.rows(),1); // CL droite
    (C1,C2)= chanVese->Correction(newphi);
    ////////////////////////////// !! FIN BOUCLE !! //////////////////////////////


    diff = (((newphi>=0).cast<double>()-0.5)*2. - ((phi>=0).cast<double>()-0.5)*2.).matrix().norm()
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
  }

  int MPI_Finalize();

  newphi = ((newphi>=0).cast<double>()-0.5)*2;
  saveVTKFile(newphi, "Results/lastsol.vtk");
  image->WriteImage(chanVese->AbsGradPhi(newphi), (namewithoutextension + "_filtered_with_contour." + extension).c_str());

  cout << "Fin de la segmentation pour l'image." << endl;

  return 0;
}
