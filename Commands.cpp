Configure Baci for debug or release Version:
./do-configure 

Compile after ./do-configure
make -j8

Run input file beam3ebtor_endmoment_endforce.dat with output name test2 and write in test3.txt:
./baci-release Input/beam3ebtor_endmoment_endforce.dat o/test2 | tee test3.txt

//Postprocess with drt-ensight-filter 
//Run input file with parallel processors: (-np=number of processors)
mpirun -np 4 ./post_drt_ensight --file=o/xxx --output=o/yyy

//Posprocess the results in the file test2:  (Postprocess in drt-ensight)
./post_drt_ensight --file=o/test2

//Postprocess with g-mesh
/*Go to Gmesh folder:*/ cd GmshOutput/
/*Run movie:*/ python ./make_movie_gmsh.py network000 output_config_network.geo 
/*Replay with VLC*/ vlc ~/tmp/x.mpg


//Include test file 
modify TestingFramework.cmake

// Configure Baci debug
 ../do-configure --debug --useACML

 //Print ids in C++
 (this->Id()==)
 
 //Add new files in baci 
 Open Cmakefilelists.txt
 create appropriate folder and add only *.cpp files in the folder.
 
 //Add new test files in baci 
 Modify TestingFramework.cmake file
 
 //Get line number and file name in C++
  std::cout<<__FILE__<<__LINE__ <<std::endl;

 // Teuchos::null hands over nothing
 // []denotes vector and () denotes matrices
 
 //Debug baci with kdbg
 mpirun -np 2 xterm -e kdbg ./baci-debug
 
 //run ctest for beam elements
 ctest -R beam
 //run ctest for statmech manager
 ctest -R statmech
 
 //svn commands
 svn stat % svn statistics
 
// Revert file solver_aztecsolver.cpp to it's previous versions
 svn revert src/solver/solver_aztecsolver.cpp
 
// Find difference in files
 svn diff src/solver/solver_aztecsolver.cpp

// Give the version number 
svn info

// copy the differences in a file mychanges.patch
svn diff>mychanges.patch

// in new baci copy the patch files
cp mychanges.patch ../baci_new/

// Apply patch in new baci
patch -p0<mychanges.patch

Note: Add to version control

// Better to use Linalg matrices since storage space is restricted.

//Structural dynamics properties
const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

//ssh -X

// Sparse Operator
  Teuchos::RCP<LINALG::SparseMatrixBase> mySparseMatrix =   Teuchos::rcp_dynamic_cast<LINALG::SparseMatrixBase>(stiff);
  std::cout<<"Mymatrix"<<*(mySparseMatrix->EpetraMatrix())<<std::endl;
  
// Print Matrix in matlab format
LINALG::PrintMatrixInMatlabFormat("matrix.m",*(mySparseMatrix->EpetraMatrix()),true);
  
// Live Mid
It's keyword for denoting orthogonal pressure. 

//Running job on cluster
qsub jobscript.sh

//compare folder through ssh in kaiser both named "src"
 diff  <(ssh mukherjee@kaiser ls -R /home/mukherjee/workspace/baci/src) <(ls -R home/mukherjee/workspace/baci/src)

//GitHub
gitk. gitgui

//random numbers
in oder to generate fix set of random numbers for every simulation run comment srand() function.
That is to say don't seed the random number generator. Line number 2317 in input generator

//fixed random numbers
   /* This part of code ensures that regardless of element type, same set of random
    * numbers are used. And that for every gauss points the random numbers remain
    * same */
   // create seed
   unsigned long seed =12411;
   // produces general pseudo random number with MT
   boost :: mt19937 rng(seed);
   // KT=0.00404531
   double std_dv=pow(2.0 * 0.00404531 / params.get<double>("delta time"),0.5);
   // first mean, second standard deviation
   boost::normal_distribution<> nd(0, std_dv);
   boost::variate_generator< boost::mt19937, boost::normal_distribution<> > dist(rng, nd);

   for (int i=0; i<randomnumbers->MyLength(); i++)
     for (int j=0; j<randomnumbers->NumVectors(); j++)
     {
       (*randomnumbers)[j][i] = dist();
 //      std::cout<<"(*randomnumbers)["<<j<<"]["<<i<<"]="<<(*randomnumbers)[j][i]<<std::endl;
     }
   rng.seed(seed);
   
        // Use same set of gauss numbers for all dofs
     for(int i=0; i<3; i++)
     {
       (*randomnumbers)[i+4][LID()]= (*randomnumbers)[i][LID()];
       (*randomnumbers)[i+8][LID()]= (*randomnumbers)[i][LID()];
     }
     //    std::cout<<"randomnumbers after="<<(*randomnumbers)<<std::endl;
 
// Commmands to run rheological simulations
1. Run a normal simulation with all the forces present.
2. Create a new file with name frequencies.txt with the values inside.
3. python CreateViscoElJobs.py ViscEl (ViscEl being the folder for rheological simulations)
4. python CopyDataFromTo.py ViscEl (copies output data from the previously ran simulation to the other folder)
5. python SetRestartJobs.py ViscEl 
6. python RunJobs.py ViscEl

//Commmands to obtain plots of rheological simulations
1. Create a directory structure same as in Cluster
2. Copy MoveData.py to the present directory
3. python ./MoveData.py ViscEl (ViscEl being the directory on Cluster)
    3.1) Give the file name TRUNK 1. file: StatMechOutput/ViscoEl
4. In case the script fails to copy, check the paths.txt file. The line of paths.txt should end with a backslash "/"

// Create python history file (.python_history) for allowing recall of commands with arrow keys. 

//dyanamic_cast comment
take an element and as it if it's like a particular elementtype. If't it is we can acees the element properties. 
If not it will return a null vector.

//Print containers of variables for a differnt variabletype
1. Check if the variabletype has a print method defined e.g. var.print().
2. and then var.print(std::cout)

//While creating new methods it is always important to use reference (&) or pointer (*)
Otherwise, the member function will copy the object and that might create problems in future. 

//Get properties of other elements using dynamic casts
//ask the truss element about the first element the first node is connected to
    DRT::Element* Element1 = Nodes()[j]->Elements()[0];
//Check via dynamic cast, if it's a beam3eb elements
DRT::ELEMENTS::Beam3eb* BeamElement = dynamic_cast<DRT::ELEMENTS::Beam3eb*>(Element1);
trefNodeAux=BeamElement->Tref()[j];

//To find processor id of an element 
this->owner

//display node ID and element id 

//debugging in mpirun
./baci-debug input.dat xxx --interactive

new shell:
gdb attach process_id

old shell:
press any key

new shell:
breakpoint, cont, ... 

//Printout Gmshfile content
 std::cout<<"gmshfilecontent"<<__FILE__<<__LINE__<<"="<<gmshfilecontent.str()<<std::endl;
 
# Epetra_Vector or Epetra_MultiVector and RCP
 Teuchos::RCP<Epetra_MultiVector> bspottriadscol= Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));
 // what this statement says is that there is a pointer named bspottriadscol. This pointer points to an object
 // of Epetra_MultiVector type. And this Epetra_MultiVector inherits its properties from a pointer to a map
 // named bspotcolmap_ and contains 4 vectors of same size with all the components of the vector being intialised 
 // to zero.
 
#Print out GMSH output
 // In method int STR::TimIntImpl::NewtonFull() in Strtimint_impl.cpp after 
 // call monitor
 " if (conman_->HaveMonitor())
  {
    conman_->ComputeMonitorValues(disn_);
  }" // add the follwing code
  //********************Begin: GMSH Output*******************************
    //This output works only with one proc!!!!
    // STEP 1: OUTPUT OF TIME STEP INDEX
    std::ostringstream filename;
    filename << "/home/meier/workspace/o/gmsh_output/";

    //filename << "o/gmsh_output/";
    if (step_<1000000)
      filename << "beams_t" << std::setw(6) << std::setfill('0') << step_+1;
    else /*(timestep>=1000000)*/
      dserror("ERROR: Gmsh output implemented for max 999.999 time steps");

    // finish filename
    filename<<".pos";

    // do output to file in c-style
    FILE* fp = NULL;

    // open file to write output data into
    fp = fopen(filename.str().c_str(), "w");

    std::stringstream gmshfileheader;
    gmshfileheader <<"General.RotationCenterGravity=0;\n";
    gmshfileheader <<"General.RotationCenterX=11;\n";
    gmshfileheader <<"General.RotationCenterY=11;\n";
    gmshfileheader <<"General.RotationCenterZ=11;\n";

    //write content into file and close it
    fprintf(fp, gmshfileheader.str().c_str());
    fclose(fp);

    //fp = fopen(filename.str().c_str(), "w");
    fp = fopen(filename.str().c_str(), "a");

    std::stringstream gmshfilecontent;
    gmshfilecontent << "View \" Step T" << step_;

    // finish step information
    gmshfilecontent << " \" {" << std::endl;

    //loop over fully overlapping column element map of proc 0
    for (int i=0;i<discret_->ElementColMap()->NumMyElements();++i)
    {

      // get pointer onto current beam element
      DRT::Element* element = discret_->lColElement(i);

      // prepare storage for nodal coordinates
      int nnodes = element->NumNode();
      LINALG::SerialDenseMatrix coord(3,2);

      // compute current nodal positions
      for (int id=0;id<3;++id)
      {
        for (int jd=0;jd<2;++jd)
        {
          double referenceposition = ((element->Nodes())[jd])->X()[id];
          std::vector<int> dofnode = discret_->Dof((element->Nodes())[jd]);
          double displacement = (*disn_)[discret_->DofColMap()->LID(dofnode[id])];
          coord(id,jd) =  referenceposition + displacement;
        }
      }

      gmshfilecontent << "SL("<< std::scientific;
      gmshfilecontent << coord(0,0) << "," << coord(1,0) << "," << coord(2,0) << ",";
      gmshfilecontent << coord(0,1) << "," << coord(1,1) << "," << coord(2,1);
      gmshfilecontent << "){" << std::scientific;
      gmshfilecontent << 1.0 << "," << 1.0 << "};" << std::endl;

    }//loop over elements

    // finish data section of this view by closing curley brackets (somehow needed to get color)
    gmshfilecontent <<"SP(0.0,0.0,0.0){0.0,0.0};"<<std::endl;
    gmshfilecontent << "};" << std::endl;

    // write content into file and close
    fprintf(fp,gmshfilecontent.str().c_str());
    fclose(fp);
    //********************End: GMSH Output:******************************* 




