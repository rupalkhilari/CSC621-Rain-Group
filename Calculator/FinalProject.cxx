#include <vtkSmartPointer.h>
#include <vtkCurvatures.h>
#include <vtkPointSource.h>
#include <vtkPoints.h>
#include <vtkKochanekSpline.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkImageViewer2.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMetaImageReader.h>
#include <vtkObjectFactory.h>
#include <vtkResliceImageViewer.h>
#include <vtkTransform.h>
#include <vtkAxesActor.h>


// added this For the spine
#include <vtkOutlineFilter.h>
#include <vtkCamera.h>
#include <vtkProperty.h>
#include <vtkMarchingCubes.h>
#include <vtkImageFlip.h>

// added this for the points
#include <vtkVertexGlyphFilter.h>
#include <vtkPointData.h>

// testing this
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkSmartPointer.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkImageCast.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
// finished testing this

#include "lib/scCalc.hh"
#include "lib/Registration.hh"
#include "lib/RegionGrowingNoThreshold.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

#include "itkFlipImageFilter.h"

using namespace std;



//seed coords
double seedX, seedY, seedZ;
bool useDatabase = false;
bool fitOnly = false;
bool useAnnotations = false;
bool useComponents = false;
bool newFile1 = true;
bool newFile2 = true;
//Helper funcion, ignore this
vtkSmartPointer<vtkActor> makeLine(double data[][3], unsigned length, double color[3])
{

  vtkSmartPointer<vtkPoints> points =
  vtkSmartPointer<vtkPoints>::New();

  for (int i = 0; i < length; ++i)
  {
    points->InsertPoint(i, data[i]);
  }

  vtkSmartPointer<vtkKochanekSpline> xSpline =
  vtkSmartPointer<vtkKochanekSpline>::New();
  vtkSmartPointer<vtkKochanekSpline> ySpline =
  vtkSmartPointer<vtkKochanekSpline>::New();
  vtkSmartPointer<vtkKochanekSpline> zSpline =
  vtkSmartPointer<vtkKochanekSpline>::New();

  vtkSmartPointer<vtkParametricSpline> spline =
  vtkSmartPointer<vtkParametricSpline>::New();
  spline->SetXSpline(xSpline);
  spline->SetYSpline(ySpline);
  spline->SetZSpline(zSpline);
  spline->SetPoints(points);

  vtkSmartPointer<vtkParametricFunctionSource> functionSource =
  vtkSmartPointer<vtkParametricFunctionSource>::New();
  functionSource->SetParametricFunction(spline);
  functionSource->Update();

  // Setup actor and mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper =
  vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(functionSource->GetOutputPort());

  vtkSmartPointer<vtkActor> actor =
  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(color[0], color[1], color[2]); //(R,G,B)

  return actor;
}

// added this - for the spine
vtkSmartPointer<vtkActor> makeSpine(char filename[]) {
  // read the image
  vtkSmartPointer<vtkMetaImageReader> reader =
    vtkSmartPointer<vtkMetaImageReader>::New();
  reader->SetFileName (filename);

  vtkSmartPointer<vtkMarchingCubes> skinExtractor =
    vtkSmartPointer<vtkMarchingCubes>::New();
  skinExtractor->SetInputConnection(reader->GetOutputPort());
  skinExtractor->SetValue(0, 800);

  vtkSmartPointer<vtkPolyDataMapper> skinMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  skinMapper->SetInputConnection(skinExtractor->GetOutputPort());
  skinMapper->ScalarVisibilityOff();

  vtkSmartPointer<vtkActor> skin =
    vtkSmartPointer<vtkActor>::New();
  skin->SetMapper(skinMapper);
  skin->GetProperty()->SetDiffuseColor(1, 1, .9412);
  skin->GetProperty()->SetOpacity(.5);

  return skin;
}


// added this:
vtkSmartPointer<vtkActor> drawPoints(std::vector<std::vector<int>> centroids, double* spacing, std::vector<std::string> vertebraePoints) {

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  unsigned length = centroids.size();
  for (unsigned i = 0; i < length; ++i)
  {
      points->InsertNextPoint(double(centroids[i][0] * spacing[0]), double(centroids[i][1] * spacing[1]), double(centroids[i][2] * spacing[2]));
  }

  /*points->InsertNextPoint (0.0, 0.0, 0.0);
  points->InsertNextPoint (1.0, 0.0, 0.0);
  points->InsertNextPoint (0.0, 1.0, 0.0);*/
 
  vtkSmartPointer<vtkPolyData> pointsPolydata =
    vtkSmartPointer<vtkPolyData>::New();
 
  pointsPolydata->SetPoints(points);
 
  vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
#else
  vertexFilter->SetInputData(pointsPolydata);
#endif
  vertexFilter->Update();
 
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->ShallowCopy(vertexFilter->GetOutput());
 
  // Setup colors
  unsigned char red[3] = {255, 0, 0};
  unsigned char green[3] = {0, 255, 0};
  unsigned char blue[3] = {0, 0, 255};
  unsigned char cyan[3] = {0, 255, 255};
  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName ("Colors");

  if (vertebraePoints.size() > 0) {
    for (unsigned i = 0; i < length; ++i)
    {
        if (vertebraePoints[i].at(0) == 'C') {
          colors->InsertNextTupleValue(cyan);
        }
        else if (vertebraePoints[i].at(0) == 'T') {
          colors->InsertNextTupleValue(red);
        }
        else if (vertebraePoints[i].at(0) == 'L') {
          colors->InsertNextTupleValue(green);
        }
        else if (vertebraePoints[i].at(0) == 'S'){
          colors->InsertNextTupleValue(blue);
        }
        else {
          colors->InsertNextTupleValue(red); // Default to red
        }
    }  
  }
  else {
    for (unsigned i = 0; i < length; ++i)
    {
        switch (i%3) {
          case 0:
            colors->InsertNextTupleValue(red);
            break;
          case 1:
            colors->InsertNextTupleValue(green);
            break;
          case 2:
            colors->InsertNextTupleValue(blue);
            break;
          default:
            colors->InsertNextTupleValue(red);
        }
    }
  }
  polydata->GetPointData()->SetScalars(colors);
 
  // Visualization
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInputConnection(polydata->GetProducerPort());
#else
  mapper->SetInputData(polydata);
#endif
 
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetPointSize(5);
 
  return actor;
}

// added this:
vtkSmartPointer<vtkActor> drawPointsDouble(double centroids[][3], double* spacing, int length, std::vector<std::string> vertebraePoints) {

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  for (unsigned i = 0; i < length; ++i)
  {
      points->InsertNextPoint(centroids[i][0], centroids[i][1], centroids[i][2]);
  }

  /*points->InsertNextPoint (0.0, 0.0, 0.0);
  points->InsertNextPoint (1.0, 0.0, 0.0);
  points->InsertNextPoint (0.0, 1.0, 0.0);*/
 
  vtkSmartPointer<vtkPolyData> pointsPolydata =
    vtkSmartPointer<vtkPolyData>::New();
 
  pointsPolydata->SetPoints(points);
 
  vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
#else
  vertexFilter->SetInputData(pointsPolydata);
#endif
  vertexFilter->Update();
 
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->ShallowCopy(vertexFilter->GetOutput());
 
  // Setup colors
  unsigned char red[3] = {255, 0, 0};
  unsigned char green[3] = {0, 255, 0};
  unsigned char blue[3] = {0, 0, 255};
  unsigned char cyan[3] = {0, 255, 255};

  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName ("Colors");
  if (vertebraePoints.size() > 0) {
    for (unsigned i = 0; i < length; ++i)
    {
        if (vertebraePoints[i].at(0) == 'C') {
          colors->InsertNextTupleValue(cyan);
        }
        else if (vertebraePoints[i].at(0) == 'T') {
          colors->InsertNextTupleValue(red);
        }
        else if (vertebraePoints[i].at(0) == 'L') {
          colors->InsertNextTupleValue(green);
        }
        else if (vertebraePoints[i].at(0) == 'S'){
          colors->InsertNextTupleValue(blue);
        }
        else {
          colors->InsertNextTupleValue(red); // Default to red
        }
    }    
  }
  else {
    for (unsigned i = 0; i < length; ++i)
    {
        switch (i%3) {
          case 0:
            colors->InsertNextTupleValue(red);
            break;
          case 1:
            colors->InsertNextTupleValue(green);
            break;
          case 2:
            colors->InsertNextTupleValue(blue);
            break;
          default:
            colors->InsertNextTupleValue(red);
        }
    }
  }
  polydata->GetPointData()->SetScalars(colors);
 
  // Visualization
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInputConnection(polydata->GetProducerPort());
#else
  mapper->SetInputData(polydata);
#endif
 
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetPointSize(5);
 
  return actor;
}


//mouse event handler
class MouseInteractorStyle3 : public vtkInteractorStyleTrackballCamera
{
public:
  static MouseInteractorStyle3* New();

  virtual void OnLeftButtonDown() 
  {
    seedX = this->Interactor->GetEventPosition()[0];
    seedY = this->Interactor->GetEventPosition()[1];
    std::cout << "Please close the window to continue...\n";
  }

};


void drawImageWithPCurve(char filename[], double centroids[][3], double* spacing, int length, int type, double color[3], std::vector<std::string> vertebraePoints) {
    //read input mhd file
  vtkSmartPointer<vtkMetaImageReader> reader =
  vtkSmartPointer<vtkMetaImageReader>::New();
  reader->SetFileName(filename);
  reader->Update();

  vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
  //reslice->SetOutputExtent(0, 9, 0, 100, 0, 0);
  reslice->SetOutputSpacing(1,1,spacing[2]);
  reslice->SetInputConnection(reader->GetOutputPort());
  reslice->Update();

  //display
  vtkSmartPointer<vtkResliceImageViewer> imageViewer =
  vtkSmartPointer<vtkResliceImageViewer>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor2 = 
  vtkSmartPointer<vtkRenderWindowInteractor>::New();

  imageViewer->SetResliceModeToAxisAligned();
  imageViewer->SetInputConnection(reslice->GetOutputPort());
  imageViewer->SetupInteractor(renderWindowInteractor2);
  imageViewer->SetColorLevel(500);
  imageViewer->SetColorWindow(2000);

  //set z coord always the most center slice 
  seedZ = (imageViewer->GetSliceMax())/2;

  imageViewer->SetSize(reader->GetWidth(), reader->GetHeight());
  imageViewer->SetSlice(seedZ);
  //imageViewer->SetSliceOrientationToXY();
  if (type == 2) {
    imageViewer->SetSliceOrientationToYZ();
    for (unsigned i = 0; i < length; i++) {
      centroids[i][0] = 400;
    }
  } else if (type == 1) {
    imageViewer->SetSliceOrientationToXZ();
  }
  imageViewer->GetRenderer()->AddActor(drawPointsDouble(centroids, spacing, length, vertebraePoints));
  imageViewer->GetRenderer()->AddActor(makeLine(centroids, length, color));

    // Setup the text and add it to the renderer
  /*vtkSmartPointer<vtkTextActor> textActor = 
  vtkSmartPointer<vtkTextActor>::New();
  textActor->SetInput ( "Hello world" );
  textActor->SetPosition2 ( 10, 40 );
  textActor->GetTextProperty()->SetFontSize ( 24 );
  textActor->GetTextProperty()->SetColor ( 1.0, 0.0, 0.0 );
  imageViewer->GetRenderer()->AddActor2D ( textActor );


  for (int i = 0; i < vertebraePoints.size(); i++) {
    // Setup the text and add it to the renderer
    vtkSmartPointer<vtkTextActor> textActor = 
    vtkSmartPointer<vtkTextActor>::New();
    textActor->SetInput ("T");
    textActor->SetPosition2 ( centroids[i][0], centroids[i][1] );
    textActor->GetTextProperty()->SetFontSize ( 12 );
    textActor->GetTextProperty()->SetColor ( 1.0, 0.0, 0.0 );
    imageViewer->GetRenderer()->AddActor2D ( textActor );    
  }*/
  imageViewer->GetRenderer()->ResetCamera();
  imageViewer->Render();
    seedZ = imageViewer->GetSlice();

  vtkSmartPointer<MouseInteractorStyle3> style =
  vtkSmartPointer<MouseInteractorStyle3>::New();
  
  renderWindowInteractor2->SetInteractorStyle(style);
  //renderWindowInteractor2->UpdateSize(100,100);
  renderWindowInteractor2->Start();
}

void drawImageWithCurve(char filename[], std::vector<std::vector<int>> centroids, double* spacing, std::vector<std::string> vertebraePoints) {
    //read input mhd file
  vtkSmartPointer<vtkMetaImageReader> reader =
  vtkSmartPointer<vtkMetaImageReader>::New();
  reader->SetFileName(filename);
  reader->Update();

  vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
  //reslice->SetOutputExtent(0, 9, 0, 100, 0, 0);
  reslice->SetOutputSpacing(1,1,spacing[2]);
  reslice->SetInputConnection(reader->GetOutputPort());
  reslice->Update();

  //display
  vtkSmartPointer<vtkResliceImageViewer> imageViewer =
  vtkSmartPointer<vtkResliceImageViewer>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor2 = 
  vtkSmartPointer<vtkRenderWindowInteractor>::New();

  imageViewer->SetResliceModeToAxisAligned();
  imageViewer->SetInputConnection(reslice->GetOutputPort());
  imageViewer->SetupInteractor(renderWindowInteractor2);
  imageViewer->SetColorLevel(500);
  imageViewer->SetColorWindow(2000);

  //set z coord always the most center slice 
  seedZ = (imageViewer->GetSliceMax())/2;

  imageViewer->SetSize(reader->GetWidth(), reader->GetHeight());
  imageViewer->SetSlice(seedZ);
  //imageViewer->SetSliceOrientationToXY();
  imageViewer->SetSliceOrientationToYZ();
  imageViewer->GetRenderer()->AddActor(drawPoints(centroids, spacing, vertebraePoints));
  imageViewer->GetRenderer()->ResetCamera();
  imageViewer->Render();
    seedZ = imageViewer->GetSlice();

  vtkSmartPointer<MouseInteractorStyle3> style =
  vtkSmartPointer<MouseInteractorStyle3>::New();
  
  renderWindowInteractor2->SetInteractorStyle(style);
  //renderWindowInteractor2->UpdateSize(100,100);
  renderWindowInteractor2->Start();
}
//for mouse event
vtkStandardNewMacro(MouseInteractorStyle3);

double calcAvg(std::vector<double> vec) {
    double sum;
    for(std::size_t i = 0; i < vec.size(); i++)
       sum += vec.at(i);
    return sum/vec.size();
}

double findPearsonCorrelation(std::vector<double> x, std::vector<double> y) {
  // Find the avg of first and second.
  double xbar = calcAvg(x);
  double ybar = calcAvg(y);
  double a = 0;
  double b = 0;
  double total_a_sq = 0;
  double total_b_sq = 0;
  double total_axb = 0;
  std::cout << "The sizes are " << x.size();
  std::cout << " and " << y.size();
  for(std::size_t i = 0; i < y.size(); i++) {
    a = (xbar - x.at(i));
    b = (ybar - y.at(i));
    total_axb += (a * b);
    total_a_sq += (a * a);
    total_b_sq += (b * b);
  }
  return total_axb/sqrt(total_a_sq * total_b_sq);
}
//This is the main method for the entire project, add your part here
int main(int argc, char *argv[])
{
  ScCalc *calculator = new ScCalc();
  double* spacing1;

  //SEED INPUT GUI
  //TODO: Ken, your class should go here
  //The output should be a double[3]
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " InputFile1\n" << "InputFile2\n";
    return EXIT_FAILURE;
  }

  //Hidden parameter -d allows the program to use pre-calculated data
  if(argc > 2)
  {
    for (int i = 2; i < argc; i++)
    {
      if (strcmp(argv[i], "-d") == 0)
      {
        useDatabase = true;
        std::cout << "-d Use Database" <<endl;
      }
      if (strcmp(argv[i], "-f") == 0)
      {
        fitOnly = true;
        std::cout << "-f Show fit Only" <<endl;
      }
      if (strcmp(argv[i], "-a") == 0) {
        useAnnotations = true;
        std::cout << "-a Using the lml file annotation in the directory" << endl;
      }
      if (strcmp(argv[i], "-c") == 0) {
        useComponents = true;
        std::cout << "-c Decompose into 2D projections" << endl;
      }
    }
  }

  //read input mhd file
  vtkSmartPointer<vtkMetaImageReader>reader =
  vtkSmartPointer<vtkMetaImageReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  spacing1 = reader->GetPixelSpacing();

  vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
  //reslice->SetOutputExtent(0, 9, 0, 100, 0, 0);
  reslice->SetOutputSpacing(1,1,spacing1[2]);
  reslice->SetInputConnection(reader->GetOutputPort());
  reslice->Update();


  //display
  vtkSmartPointer<vtkResliceImageViewer> imageViewer =
  vtkSmartPointer<vtkResliceImageViewer>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor2 = 
  vtkSmartPointer<vtkRenderWindowInteractor>::New();

  // added this
  vtkSmartPointer<vtkImageFlip> flipXFilter = vtkSmartPointer<vtkImageFlip>::New();
  flipXFilter->SetFilteredAxis(1); // flip y axis
  flipXFilter->SetInputConnection(reslice->GetOutputPort()); 
  flipXFilter->Update();


  imageViewer->SetResliceModeToAxisAligned();
  imageViewer->SetInputConnection(flipXFilter->GetOutputPort());
  imageViewer->SetupInteractor(renderWindowInteractor2);
  imageViewer->SetColorLevel(500);
  imageViewer->SetColorWindow(2000);

  //set z coord always the most center slice 
  seedZ = (imageViewer->GetSliceMax())/2;

  imageViewer->SetSize(reader->GetWidth(), reader->GetHeight());
  imageViewer->SetSlice(seedZ);
  //imageViewer->SetSliceOrientationToXY();
  imageViewer->SetSliceOrientationToXY();
  imageViewer->Render();

  seedZ = imageViewer->GetSlice();

  vtkSmartPointer<MouseInteractorStyle3> style =
  vtkSmartPointer<MouseInteractorStyle3>::New();
  
  renderWindowInteractor2->SetInteractorStyle(style);
  //renderWindowInteractor2->UpdateSize(100,100);
  renderWindowInteractor2->Start();

  seedZ = imageViewer->GetSlice();

  double seed1[3] = {seedX, seedY, seedZ};

  std::cout << "Seed1 Set at: " << seed1[0] <<" "<< seed1[1] <<" "<<seed1[2] <<endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  
  //debug log
  std::cout << "Spacing Set at: " << spacing1[0] << " " << spacing1[1] << " " << spacing1[2] <<endl;

  //Check if files are in database
  DIR* dataFolder;
  struct stat st = {0};

  if (useDatabase)
  {
    if (stat("preComputedData", &st) == -1) 
    {
      mkdir("preComputedData", 0777);
    }

    dataFolder = opendir("preComputedData");

    struct dirent *ent; 
    while((ent = readdir(dataFolder)) != NULL) 
    { 
      if (strcmp(argv[1], ent->d_name) == 0)
      {
        newFile1 = false;
        std::cout << "File 1: " << argv[1] << " found in database." << endl;
      }
    } 

  }

  //SEGMENTATION
  std::vector<std::vector<int>> centroids1;
  RegionGrowingNoThreshold region_growing;

  char* fullName1;
  fullName1 = (char*)malloc(strlen(argv[1]) + 18 + 1);
  strcpy(fullName1, "./preComputedData/");
  strcat(fullName1, argv[1]);

  calculator->spacing1 = spacing1;
  if (useDatabase && !newFile1)
  {
    //Load from database
    centroids1 = calculator->loadVector(fullName1);
    std::cout << "IMAGE1 LOADED FROM DATABASE..." <<endl;
  }
  else if (useAnnotations) {
    // Use the centroids from the annotation information instead of performing segmentation.
    // read the lml file.
    centroids1 = calculator->loadAnnotationData(argv[1]);
    std::cout << "LOADING CENTROIDS from .lml file" << endl;
      // Loading the color map
  calculator->loadColorMap(argv[1]);
    std::cout << "LOADING THE VERTEBRAE MAPPING. " << endl;
  }
  else
  {
    //Calculate fresh
    std::cout << "CALCULATING SEGMENTATION IMAGE1..." <<endl;
    centroids1 = region_growing.GetCentroids(argv[1], seed1[0], seed1[1], seed1[2]);
    std::cout << "IMAGE1 SEGMENTATION COMPLETE" <<endl;
    std::cout << "Writing out Image1 with centroids" << endl;
    //region_growing.WriteImageWithCentroids(argv[1], "/tmp/results/output1.mhd", centroids1);
    if (useDatabase)
    {
      calculator->saveVector(centroids1, fullName1);
    }
  }

  //calculator->printVector(centroids1);
  //calculator->printVector(centroids2);
  
  //Hardcoded segmentation output
  double spiral[7][3] = {{0.0, 0.0, 0.0},
  {1.0, 1.0, 0.0},
  {0.5, 2.0, 1.0},
  {0.0, 3.0, 0.0},
  {1.0, 4.0, 0.0},
  {0.5, 5.0, 1.0},
  {0.0, 6.0, 0.0}};
  unsigned spLength = 7;

  double spiral2[7][3] = {{0.0, 1.0, 0.0},
  {1.0, 2.0, 0.0},
  {0.5, 3.0, 1.0},
  {0.0, 4.0, 0.0},
  {1.0, 5.0, 0.0},
  {0.5, 6.0, 1.0},
  {0.0, 7.0, 0.0}};
  unsigned spLength2 = 7;

  //REGISTRATION
  Registration *reg = new Registration();

  //Mutli Res Image Registration
  double (*trans)[4] = new double[4][4];
  //test with 1 iteration of optimizer
  //reg->multiResRegistration(argv[1], argv[2], trans, 1);

  //Rigid 3D Registration
  double trans2[4][4]; // to be populated by registration algorithm
  // test with 1 iteration of optimizer
  //reg->rigidAlign(argv[1], argv[2], trans2, 1); 

  //FINAL RESULT CALCULATION
  //TODO: Juris will need to improve this to produce reasonable results
  calculator->spacing1 = spacing1;
  calculator->loadSpine1(centroids1);
  //calculator->printSpine1();
 // calculator->loadTransofrm(trans);
  //calculator->printTransform();
  //calculator->transformSpine1();
 //calculator->printSpine1();{
  if (useComponents == true) {
    calculator->crateSpineFitX(1);
    calculator->crateSpineFitY(2);
  }
  else {
     calculator->crateSpineFit(1);
  }
 
  calculator->printAngles();
  // calculator->multidimfit();




  //Set colors for spine
  double color1[3] = {1, 0, 0};
  double color2[3] = {0, 1, 0};
  double color3[3] = {0, 0, 1};
  double *colorA  = color3;
  double *colorB  = color3;

  //Create axis
  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  transform->Translate(0.0, 0.0, 0.0);

  vtkSmartPointer<vtkAxesActor> axes =
    vtkSmartPointer<vtkAxesActor>::New();

  axes->SetUserTransform(transform);

  // properties of the axes labels can be set as follows
  // this sets the x axis label to red
  // axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1,0,0);
 
  // the actual text of the axis label can be changed:
  // axes->SetXAxisLabelText("test");
  axes->SetTotalLength(50,50,50);
  axes->AxisLabelsOff();

  // added this
  vtkSmartPointer<vtkCamera> camera =  vtkSmartPointer<vtkCamera>::New();
  camera->SetPosition(0, 0, 100);
  camera->SetFocalPoint(0, 0, 0);
  camera->Roll(90.0);
  camera->Pitch(90.0);
  std::cout<<"The current viewup is " << *(camera->GetViewUp());
  
  // Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
  vtkSmartPointer<vtkRenderer>::New();

  // setup the camera
  renderer->SetActiveCamera(camera);
  vtkSmartPointer<vtkRenderWindow> renderWindow =
  vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderer->AddActor(axes);
  if (!fitOnly)
  {
    // red curve
  	renderer->AddActor(makeLine(calculator->spine1,calculator->spine1Length,color1));
  } else {colorA = color1; colorB = color2;}

  // blue curve
  renderer->AddActor(makeLine(calculator->fit1,calculator->spine1Length,colorA));
  if (useComponents) {
    renderer->AddActor(makeLine(calculator->fit2, calculator->spine1Length, color2));
    renderer->AddActor(drawPointsDouble(calculator->fit2, spacing1, calculator->spine1Length, calculator->vertebraePoints));
    drawImageWithPCurve(argv[1], calculator->fit2, spacing1, calculator->spine1Length, 2, color2, calculator->vertebraePoints);
  }
  renderer->AddActor(makeSpine(argv[1]));
  // Adding points on the main line
  renderer->AddActor(drawPoints(centroids1, spacing1, calculator->vertebraePoints));

  // Adding points on the projected lines.
  renderer->AddActor(drawPointsDouble(calculator->fit1, spacing1, calculator->spine1Length, calculator->vertebraePoints));

  //Ouput final view
  renderWindow->Render();
  renderWindowInteractor->Start();


  // Drawing image with curve
  //drawImageWithCurve(argv[1], centroids1, spacing1);
  drawImageWithPCurve(argv[1], calculator->fit1, spacing1, calculator->spine1Length, 1, colorA, calculator->vertebraePoints);

  // Finding the correlation
  double r_kappa_ca_X = findPearsonCorrelation(calculator->kappaX, calculator->caX);
  std::cout << "The X correlation is " << r_kappa_ca_X << endl;
  double r_kappa_ca_Y = findPearsonCorrelation(calculator->kappaY, calculator->caY);
  std::cout << "The Y correlation is " << r_kappa_ca_Y << endl;
  double r_max_kappa_X = findPearsonCorrelation(calculator->maxAnglesY, calculator->kappaY);
  std::cout << "The maxAngle correlation is " << r_max_kappa_X << endl;
  return EXIT_SUCCESS;
}
