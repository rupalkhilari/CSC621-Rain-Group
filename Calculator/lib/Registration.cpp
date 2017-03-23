#include "Registration.hh"
#include <ctime>

using namespace std;

template <typename TRegistration>

class RegistrationInterfaceCommand : public itk::Command
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );
protected:
  RegistrationInterfaceCommand() {};
public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
  typedef   OptimizerType *                            OptimizerPointer;
  void Execute(itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
      return;
      }
    RegistrationPointer registration = static_cast<RegistrationPointer>( object );
    if(registration == ITK_NULLPTR)
      {
      return;
      }
    OptimizerPointer optimizer = static_cast< OptimizerPointer >(registration->GetModifiableOptimizer() );
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : "
              << registration->GetCurrentLevel()  << std::endl;
    std::cout << std::endl;
    if ( registration->GetCurrentLevel() == 0 )
      {
      optimizer->SetMaximumStepLength( 16.00 );
      optimizer->SetMinimumStepLength( 0.01 );
      }
    else
      {
      optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / 4.0 );
      optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 10.0 );
      }
    }
  void Execute(const itk::Object * , const itk::EventObject & ) ITK_OVERRIDE
    { return; }
};


//  The following section of code implements an observer
//  that will monitor the evolution of the registration process.
//
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef   itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef   const OptimizerType *                    OptimizerPointer;
  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
      Execute( (const itk::Object *)caller, event);
    }
  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
      OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
      if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};

void Registration::multiResRegistration(string fixedImageInput, string movingImageInput, double transformParameters[][4], int maxNumberOfIterations) {

  cout<<"Working on multi-res image registration... (please be patient)"<<endl;

  const    unsigned int    Dimension = 3;
  typedef  unsigned short  PixelType;
  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;
  typedef   float                                    InternalPixelType;
  typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
  typedef itk::TranslationTransform< double, Dimension > TransformType;
  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
  typedef itk::LinearInterpolateImageFunction<
                                    InternalImageType,
                                    double             > InterpolatorType;
  typedef itk::MattesMutualInformationImageToImageMetric<
                                    InternalImageType,
                                    InternalImageType >   MetricType;
  typedef itk::MultiResolutionImageRegistrationMethod<
                                    InternalImageType,
                                    InternalImageType >   RegistrationType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   FixedImagePyramidType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   MovingImagePyramidType;
  //  All the components are instantiated using their \code{New()} method
  //  and connected to the registration object as in previous example.
  //
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  MetricType::Pointer         metric        = MetricType::New();
  FixedImagePyramidType::Pointer fixedImagePyramid =
      FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingImagePyramid =
      MovingImagePyramidType::New();
  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetInterpolator(  interpolator  );
  registration->SetMetric( metric  );
  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );
  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  fixedImageReader->SetFileName(  fixedImageInput );
  movingImageReader->SetFileName( movingImageInput );

  typedef itk::CastImageFilter<
                        FixedImageType, InternalImageType > FixedCastFilterType;
  typedef itk::CastImageFilter<
                        MovingImageType, InternalImageType > MovingCastFilterType;
  FixedCastFilterType::Pointer fixedCaster   = FixedCastFilterType::New();
  MovingCastFilterType::Pointer movingCaster = MovingCastFilterType::New();
  fixedCaster->SetInput(  fixedImageReader->GetOutput() );
  movingCaster->SetInput( movingImageReader->GetOutput() );
  registration->SetFixedImage(    fixedCaster->GetOutput()    );
  registration->SetMovingImage(   movingCaster->GetOutput()   );
  fixedCaster->Update();
  registration->SetFixedImageRegion(
       fixedCaster->GetOutput()->GetBufferedRegion() );
  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );
  initialParameters[0] = 0.0;  // Initial offset in mm along X
  initialParameters[1] = 0.0;  // Initial offset in mm along Y
  initialParameters[2] = 0.0;  // Initial offset in mm along Z
  registration->SetInitialTransformParameters( initialParameters );
  metric->SetNumberOfHistogramBins( 128 );
  metric->SetNumberOfSpatialSamples( 50000 );
  metric->ReinitializeSeed( 76926294 );

  optimizer->SetNumberOfIterations( maxNumberOfIterations );
  optimizer->SetRelaxationFactor( 0.9 );
  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
  CommandType::Pointer command = CommandType::New();
  registration->AddObserver( itk::IterationEvent(), command );
  registration->SetNumberOfLevels( 3 );
  try
    {
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    }

  ParametersType finalParameters = registration->GetLastTransformParameters();
  double TranslationAlongX = finalParameters[0];
  double TranslationAlongY = finalParameters[1];
  double TranslationAlongZ = finalParameters[2];
  unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  double bestValue = optimizer->GetValue();
  // Print out results
  //
  std::cout << "Result = " << std::endl;
  std::cout << " Translation X = " << TranslationAlongX  << std::endl;
  std::cout << " Translation Y = " << TranslationAlongY  << std::endl;
  std::cout << " Translation Z = " << TranslationAlongZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;
  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;
  TransformType::Pointer finalTransform = TransformType::New();
  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );
  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetTransform( finalTransform );
  resample->SetInput( movingImageReader->GetOutput() );
  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  PixelType backgroundGrayLevel = 100;

  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( backgroundGrayLevel );
  typedef  unsigned char  OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter<
                        FixedImageType,
                        OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();
  writer->SetFileName( "Reg_MultiRes_Result.mhd" );

  cout<<endl;
  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();
  //
  // Generate checkerboards before and after registration
  //
  typedef itk::CheckerBoardImageFilter< FixedImageType > CheckerBoardFilterType;
  CheckerBoardFilterType::Pointer checker = CheckerBoardFilterType::New();
  checker->SetInput1( fixedImage );
  checker->SetInput2( resample->GetOutput() );
  caster->SetInput( checker->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  resample->SetDefaultPixelValue( 0 );

  // Before registration
  TransformType::Pointer identityTransform = TransformType::New();
  identityTransform->SetIdentity();
  resample->SetTransform( identityTransform );
  writer->SetFileName( "Reg_Tiled_Before.mhd" );
  writer->Update();

  // After registration
  resample->SetTransform( finalTransform );
  writer->SetFileName( "Reg_Tiled_After.mhd" );
  writer->Update();

  // interpret final transformation parameters as 4x4 matrix
  // tokenize the parameters
  stringstream ss;
  ss << finalParameters << endl;
  string paramsString = ss.str();
  paramsString = paramsString.substr(1, paramsString.length() - 3);
  vector<double> params_v;
  char delim = ',';
  size_t start = paramsString.find_first_not_of(delim), end = start;
  while (start != string::npos)
  {
      end = paramsString.find(delim, start);
      params_v.push_back(atof(paramsString.substr(start, end - start).c_str()));
      start = paramsString.find_first_not_of(delim, end);
  }

  // populate the translation matrix
  transformParameters[0][0] = params_v[0];
  transformParameters[1][0] = params_v[1];
  transformParameters[2][0] = params_v[2];

  transformParameters[0][1] = params_v[3];
  transformParameters[1][1] = params_v[4];
  transformParameters[2][1] = params_v[5];

  transformParameters[0][2] = params_v[6];
  transformParameters[1][2] = params_v[7];
  transformParameters[2][2] = params_v[8];

  transformParameters[3][0] = params_v[9];
  transformParameters[3][1] = params_v[10];
  transformParameters[3][2] = params_v[11];

  transformParameters[0][3] = transformParameters[1][3]
          = transformParameters[2][3] = 0.0;
  transformParameters[3][3] = 1.0;

  // display the matrix
  cout << "   Transform matrix:" << endl;
  for (int i = 0; i < 4; i++)
  {
      for (int j = 0; j < 4; j++)
          cout << setw(15) << transformParameters[j][i] << " ";
      cout << endl;
  }
}

void Registration::rigidAlign(string fixedImageInput, string movingImageInput, double transformParameters[][4], int maxNumberOfIterations) {

  cout << "Starting Rigid Alignment now... (please be patient)" << endl;

  //Metrics for time
  clock_t begint = clock();

  const unsigned int                          Dimension = 3;
  typedef  float                              PixelType;
  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  typedef itk::VersorRigid3DTransform< double > TransformType;

  typedef itk::RegularStepGradientDescentOptimizerv4<double>    OptimizerType;
  typedef itk::MeanSquaresImageToImageMetricv4<
                                            FixedImageType,
                                            MovingImageType >   MetricType;
  typedef itk::ImageRegistrationMethodv4<
                                      FixedImageType,
                                      MovingImageType,
                                      TransformType >           RegistrationType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );

  TransformType::Pointer  initialTransform = TransformType::New();

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  fixedImageInput );
  movingImageReader->SetFileName( movingImageInput );

  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );

  typedef itk::CenteredTransformInitializer<
    TransformType,
    FixedImageType,
    MovingImageType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer =
    TransformInitializerType::New();

  initializer->SetTransform(   initialTransform );
  initializer->SetFixedImage(  fixedImageReader->GetOutput() );
  initializer->SetMovingImage( movingImageReader->GetOutput() );

  initializer->MomentsOn();

  initializer->InitializeTransform();

  typedef TransformType::VersorType  VersorType;
  typedef VersorType::VectorType     VectorType;
  VersorType     rotation;
  VectorType     axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  const double angle = 0;
  rotation.Set(  axis, angle  );
  initialTransform->SetRotation( rotation );

  registration->SetInitialTransform( initialTransform );

  // set algorithmic parameters here; these could be passed as function
  // arguments instead
  const double translationScale = 1.0 / 1000.0;
  double learningRate = 0.2;
  double minStepLength = 0.001;

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( initialTransform->GetNumberOfParameters() );
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  optimizer->SetScales( optimizerScales );
  optimizer->SetNumberOfIterations( maxNumberOfIterations );
  optimizer->SetLearningRate( learningRate );
  optimizer->SetMinimumStepLength( minStepLength );
  optimizer->SetReturnBestParametersAndValue(true);

  // cout << endl << endl;
  // cout << "Parameters = " << endl;
  // cout << " optimizerScales  = " << optimizerScales << endl;
  // cout << " translationScale = " << translationScale << endl;
  // cout << " learningRate     = " << learningRate << endl;
  // cout << " minStepLength    = " << minStepLength << endl << endl;
  // cout << " Output from optimizer iterations:" << endl;

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  const unsigned int numberOfLevels = 1;

  RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize( 1 );
  shrinkFactorsPerLevel[0] = 1;

  RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize( 1 );
  smoothingSigmasPerLevel[0] = 0;

  registration->SetNumberOfLevels( numberOfLevels );
  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

  try
    {
    registration->Update();
    cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << endl;
    }
  catch( itk::ExceptionObject & err )
    {
    cerr << "ExceptionObject caught !" << endl;
    cerr << err << endl;
    return;
    }

  const TransformType::ParametersType finalParameters =
                            registration->GetOutput()->Get()->GetParameters();

  const double versorX              = finalParameters[0];
  const double versorY              = finalParameters[1];
  const double versorZ              = finalParameters[2];
  const double finalTranslationX    = finalParameters[3];
  const double finalTranslationY    = finalParameters[4];
  const double finalTranslationZ    = finalParameters[5];
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const double bestValue = optimizer->GetValue();

  // Print out results
  cout << endl << endl;
  cout << "Result = " << endl;
  cout << " versor X      = " << versorX  << endl;
  cout << " versor Y      = " << versorY  << endl;
  cout << " versor Z      = " << versorZ  << endl;
  cout << " Translation X = " << finalTranslationX  << endl;
  cout << " Translation Y = " << finalTranslationY  << endl;
  cout << " Translation Z = " << finalTranslationZ  << endl;
  cout << " Iterations    = " << numberOfIterations << endl;
  cout << " Metric value  = " << bestValue          << endl;

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetFixedParameters( registration->GetOutput()->Get()->GetFixedParameters() );
  finalTransform->SetParameters( finalParameters );

  TransformType::MatrixType matrix = finalTransform->GetMatrix();
  TransformType::OffsetType offset = finalTransform->GetOffset();
  // cout << endl << "Matrix = " << endl << matrix << endl;
  // cout << "Offset = " << offset << endl;

  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetTransform( finalTransform );
  resampler->SetInput( movingImageReader->GetOutput() );

  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  resampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 100 );

  typedef  unsigned char                                          OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension >                OutputImageType;
  typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >                 WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( "Reg_Rigid_Result.mhd" );

  cout<<endl;

  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  typedef itk::SubtractImageFilter<
                                  FixedImageType,
                                  FixedImageType,
                                  FixedImageType > DifferenceFilterType;
  DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

  typedef itk::RescaleIntensityImageFilter<
                                  FixedImageType,
                                  OutputImageType >   RescalerType;
  RescalerType::Pointer intensityRescaler = RescalerType::New();

  intensityRescaler->SetInput( difference->GetOutput() );
  intensityRescaler->SetOutputMinimum(   0 );
  intensityRescaler->SetOutputMaximum( 255 );

  difference->SetInput1( fixedImageReader->GetOutput() );
  difference->SetInput2( resampler->GetOutput() );

  resampler->SetDefaultPixelValue( 1 );

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput( intensityRescaler->GetOutput() );

  // Compute the difference image between the
  // fixed and resampled moving image.
  writer2->SetFileName( "Reg_Diff_Before.mhd" );
  writer2->Update();

  typedef itk::IdentityTransform< double, Dimension > IdentityTransformType;
  IdentityTransformType::Pointer identity = IdentityTransformType::New();

  // Compute the difference image between the
  // fixed and moving image before registration.
  resampler->SetTransform( identity );
  writer2->SetFileName( "Reg_Diff_After.mhd" );
  writer2->Update();

  clock_t endt = clock();
  double elapsed_secs = double(endt - begint) / CLOCKS_PER_SEC;

  // interpret final transformation parameters as 4x4 matrix
  // tokenize the parameters
  stringstream ss;
  ss << finalParameters << endl;
  string paramsString = ss.str();
  paramsString = paramsString.substr(1, paramsString.length() - 3);
  vector<double> params_v;
  char delim = ',';
  size_t start = paramsString.find_first_not_of(delim), end = start;
  while (start != string::npos)
  {
      end = paramsString.find(delim, start);
      params_v.push_back(atof(paramsString.substr(start, end - start).c_str()));
      start = paramsString.find_first_not_of(delim, end);
  }

  // populate the translation matrix
  transformParameters[0][0] = params_v[0];
  transformParameters[1][0] = params_v[1];
  transformParameters[2][0] = params_v[2];

  transformParameters[0][1] = params_v[3];
  transformParameters[1][1] = params_v[4];
  transformParameters[2][1] = params_v[5];

  transformParameters[0][2] = params_v[6];
  transformParameters[1][2] = params_v[7];
  transformParameters[2][2] = params_v[8];

  transformParameters[3][0] = params_v[9];
  transformParameters[3][1] = params_v[10];
  transformParameters[3][2] = params_v[11];

  transformParameters[0][3] = transformParameters[1][3]
          = transformParameters[2][3] = 0.0;
  transformParameters[3][3] = 1.0;

  // display the matrix
  cout << "   Transform matrix:" << endl;
  for (int i = 0; i < 4; i++)
  {
      for (int j = 0; j < 4; j++)
          cout << setw(15) << transformParameters[j][i] << " ";
      cout << endl;
  }

  cout<<"The time elapsed was: " << elapsed_secs << endl;
}

/**
 * Performs 3D affine registration using mean squares, with a default max
 * number of iterations of 300.
 * @param fixedImageInput        - a metaimage header
 * @param movingImageInput       - a metaimage header
 * @param transformParameters    - array[4][4] to write transform matrix
 */
void Registration::affineAlign(string fixedImageInput, string movingImageInput,
        double transformParameters[][4])
{
    return Registration::affineAlign(fixedImageInput, movingImageInput,
            transformParameters, 1);
}

/**
 * Performs 3D affine registration using mean squares.
 * @param fixedImageInput        - a metaimage header
 * @param movingImageInput       - a metaimage header
 * @param transformParameters    - array[4][4] to write transform matrix
 * @param maxNumberOfIterations  - max number of iterations; default is 300
 */
 void Registration::affineAlign(string fixedImageInput, string movingImageInput,
         double transformParameters[][4], int maxNumberOfIterations)
 {
 // Define image types
    const unsigned int Dimension = 3;
    typedef int16_t PixelType; // this determines the size of a pixel
    typedef itk::Image< PixelType, Dimension > FixedImageType;
    typedef itk::Image< PixelType, Dimension > MovingImageType;

	//  The transform type is instantiated using the code below. The template
    //  parameters of this class are the representation type of the space
    //  coordinates and the space dimension.
    typedef itk::AffineTransform< double, Dimension > TransformType;
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::MeanSquaresImageToImageMetric< FixedImageType, MovingImageType > MetricType;
    typedef itk::LinearInterpolateImageFunction< MovingImageType, double > InterpolatorType;
    typedef itk::ImageRegistrationMethod< FixedImageType, MovingImageType > RegistrationType;
    MetricType::Pointer metric = MetricType::New();
    OptimizerType::Pointer optimizer = OptimizerType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    RegistrationType::Pointer registration = RegistrationType::New();
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetInterpolator(interpolator);

	//  The transform object is constructed below and passed to the registration
    //  method.
    TransformType::Pointer transform = TransformType::New();
    registration->SetTransform(transform);

	// Initialize readers for the two input files
	typedef itk::ImageFileReader< FixedImageType > FixedImageReaderType;
    itk::MetaImageIO::Pointer fixedImageIO = itk::MetaImageIO::New();
    FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
    fixedImageReader->SetImageIO(fixedImageIO);
    fixedImageReader->SetFileName(fixedImageInput);
    fixedImageReader->SetUseStreaming(true);
    fixedImageIO->SetUseStreamedReading(true);

    typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
    itk::MetaImageIO::Pointer movingImageIO = itk::MetaImageIO::New();
    MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
    movingImageReader->SetImageIO(movingImageIO);
    movingImageReader->SetFileName(movingImageInput);
    movingImageReader->SetUseStreaming(true);
    movingImageIO->SetUseStreamedReading(true);


	// configure registration algorithm
	registration->SetFixedImage(fixedImageReader->GetOutput());
    registration->SetMovingImage(movingImageReader->GetOutput());
    fixedImageReader->Update();
    registration->SetFixedImageRegion(fixedImageReader->GetOutput()->GetBufferedRegion());

    // Get & print useful statistics about the input images
    //fixedImageReader->GenerateOutputInformation();
    //movingImageReader->GenerateOutputInformation();
    FixedImageType::SizeType fixedImageSize;
    MovingImageType::SizeType movingImageSize;
    fixedImageSize = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
    movingImageSize = movingImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

    cout << "dimensions of input images:\n";
    cout << "   fixed  [x]:" << fixedImageSize[0] << " [y]:" << fixedImageSize[1]
            << " [z]:" << fixedImageSize[2] << "\n"; // dimensions of input 1
    cout << "   moving [x]:" << movingImageSize[0] << " [y]:" << movingImageSize[1]
            << " [z]:" << movingImageSize[2] << "\n"; // dimensions of input 2

	 //  In this example, we again use the
    //  CenteredTransformInitializer helper class in order to compute
    //  a reasonable value for the initial center of rotation and the
    //  translation. The initializer is set to use the center of mass of each
    //  image as the initial correspondence correction.
    typedef itk::CenteredTransformInitializer< TransformType, FixedImageType,
            MovingImageType > TransformInitializerType;
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTransform(transform);
    initializer->SetFixedImage(fixedImageReader->GetOutput());
    initializer->SetMovingImage(movingImageReader->GetOutput());
    initializer->MomentsOn();
    initializer->InitializeTransform();

    //  Now we pass the parameters of the current transform as the initial
    //  parameters to be used when the registration process starts.
    registration->SetInitialTransformParameters(transform->GetParameters());

    //  Keeping in mind that the scale of units in scaling, rotation and
    //  translation are quite different, we take advantage of the scaling
    //  functionality provided by the optimizers. We know that the first $N
    //  \times N$ elements of the parameters array correspond to the rotation
    //  matrix factor, and the last $N$ are the components of the translation to
    //  be applied after multiplication with the matrix is performed.
    double translationScale = 1.0 / 1000.0;

	// configure optimizer
	typedef OptimizerType::ScalesType OptimizerScalesType;
    OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
    optimizerScales[0] = 1.0;
    optimizerScales[1] = 1.0;
    optimizerScales[2] = 1.0;
    optimizerScales[3] = 1.0;
    optimizerScales[4] = 1.0;
    optimizerScales[5] = 1.0;
    optimizerScales[6] = 1.0;
    optimizerScales[7] = 1.0;
    optimizerScales[8] = 1.0;
    optimizerScales[9] = translationScale;
    optimizerScales[10] = translationScale;
    optimizerScales[11] = translationScale;
    optimizer->SetScales(optimizerScales);

 //  We also set the usual parameters of the optimization method. In this
    //  case we are using an RegularStepGradientDescentOptimizer. Below, we
    //  define the optimization parameters like initial step length, minimal
    //  step length and number of iterations. These last two act as stopping
    //  criteria for the optimization.
    double steplength = 0.1;

  optimizer->SetMaximumStepLength(steplength);
    optimizer->SetMinimumStepLength(0.0001);
    optimizer->SetNumberOfIterations(maxNumberOfIterations);
    optimizer->MinimizeOn(); // set the optimizer to do minimization

    // Create the Command observer and register it with the optimizer.
	// TODO: We might not need this.
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);

	// Trigger execution of the registration method by calling the Update()
    // method.
    try
    {
        cout << "\nExecuting registration. This will take awhile."
                << endl;
        registration->Update();
        cout << "\nFinished registration.\nOptimizer stop condition: "
                << registration->GetOptimizer()->GetStopConditionDescription()
                << endl;
    } catch (itk::ExceptionObject & err)
    {
        cerr << "ExceptionObject caught !" << endl;
        cerr << err << endl;
        return;
    }

    // After completion of optimization, recover parameters from the
    // registration method using GetLastTransformParameters().
    OptimizerType::ParametersType finalParameters =
            registration->GetLastTransformParameters();
    const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
    const double bestValue = optimizer->GetValue();

    // Print results
    cout << "Result:" << endl;
    cout << "   Iterations   = " << numberOfIterations << endl;
    cout << "   Metric value = " << bestValue << endl;

    // interpret final transformation parameters as 4x4 matrix
    // tokenize the parameters
    stringstream ss;
    ss << finalParameters << endl;
    string paramsString = ss.str();
    paramsString = paramsString.substr(1, paramsString.length() - 3);
    vector<double> params_v;
    char delim = ',';
    size_t start = paramsString.find_first_not_of(delim), end = start;
    while (start != string::npos)
    {
        end = paramsString.find(delim, start);
        params_v.push_back(atof(paramsString.substr(start, end - start).c_str()));
        start = paramsString.find_first_not_of(delim, end);
    }

    // populate the translation matrix
    transformParameters[0][0] = params_v[0];
    transformParameters[1][0] = params_v[1];
    transformParameters[2][0] = params_v[2];

    transformParameters[0][1] = params_v[3];
    transformParameters[1][1] = params_v[4];
    transformParameters[2][1] = params_v[5];

    transformParameters[0][2] = params_v[6];
    transformParameters[1][2] = params_v[7];
    transformParameters[2][2] = params_v[8];

    transformParameters[3][0] = params_v[9];
    transformParameters[3][1] = params_v[10];
    transformParameters[3][2] = params_v[11];

    transformParameters[0][3] = transformParameters[1][3]
            = transformParameters[2][3] = 0.0;
    transformParameters[3][3] = 1.0;

    // display the matrix
    cout << "   Transform matrix:" << endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
            cout << setw(15) << transformParameters[j][i] << " ";
        cout << endl;
    }
}
