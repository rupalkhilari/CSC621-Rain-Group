#include "scCalc.hh" 
#include <map>
#include <algorithm>

using namespace std;



void ScCalc::printVector(vector<vector<int> > vec)
{
  unsigned length = vec.size();
  for (unsigned i = 0; i < length; ++i)
  {
      cout << "Point " << i << ": [";
      for (unsigned j = 0; j < 3; ++j)
      {
          cout << vec[i][j];
          if (j != 2) cout << ", ";
      }
      cout << "]" << endl;
  }
}

void ScCalc::saveVector(vector<vector<int> > vec, char* fileName)
{
  unsigned length = vec.size();
  ofstream f(fileName);
  for (unsigned i = 0; i < length; ++i)
  {
      for (unsigned j = 0; j < 3; ++j)
      {
          f << vec[i][j];
          if (j != 2) f << ',';
      }
      f << '\n';
  }
}

vector<vector<int>> ScCalc::loadVector(char* fileName)
{
    FILE *fp;
    vector<vector<int> > vec;

    string sFileName = fileName;
    ifstream fileStream(sFileName);
    if (!fileStream.is_open())
    {
        cout << "Exiting unable to open file " << fileName << endl;
    }

    string line;
    while(getline(fileStream, line, '\n')) {
        stringstream ss(line);
        vector<int> numbers;
        string in_line;
        while(getline (ss, in_line, ',')) 
        {
          int i = stoi(in_line, 0);
          numbers.push_back(i);
        }
        vec.push_back(numbers);
    }

    return vec;
}

void ScCalc::loadSpine1(vector<vector<int>> spine)
{
    unsigned length = spine.size();
    double (*newSpine)[3] = new double[length][3];
    for (unsigned i = 0; i < length; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
        {
            newSpine[i][j] = spacing1[j] * spine[i][j];
         }
    }

    spine1 = newSpine; 
    spine1Length = length;
}


void ScCalc::loadSpine1(double spine[][3], unsigned length)
{
    loadSpineX(spine, length, 0);
}

void ScCalc::loadSpineX(double spine[][3], unsigned length, unsigned spineNumber)
{
    double (*newSpine)[3] = new double[length][3];
    double* spacing;

    if (!spineNumber) spacing = spacing1;

    for (unsigned i = 0; i < length; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
        {
            newSpine[i][j] = spacing[j] * spine[i][j];
        }
    }

    if (!spineNumber) {spine1 = newSpine; spine1Length = length;}
}

void ScCalc::loadTransofrm(double matrix[][4])
{

    for (unsigned i = 0; i < 4; ++i)
    {
        for (unsigned j = 0; j < 4; ++j)
        {
            transMatrix[i][j] = matrix[j][i];
        }
    }
}

void ScCalc::transformSpine1()
{
    double (*newSpine)[3] = new double[spine1Length][3];

    for (unsigned i = 0; i < spine1Length; ++i)
    {
        for(int j=0; j<3; ++j) 
        {
            newSpine[i][j] = 0;

            for(int k=0; k<3; ++k)
            {
                newSpine[i][j]+= (spine1[i][k]) * (transMatrix[j][k]);
            }

            newSpine[i][j]+=transMatrix[j][3];
        }
    }

    spine1 = newSpine;
}

void ScCalc::compareSpines()
{
    double (*ptr)[3] = spine1;

    for (unsigned i = 0; i < spine1Length - 9; ++i)
    {
        cout << "derivative = " << curveDerivative(ptr);
        cout << " second derivative = " << curveSecDerivative(ptr++);
        cout << endl;
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            cout << transMatrix[i][j] << " ";
        }

        cout << endl;
    }
}

void ScCalc::printSpine1()
{
    for (unsigned i = 0; i < spine1Length; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            cout << spine1[i][j] << " ";
        }
        cout << endl;
    }
}

void ScCalc::printTransform()
{
    for (unsigned i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            cout << transMatrix[i][j] << " ";
        }
        cout << endl;
    }
}

void ScCalc::printAngles()
{
    cout << endl;
    cout << "=======================================" << endl;
    cout << "= Spine Curvature Summary (Degrees)   =" << endl;
    cout << "=======================================" << endl;
    double max1 = 0;

    unsigned length;
    length = spine1Length;
    int counter = 0;
    for (unsigned i = 0; i < length - 3; ++i)
    {
        cout << ++counter <<" ";
        if (i < spine1Length - 3) 
        {
            cout << "Spine1: " << angles1[i] << ", ";
            if (angles1[i] > max1) max1 = angles1[i];
        } else cout << "               ";
        cout << endl;
    }
    cout << "=======================================" << endl;
    cout << "Spine1 MAX: " << max1 << endl;
}

double ScCalc::firstDerivative(double a[5])
{
    double firstDer(0.0);
    double t(0.001);

    // no loop here as 5 points is just enough to calculate the derivative
    // of one elem (the middle one). If you do use a loop, it has to start
    // at 2 and end 2 elems before the end (so i - 2 and i + 2 are always
    // valid indexes for the array.)

    size_t i(2);

    firstDer = (   4.0 / 3.0 * (a[i + 1] - a[i - 1]) / (2.0 * t)
                 - 1.0 / 3.0 * (a[i + 2] - a[i - 2]) / (4.0 * t) );
    // Rather than use (double)4 you can just use 4.0. The .0 tells the
    // compiler you mean a double (if you need a float instead, append
    // an f, like 3.0f )

    return firstDer;
}

double ScCalc::curveDerivative(double points[5][3])
{
	double *x = new double[5];
	double *y = new double[5];
	double *z = new double[5];

	for (int i = 0; i < 5; ++i)
	{
		x[i] = points[i][0];
		y[i] = points[i][1];
		z[i] = points[i][2];
	}

	double dX = firstDerivative(x);
	double dY = firstDerivative(y);
	double dZ = firstDerivative(z);

	return sqrt(dX * dX + dY * dY + dZ * dZ);
}

double ScCalc::curveSecDerivative(double points[9][3])
{
	double *dOne = new double[9];

	double (*ptr)[3] = points;

	for (int i = 0; i < 5; ++i)
	{
		dOne[i] = curveDerivative(ptr++);
	}

	return firstDerivative(dOne);
}

// std::vector<T> polyfit( const std::vector<T>& oX, 
//     const std::vector<T>& oY, int nDegree )
// {
//     using namespace boost::numeric::ublas;
 
//     if ( oX.size() != oY.size() )
//         throw std::invalid_argument( "X and Y vector sizes do not match" );
 
//     // more intuative this way
//     nDegree++;
    
//     size_t nCount =  oX.size();
//     matrix<T> oXMatrix( nCount, nDegree );
//     matrix<T> oYMatrix( nCount, 1 );
    
//     // copy y matrix
//     for ( size_t i = 0; i < nCount; i++ )
//     {
//         oYMatrix(i, 0) = oY[i];
//     }
 
//     // create the X matrix
//     for ( size_t nRow = 0; nRow < nCount; nRow++ )
//     {
//         T nVal = 1.0f;
//         for ( int nCol = 0; nCol < nDegree; nCol++ )
//         {
//             oXMatrix(nRow, nCol) = nVal;
//             nVal *= oX[nRow];
//         }
//     }
 
//     // transpose X matrix
//     matrix<T> oXtMatrix( trans(oXMatrix) );
//     // multiply transposed X matrix with X matrix
//     matrix<T> oXtXMatrix( prec_prod(oXtMatrix, oXMatrix) );
//     // multiply transposed X matrix with Y matrix
//     matrix<T> oXtYMatrix( prec_prod(oXtMatrix, oYMatrix) );
 
//     // lu decomposition
//     permutation_matrix<int> pert(oXtXMatrix.size1());
//     const std::size_t singular = lu_factorize(oXtXMatrix, pert);
//     // must be singular
//     BOOST_ASSERT( singular == 0 );
 
//     // backsubstitution
//     lu_substitute(oXtXMatrix, pert, oXtYMatrix);
 
//     // copy the result to coeff
//     return std::vector<T>( oXtYMatrix.data().begin(), oXtYMatrix.data().end() );
// }

// added this
void ScCalc::getGeometricCurvature(double points[][3], unsigned spLength, unsigned order, double *xStore, double *yStore, double *dz) {
    // we have the coefficients. Compute the first derivative at the point
    double (*fder)[3];
    double (*sder)[3];
    fder = new double[spLength][3];
    sder = new double[spLength][3];
    //double *kappa;
    cout << "Calculating the geometric curvature !!! " << endl;
    for (int i = 0; i < spLength; ++i) 
    {
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;
        for (int j = order; j > 1; j--) 
        {
            sumX += (j * xStore[j]) * pow(dz[i], j-1);
            sumY += (j * yStore[j]) * pow(dz[i], j-1);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;
        fder[i][0] = sumX;
        fder[i][1] = sumY;
        fder[i][2] = dz[i]; // verify this
    }

    // calculate the second derivative
    for (unsigned i = 0; i < spLength; ++i) 
    {
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;

        for (unsigned j = order; j > 3; --j) 
        {
            sumX += (j * (j-1) * xStore[j]) * pow(dz[i], j-2);
            sumY += (j * (j-1) * yStore[j]) * pow(dz[i], j-2);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;
        sder[i][0] = sumX;
        sder[i][1] = sumY;
        sder[i][2] = dz[i]; // verify this
    }

    // Calcuating the kappa value at each point
     // calculate the second derivative
    for (unsigned i = 0; i < spLength; ++i) 
    {
        cout << sqrt(pow((sder[i][2] * fder[i][1]) - (sder[i][1] * fder[i][2]), 2) +
            pow((sder[i][0] * fder[i][2]) - (sder[i][2] * fder[i][0]), 2) + 
            pow((sder[i][1] * fder[i][0]) - (sder[i][0] * fder[i][1]), 2)) /
            pow(pow(fder[i][0], 2) + pow(fder[i][1], 2) + pow(fder[i][2], 2), 1.5);
        cout << endl;
    }   

}

// Added this
// Assuming parameter z
void ScCalc::getGeometricCurvatureMod(double points[][3], unsigned spLength, unsigned order, double *xStore, double *yStore, double *dz, int type) {
    // we have the coefficients. Compute the first derivative at the point
    double (*fder)[3];
    double (*sder)[3];
    fder = new double[spLength][3];
    sder = new double[spLength][3];
    //double *kappa;
    cout << "Calculating the geometric curvature !!! " << endl;
    int component = 0;
    if (xStore[0] == xStore[1] && xStore[1] == xStore[2]) {
        component = 1;
    }
    if (yStore[0] == yStore[1] && yStore[1] == yStore[2]) {
        component = 2;
    }

    for (int i = 0; i < spLength; ++i) 
    {
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;
        for (int j = order; j > 1; j--) 
        {
            sumX += (j * xStore[j]) * pow(dz[i], j-1);
            sumY += (j * yStore[j]) * pow(dz[i], j-1);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;
        if (component == 1) {
            fder[i][0] = 1;
        }
        else {
            fder[i][0] = sumX;
        }
        fder[i][1] = sumY;
        fder[i][2] = dz[i]; // verify this
    }

    // calculate the second derivative
    for (unsigned i = 0; i < spLength; ++i) 
    {
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;

        for (unsigned j = order; j > 3; --j) 
        {
            sumX += (j * (j-1) * xStore[j]) * pow(dz[i], j-2);
            sumY += (j * (j-1) * yStore[j]) * pow(dz[i], j-2);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;
        if (component == 1) {
            sder[i][0] = 0;
        }
        else {
            sder[i][0] = sumX;
        }
        sder[i][1] = sumY;
        sder[i][2] = dz[i]; // verify this
    }

    // Calcuating the kappa value at each point
    cout << "Composite formula : " << endl;
    for (unsigned i = 0; i < spLength; i++) {
        double kappa = (fder[i][0]*sder[i][1]) - (fder[i][1]*sder[i][0]) / 
        pow(pow(fder[i][0], 2) + pow(fder[i][1], 2), 1.5);
        cout << kappa << endl;
    }

    cout << "X component formula : " << endl;
    for (unsigned i = 0; i < spLength; i++) {
        double kappa = abs(sder[i][0]) / pow(pow(fder[i][0], 2) + 1, 1.5);
        if (type == 1) {
            cout << kappa << endl;
            kappaX.push_back(kappa);
        }
    }

    cout << "Y component formula : " << endl;
    for (unsigned i = 0; i < spLength; i++) {
        double kappa = abs(sder[i][1]) / pow(pow(fder[i][1], 2) + 1, 1.5);
        if (type == 2) {
            cout << kappa << endl;
            kappaY.push_back(kappa);
        }
    }

    // 
    /*for (unsigned i = 0; i < spLength; ++i) 
    {
        cout << sqrt(pow((sder[i][2] * fder[i][1]) - (sder[i][1] * fder[i][2]), 2) +
            pow((sder[i][0] * fder[i][2]) - (sder[i][2] * fder[i][0]), 2) + 
            pow((sder[i][1] * fder[i][0]) - (sder[i][0] * fder[i][1]), 2)) /
            pow(pow(fder[i][0], 2) + pow(fder[i][1], 2) + pow(fder[i][2], 2), 1.5);
        cout << endl;
    } */    
}

void ScCalc::computeGeometricCurvature2D(double points[][3], unsigned spLength, unsigned order, double *coefficients, double *x) {
       // we have the coefficients. Compute the first derivative at the point
    double (*fder)[2];
    double (*sder)[2];
    fder = new double[spLength][2];
    sder = new double[spLength][2];
    //double *kappa;
    cout << "Calculating the new geometric curvature !!! " << endl;

    /* Considering the parameteric equations in terms of t, we substitule x = t,
    thus fder[i][0] = 1 and sder[i][0] = 0. Taking derivative wrt t
    */
    for (int i = 0; i < spLength; ++i) 
    {
        double sumY = 0;
        for (int j = order-1; j > 1; j--) 
        {
            sumY += ((j * coefficients[j]) * pow(x[i], j-1));
        }
        fder[i][0] = 1;
        fder[i][1] = sumY;
    }

    // calculate the second derivative
    for (unsigned i = 0; i < spLength; ++i) 
    {
        double sumY = 0;

        for (unsigned j = order-1; j > 1; --j) 
        {
            sumY += (j * (j-1) * coefficients[j]) * pow(x[i], j-2);
        }
        sder[i][0] = 0;
        sder[i][1] = sumY;
    }

    // Calcuating the kappa value at each point
    double *kappas = new double[spLength];
    cout << "Composite formula : " << endl;
    for (unsigned i = 0; i < spLength; i++) {
        double kappa = abs(fder[i][0]*sder[i][1]) - (fder[i][1]*sder[i][0]) / 
        pow(pow(fder[i][0], 2) + pow(fder[i][1], 2), 1.5);
        cout << kappa << endl;
        kappas[i] = kappa;
    }

    cout << "The new curvature angles are : " << endl;
    for (unsigned i = 0; i < spLength-1; i++) {
        cout << getCurvatureAngle(points[i][0], points[i][1], 0, points[i+1][0], points[i+1][1], 0, kappas[i]) << endl;
    }

} 
// Assuming parameter z
void ScCalc::getCurvatureAngleMod(double points[][3], unsigned spLength, unsigned order, double *xStore, double *yStore, double *dz, int type) {
    // we have the coefficients. Compute the first derivative at the point
    double (*fder)[3];
    double (*sder)[3];
    fder = new double[spLength][3];
    sder = new double[spLength][3];
    //double *kappa;
    cout << "Calculating the geometric curvature !!! " << endl;
    for (int i = 0; i < spLength-1; ++i) 
    {
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;
        for (int j = order; j > 1; j--) 
        {
            sumX += (j * xStore[j]) * pow((dz[i]+dz[i+1])/2, j-1);
            sumY += (j * yStore[j]) * pow((dz[i]+dz[i+1])/2, j-1);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;
        fder[i][0] = sumX;
        fder[i][1] = sumY;
        fder[i][2] = dz[i]; // verify this
    }

    // calculate the second derivative
    for (unsigned i = 0; i < spLength-1; ++i) 
    {
        double sumX = 0;
        double sumY = 0;
        double sumZ = 0;

        for (unsigned j = order; j > 3; --j) 
        {
            sumX += (j * (j-1) * xStore[j]) * pow((dz[i]+dz[i+1])/2, j-2);
            sumY += (j * (j-1) * yStore[j]) * pow((dz[i]+dz[i+1])/2, j-2);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;
        sder[i][0] = sumX;
        sder[i][1] = sumY;
        sder[i][2] = dz[i]; // verify this
    }

    // Calcuating the kappa value at each point
    /*cout << "Composite formula : " << endl;
    for (unsigned i = 0; i < spLength-1; i++) {
        double kappa = abs((fder[i][0]*sder[i][1]) - (fder[i][1]*sder[i][0])) / 
        pow(pow(fder[i][0], 2) + pow(fder[i][1], 2), 1.5);
        cout << "At point " << i << ", K = " << kappa << endl;
        double ca = getCurvatureAngle(points[i][0], points[i][1], points[i][2], points[i+1][0], points[i+1][1], points[i+1][2], kappa);
        cout << " curvature angle = " << ca << endl;
    }

    cout << "X component formula : " << endl;
    for (unsigned i = 0; i < spLength-1; i++) {
        double kappa = abs(sder[i][0]) / pow(pow(fder[i][0], 2) + 1, 1.5);
        cout << "At point " << i << ", KX = " << kappa << endl;
        double ca = getCurvatureAngle(points[i][0], points[i][1], points[i][2], points[i+1][0], points[i+1][1], points[i+1][2], kappa);
        cout << " curvature angle = " << ca << endl;
        if (type == 1) {
            caX.push_back(ca);
        }
    }

    cout << "Y component formula : " << endl;
    for (unsigned i = 0; i < spLength-1; i++) {
        double kappa = abs(sder[i][1]) / pow(pow(fder[i][1], 2) + 1, 1.5);
        cout << "At point " << i << ", KY = " << kappa << endl;
        double ca = getCurvatureAngle(points[i][0], points[i][1], points[i][2], points[i+1][0], points[i+1][1], points[i+1][2], kappa);
        cout << " curvature angle = " << ca << endl;
        if (type == 2) {
            caY.push_back(ca);
        }
    }*/

    // 
    /*for (unsigned i = 0; i < spLength; ++i) 
    {
        cout << sqrt(pow((sder[i][2] * fder[i][1]) - (sder[i][1] * fder[i][2]), 2) +
            pow((sder[i][0] * fder[i][2]) - (sder[i][2] * fder[i][0]), 2) + 
            pow((sder[i][1] * fder[i][0]) - (sder[i][0] * fder[i][1]), 2)) /
            pow(pow(fder[i][0], 2) + pow(fder[i][1], 2) + pow(fder[i][2], 2), 1.5);
        cout << endl;
    } */    
}

double ScCalc::getCurvatureAngle(double x1, double y1, double z1, double x2, double y2, double z2, double kappa) {

        // Finding the distance between consecutive points
        double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
        return distance * kappa;
        //cout << "The distance between (" << points[i+1][0] <<"," << points[i+1][1] "," << points[i+1][2] << ")";
        //cout << " and (" << points[i][0] <<"," << points[i][1] "," << points[i][2] << ")" << " is " << distance << endl;
}

void ScCalc::crateSpineFit(int spineNum, unsigned order)
{
    double (*spine)[3];
    double (*fit)[3];
    unsigned spLength;
    double *maXanX;
    if (spineNum == 1) {
        spine = spine1;
        spLength = spine1Length;
        fit1 = new double[spLength][3];
        fit = fit1;
        angles1 = new double[spLength - 3];
        maXanX = angles1;
    }
    double *dx = new double[spLength];
    double *dy = new double[spLength];
    double *dz = new double[spLength];

    double *xStore = new double[order];
    double *yStore = new double[order];
    cout << "The spine length is " << spLength;
    for (unsigned i = 0; i < spLength; ++i) 
    {
        dx[i] = spine[i][0];
        dy[i] = spine[i][1];
        dz[i] = spine[i][2];
    }

    polynomialfit(spLength, order, dz, dx, xStore);
    polynomialfit(spLength, order, dz, dy, yStore);
    /*xStore[0] = 0;
    xStore[1] = 0;
    xStore[2] = 0;
    xStore[3] = 0;
    xStore[4] = 0;
    xStore[5] = 0;
    xStore[6] = 0;*/
    cout << endl;
    cout << "Spine" << spineNum << " polynomial curve:" << endl;
    cout << endl;

    for (int i = 0; i < order; ++i)
    {
        cout << i << " order term ";
        cout << "X: " << xStore[i];
        cout << " Y: " << yStore[i] << endl;
    }

    for (unsigned i = 0; i < spLength; ++i) 
    {
        double sumX = 0;
        double sumY = 0;

        for (unsigned j = 0; j < order; ++j) 
        {
            sumX += xStore[j] * pow(dz[i], j);
            sumY += yStore[j] * pow(dz[i], j);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;

        fit[i][0] = sumX;
        fit[i][1] = sumY;
        fit[i][2] = dz[i];
    }

    double *anX = new double[spLength];
    double *anY = new double[spLength];
    

    anglesX(fit, anX, spLength);
    anglesY(fit, anY, spLength);
    std::cout << "Geometric curvature for " << std::endl;
    getGeometricCurvatureMod(fit, spLength, order, xStore, yStore, dz, 0);
    getCurvatureAngleMod(fit, spLength, order, xStore, yStore, dz, 0);
    analyze3DAngles(getMax3DanglesMod(fit, maXanX, spLength), 0);
}

void ScCalc::crateSpineFitX(int spineNum, unsigned order)
{
    std::cout<<"Analyzing angles in X" << std::endl;
    double (*spine)[3];
    double (*fit)[3];
    unsigned spLength;
    double *maXanX;
    if (spineNum == 1) {
        spine = spine1;
        spLength = spine1Length;
        fit1 = new double[spLength][3];
        fit = fit1;
        angles1 = new double[spLength - 3];
        maXanX = angles1;
    }
    double *dx = new double[spLength];
    double *dy = new double[spLength];
    double *dz = new double[spLength];

    double *xStore = new double[order];
    double *yStore = new double[order];
    cout << "The spine length is " << spLength;
    for (unsigned i = 0; i < spLength; ++i) 
    {
        dx[i] = spine[i][0];
        dy[i] = spine[i][1];
        dz[i] = spine[i][2];
    }

    polynomialfit(spLength, order, dz, dx, xStore);

    yStore[0] = 0;
    yStore[1] = 0;
    yStore[2] = 0;
    yStore[3] = 0;
    yStore[4] = 0;
    yStore[5] = 0;
    yStore[6] = 0;
    cout << endl;
    cout << "Spine" << spineNum << " polynomial curve:" << endl;
    cout << endl;

    for (int i = 0; i < order; ++i)
    {
        cout << i << " order term ";
        cout << "X: " << xStore[i];
        cout << " Y: " << yStore[i] << endl;
    }

    for (unsigned i = 0; i < spLength; ++i) 
    {
        double sumX = 0;
        double sumY = 0;

        for (unsigned j = 0; j < order; ++j) 
        {
            sumX += xStore[j] * pow(dz[i], j);
            sumY += yStore[j] * pow(dz[i], j);
        }

        cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;

        fit[i][0] = sumX;
        fit[i][1] = sumY;
        fit[i][2] = dz[i];
    }

    double *anX = new double[spLength];
    double *anY = new double[spLength];
    

    //anglesX(fit, anX, spLength);
    //anglesY(fit, anY, spLength);
    std::cout << "Geometric curvature for X " << std::endl;
    getGeometricCurvatureMod(fit, spLength, order, xStore, yStore, dz, 1);
    getCurvatureAngleMod(fit, spLength, order, xStore, yStore, dz, 1);
    computeGeometricCurvature2D(fit, spLength, order, xStore, dz);
    compX = getMax3DanglesMod(fit, maXanX, spLength);
    cout << "Printing the maxAnX for X before components" << endl;
    for(unsigned i; i < spLength; i++) {
        cout  << maXanX[i] << endl;
        maxAnglesX.push_back(maXanX[i]);
    }

    analyze3DAngles(compX, 1);
}

void ScCalc::crateSpineFitY(int spineNum, unsigned order)
{
    std::cout<<"Analyzing angles in Y" << std::endl;
    double (*spine)[3];
    double (*fit)[3];
    unsigned spLength;
    double *maXanX;
    if (spineNum == 2) {
        spine = spine1;
        spLength = spine1Length;
        fit2 = new double[spLength][3];
        fit = fit2;
        angles2 = new double[spLength - 3];
        maXanX = angles2;
    }
    double *dx = new double[spLength];
    double *dy = new double[spLength];
    double *dz = new double[spLength];

    double *xStore = new double[order];
    double *yStore = new double[order];
    cout << "The spine length is " << spLength;
    for (unsigned i = 0; i < spLength; ++i) 
    {
        dx[i] = spine[i][0];
        dy[i] = spine[i][1];
        dz[i] = spine[i][2];
    }

    //polynomialfit(spLength, order, dz, dx, xStore);
    polynomialfit(spLength, order, dz, dy, yStore);
    xStore[0] = 0;
    xStore[1] = 0;
    xStore[2] = 0;
    xStore[3] = 0;
    xStore[4] = 0;
    xStore[5] = 0;
    xStore[6] = 0;
    cout << endl;
    cout << "Spine" << spineNum << " polynomial curve:" << endl;
    cout << endl;

    for (int i = 0; i < order; ++i)
    {
        cout << i << " order term ";
        cout << "X: " << xStore[i];
        cout << " Y: " << yStore[i] << endl;
    }

    for (unsigned i = 0; i < spLength; ++i) 
    {
        double sumX = 0;
        double sumY = 0;

        for (unsigned j = 0; j < order; ++j) 
        {
            sumX += xStore[j] * pow(dz[i], j);
            sumY += yStore[j] * pow(dz[i], j);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;

        fit[i][0] = sumX;
        fit[i][1] = sumY;
        fit[i][2] = dz[i];
    }

    double *anX = new double[spLength];
    double *anY = new double[spLength];
    

    //anglesX(fit, anX, spLength);
    //anglesY(fit, anY, spLength);
    std::cout << "Geometric curvature for Y" << std::endl;
    getGeometricCurvatureMod(fit, spLength, order, xStore, yStore, dz, 2);
    computeGeometricCurvature2D(fit, spLength, order, yStore, dz);
    getCurvatureAngleMod(fit, spLength, order, xStore, yStore, dz, 2);
    compY = getMax3DanglesMod(fit, maXanX, spLength);
    cout << "Printing the maxAnX for Y before components" << endl;
    for(unsigned i; i < spLength; i++) {
        cout << maXanX[i] << endl;
        maxAnglesY.push_back(maXanX[i]);
    }
    analyze3DAngles(compY, 2);
}

bool ScCalc::polynomialfit(int obs, int degree, 
           double *dx, double *dy, double *store) /* n, p */
{
  gsl_multifit_linear_workspace *ws;
  gsl_matrix *cov, *X;
  gsl_vector *y, *c;
  double chisq;
 
  int i, j;
 
  X = gsl_matrix_alloc(obs, degree);
  y = gsl_vector_alloc(obs);
  c = gsl_vector_alloc(degree);
  cov = gsl_matrix_alloc(degree, degree);
 
  for(i=0; i < obs; i++) {
    for(j=0; j < degree; j++) {
      gsl_matrix_set(X, i, j, pow(dx[i], j));
    }
    gsl_vector_set(y, i, dy[i]);
  }
 
  ws = gsl_multifit_linear_alloc(obs, degree);
  gsl_multifit_linear(X, y, c, cov, &chisq, ws);
 
  /* store result ... */
  for(i=0; i < degree; i++)
  {
    store[i] = gsl_vector_get(c, i);
  }
 
  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  return true; /* we do not "analyse" the result (cov matrix mainly)
          to know if the fit is "good" */
}

void ScCalc::getMax3Dangles(double points[][3], double angles[], int npoints)
{
    for (int i = 1; i < npoints - 2; ++i) 
    {
        double *vec1 = new double[3];
        for (int k = 0; k < 3; ++k) vec1[k] = points[i][k] - points[i-1][k];
        angles[i - 1] = 0;

        for (int j = i + 2; j < npoints; ++j)
        {

            double *vec2 = new double[3];
            for (int k = 0; k < 3; ++k) vec2[k] = points[j][k] - points[j - 1][k];

            double dotProduct = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
            double mag1 = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]);
            double mag2 = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]);

            double newAngle = acos(dotProduct / (mag1 * mag2)) * 180/3.1415926358979323;
            if (newAngle > angles[i - 1]) {
                //cout << "Max so far: " << i << " and " << j << " with angle " << newAngle << endl;
              angles[i - 1] = newAngle;
            }
            else {
                //cout << "Broke: " << i << " and " << j << " with angle " << newAngle << endl;
                break;
            }
        }
    }
}


int getClosestIntersectionPoint(double points[][3], int npoints, double* p1, double* p2, double* p3, double* p4) {

    double x1 = 0;
    double x2 = 0;
    double x3 = 0;
    double x4 = 0;
    double y1 = 0;
    double y2 = 0;
    double y3 = 0;
    double y4 = 0;
    if (p1[0] == 0 &&  p2[0] == 0 && p2[0] == p3[0]) {
         x1 = p1[1];
         x2 = p2[1];
         x3 = p3[1];
         x4 = p4[1];
         y1 = p1[2];
         y2 = p2[2];
         y3 = p3[2];
         y4 = p4[2];
    }
    else {
         x1 = p1[0];
         x2 = p2[0];
         x3 = p3[0];
         x4 = p4[0];
         y1 = p1[2];
         y2 = p2[2];
         y3 = p3[2];
         y4 = p4[2];
    }
    cout << x1 << "  " << x2 << "  " << y1 << endl;
    double px = (((x1*y2 - y1*x2)*(x3-x4)) - ((x1-x2)*(x3*y4 - y3*y4)))/((x1-x2)*(y3-y4) - (y1-y2)*(x3-x4));
    double py = (((x1*y2 - y1*x2)*(y3-y4)) - ((y1-y2)* (x3*y4 - y3*x4))) /
    ((x1-x2)*(y3-y4) - (y1-y2)*(x3-x4));
    // cout << "Point of intersection " << px << " and " << py << endl;

    // Find the closest one the intersection point, use only the z for now.
    int minIndex = 0;
    double minDistance = abs(py - points[0][2]);
    for (unsigned i = 1; i < npoints; i ++ ) {
        if (abs(py - points[i][2]) < minDistance) {
            minDistance = abs( py - points[i][2]);
            minIndex = i;
        }
    }
    return minIndex;
}


// added this
// Finding the maxangles within the section. 
// Compare with 3D angles and the 2D results with planar projections
vector<SpineComponent> ScCalc::getMax3DanglesMod(double points[][3], double angles[], int npoints)
{
  map<int, map<int, double>> zAngleMap;
  bool detectedInflection = false;
  vector<SpineComponent> components;

    for (int i = 1; i < npoints - 2; ++i) 
    {
        double *vec1 = new double[3];
        for (int k = 0; k < 3; ++k) vec1[k] = points[i][k] - points[i-1][k];
        angles[i - 1] = 0;
        int j = 0;
        double newAngle = 0.0;
        for (j = i + 2; j < npoints; ++j)
        {

            double *vec2 = new double[3];
            for (int k = 0; k < 3; ++k) vec2[k] = points[j][k] - points[j - 1][k];

            double dotProduct = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
            double mag1 = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]);
            double mag2 = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]);

            newAngle = acos(dotProduct / (mag1 * mag2)) * 180/3.1415926358979323;
            if (newAngle > angles[i - 1]) {

               cout << "Max so far: " << i << " and " << j << " with angle " << newAngle << endl;
              angles[i - 1] = newAngle;
              zAngleMap[i][j] = newAngle;
            }
            else {
                cout << "Broke: " << i << " and " << j << " with angle " << newAngle << endl;
                detectedInflection = true;
                break;
            }
        }
        if (!detectedInflection) {
          zAngleMap[i][j-1] = newAngle;
        }
        detectedInflection = false;
    }

    // Print out the map to figure out what the extremes might be
    for(auto it = zAngleMap.begin(); it != zAngleMap.end(); ++it) {
      for(auto it1 = it->second.begin(); it1 != it->second.cend(); ++it1) {
        std::cout << it->first-1 << " and " << it->first << "  " << it1->first-1 << " and "  << it1->first << " \tAngle: " << it1->second << "\n";
        SpineComponent sp;
        sp.index1 = it->first;
        sp.index2 = it1->first;
        sp.angle = it1->second;
        // find the intersection point of both points and match it to the closest vertebra centroid.
        sp.apex = getClosestIntersectionPoint(points, npoints, points[it->first], points[it->first - 1], points[it1->first], points[it1->first - 1]);
        cout << sp.apex << endl;
        components.push_back(sp);
      }
    }
    return components;
}


bool compareByAngle(const SpineComponent &a, const SpineComponent &b) {
    return a.angle > b.angle;
}
// added this
void ScCalc::analyze3DAngles(std::vector<SpineComponent> comp, int type) {
  // find the largest angle - record the endpoints.
  // sort the angles based on
  std::sort(comp.begin(), comp.end(), compareByAngle);
  for(int i = 0; i < comp.size(); i++ ) {
    std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
  }

  // filter the results to print sectional results.
  std::cout << "Printing out the final angles " << endl;
  int deviation = 1;
  std::vector<SpineComponent> results;
  for(int i = 0; i < comp.size(); i++) {
    // Check if the span is already covered by another curve
    bool includeInResults = true;
    for (int j = 0; j < results.size(); j++) {
      if (comp[i].index1 >= results[j].index1-deviation && comp[i].index2 <= results[j].index2+deviation) {
        includeInResults = false;
        break;
      }
    }
    if (includeInResults == true) {
      results.push_back(comp[i]);
      std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << " Apex at " << comp[i].apex+1;
      //std::cout << vertebraePoints.size() << endl;
      std::cout << " i.e. between : " << vertebraePoints[comp[i].index1-1] << " and " << vertebraePoints[comp[i].index2-1] << " Apex at " << vertebraePoints[comp[i].apex] << endl;
    }
  }

  // Do the lenke classification here.
  doLenkeClassificationReloaded(comp, results, type);
  doSimpleClassification(comp, results);
}

void ScCalc::doLenkeClassificationReloaded(std::vector<SpineComponent> comp, std::vector<SpineComponent> results, int type) {
    if (type == 0) {
        doLenkeClassification(comp, results);
    }
    else {
        if (type == 1) {
            std::string primaryCurve = "";
            // Find some important indicies:
            int t5Index = 0;
            int t12Index = 0;
            int s1ClosestIndex = 0;
            bool s1Assigned = false;
            for(int i = 0; i < vertebraePoints.size(); i++) {
                if (vertebraePoints[i].compare("T5_center") == 0)
                    t5Index = i;
                else if (vertebraePoints[i].compare("T12_center") == 0) 
                    t12Index = i;
                else if (vertebraePoints[i].compare("S1_center") == 0) {
                    s1ClosestIndex = i;
                    s1Assigned = true;
                }
            }
            int offset = 1;
            if (!s1Assigned) {
                if (t5Index < t12Index) {
                    s1ClosestIndex = vertebraePoints.size()-1;
                    offset =-1;
                }
                else
                    s1ClosestIndex = 0;
            }

            // Step 1: Find the apex (CORONAL)
            // Find the average X of all vertibrae, then find one with the greatest deviation from the average.
            double sum = 0;
            double avg = 0;
            for (int i = 0; i <  vertebraePoints.size(); i++) {
                sum += fit1[i][0];
            }
            avg = sum / vertebraePoints.size();

            std::string apexVert = "";
            int apexVertIndex = 0;
            for (int i = 0; i <  vertebraePoints.size(); i++) {
                if (vertebraePoints[i] == "T3_center") {
                    apexVertIndex = i;
                    break;
                }
            }

            for (int i = apexVertIndex; i <  vertebraePoints.size(); i++) {
                std::cout << "Distance is " << abs(avg - fit1[i][0]) << endl;
                if (abs(avg - fit1[i][0]) > abs(avg - fit1[apexVertIndex][0])) {
                    apexVertIndex = i;
                }
            }
            apexVert = vertebraePoints[apexVertIndex];
            std::cout << "The apex vert is " << apexVert << endl;


            curveType = "Undefined";
            for (int i = 0; i < results.size(); i++) {
                if (i == 0 && vertebraePoints[results[i].index1].substr(0, 1) == "T" && vertebraePoints[results[i].index2].substr(0, 1) == "T") {
                    // This type 1 to 4 with MT as major curve
                    // Include criteria to check 1 to 4
                    std::cout << "Type 1,2,3,4";
                    curveType = "1";
                }

                if ((i == 0 && vertebraePoints[results[i].index1].substr(0, 1) == "T" && vertebraePoints[results[i].index2].substr(0, 1) == "L") || 
                    (i == 0 && vertebraePoints[results[i].index1].substr(0, 1) == "L" && vertebraePoints[results[i].index2].substr(0, 1) == "T") || 
                    (i == 0 && vertebraePoints[results[i].index1].substr(0, 1) == "L" && vertebraePoints[results[i].index2].substr(0, 1) == "L")) {
                    // Type 4 or 5 or 6 
                    std::cout<<"Type 4,5,6 ";
                    curveType = "6";
                    apexVert = vertebraePoints[apexVertIndex];
                    std::cout << "The apex vert is " << apexVert << endl;
                    // Classify to the appropriate curve type
                    if (apexVert.substr(0,2) == "T6" || apexVert.substr(0,2) == "T7" || apexVert.substr(0,2) == "T8" || 
                        apexVert.substr(0,2) == "T9" || apexVert.substr(0,3) == "T10" || apexVert.substr(0,3) == "T11") {
                        curveType = "1";
                    }
                }
            }

            // Step 2: Find lumbar modifier (CORONAL)
            // Find the CSVL
            sum = 0;
            avg = 0;
            int lumbarCount = 0;
            for (int i = 0; i <  vertebraePoints.size(); i++) {
                if (vertebraePoints[i].substr(0, 1) == "L"){
                    sum += fit1[i][0];
                    lumbarCount += 1;
                }
            }
            avg = sum / lumbarCount;

            std::cout<< fit1[s1ClosestIndex][0] << " is the xcoordinate of s1" << std::endl;
            int maxLumbarIndex = 0;
            // find the lumbar vertebrae that is furthest from the CSVL
            for(int i = 0; i < vertebraePoints.size(); i++) {
                if (vertebraePoints[i].substr(0, 1) == "L"){
                    if (abs(avg - fit1[i][0]) > abs(avg - fit1[maxLumbarIndex][0])) {
                        maxLumbarIndex = i;
                    }
                }
            }
            int vertSize = 0;
            // Assume the size of the vertebrae is the distance between two consecutive vertebrae
            vertSize = abs(fit1[s1ClosestIndex][0] - fit1[s1ClosestIndex + offset][0]);
            // Assign the appropriate lumbar modifier
            if (fit1[maxLumbarIndex][0] - fit1[s1ClosestIndex][0] > (vertSize/2)) {
                lumbarModifier = " Modifier C";
            }
            else if (fit1[maxLumbarIndex][0] - fit1[s1ClosestIndex][0] < (vertSize/2)) {
                lumbarModifier = " Modifier A ";
            }
            else {
                lumbarModifier = " Modifier B";
            }
        }
        else if (type == 2) {
            // Find some important indicies:
            int t2Index = 0;
            int t5Index = 0;
            int t10Index = 0;
            int t12Index = 0;
            int l2Index = 0;
            int s1ClosestIndex = 0;
            bool s1Assigned = false;
            for(int i = 0; i < vertebraePoints.size(); i++) {
                if (vertebraePoints[i].compare("T5_center") == 0)
                    t5Index = i;
                else if (vertebraePoints[i].compare("T12_center") == 0) 
                    t12Index = i;
                else if (vertebraePoints[i].compare("S1_center") == 0) {
                    s1ClosestIndex = i;
                    s1Assigned = true;
                }
            }

          // Step 1
          bool isPT = false;
          bool isMT = false;
          bool isTL = false;
          if (curveType == "1") {
            // Definitely Type 1 to 4
            isMT = true;
            // check for minor curves
            for(int i = 0; i < comp.size(); i++ ) {
              if ((comp[i].index1 == t2Index && comp[i].index2 == t5Index) || 
               (comp[i].index1 == t5Index && comp[i].index2 == t2Index)) {
                  std::cout << "The angle is " << comp[i].angle << " between t2 and t5" << std::endl;
                  std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
                  if (comp[i].angle >= 20) {
                      isPT = true;
                  }   
              }

              if ((comp[i].index1 == t10Index && comp[i].index2 == l2Index) || 
               (comp[i].index1 == l2Index && comp[i].index2 == t10Index)) {
                  std::cout << "The angle is " << comp[i].angle << " between t10 and l2" << std::endl;
                  std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
                  if (comp[i].angle >= 20) {
                      isTL = true;
                  }
              }
            }
            std::cout << "The CURVE TYPE IS ";
            if (isPT && isTL) {
              curveType =  "Type 4";
            }
            else if (isPT && !isTL) {
              curveType = "Type 2";
            }
            else if (!isPT && isTL) {
              curveType = "Type 3";
            }
            else if (!isPT && !isTL) {
              curveType = "Type 1";
            }
          }
          else if (curveType == "6") {
            // Definitely Type 4 to 6
            isTL = true;
            // check for minor curves
            for(int i = 0; i < comp.size(); i++ ) {
              if ((comp[i].index1 == t2Index && comp[i].index2 == t5Index) || 
               (comp[i].index1 == t5Index && comp[i].index2 == t2Index)) {
                  std::cout << "The angle is " << comp[i].angle << " between t2 and t5" << std::endl;
                  std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
                  if (comp[i].angle >= 20) {
                      isPT = true;
                  }   
              }

              if ((comp[i].index1 == t10Index && comp[i].index2 == l2Index) || 
               (comp[i].index1 == l2Index && comp[i].index2 == t10Index)) {
                  std::cout << "The angle is " << comp[i].angle << " between t10 and l2" << std::endl;
                  std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
                  if (comp[i].angle >= 20) {
                      isMT = true;
                  }
              }
            std::cout << "The CURVE TYPE IS ";
            if (isPT && isMT) {
              curveType = "Type 4";
            }
            else if (!isPT && isMT) {
              curveType = "Type 6";
            }
            else if (!isPT && !isMT) {
              curveType = "Type 5";
            }

          }
        }

                  // Step 1.5: Find the non-structural curves
          for(int i = 0; i < comp.size(); i++ ) {
            if ((comp[i].index1 == t2Index && comp[i].index2 == t5Index) || 
             (comp[i].index1 == t5Index && comp[i].index2 == t2Index)) {
                std::cout << "The angle is " << comp[i].angle << " between t2 and t5" << std::endl;
                std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
                if (comp[i].angle >= 20) {
                    sagittalModifier = " - (hypokhyphosis)";
                }
                else if (comp[i].angle >= 10 && comp[i].angle <= 40) {
                    sagittalModifier = " N (normal)";
                }
                else {
                    sagittalModifier = " + (hyperkhyphosis)";
                }   
            }
          }

            // Step 3: Find sagittal modifier (SAGITTAL)
            // Measure T5 to T12
          for(int i = 0; i < comp.size(); i++ ) {
            if ((comp[i].index1 == t5Index && comp[i].index2 == t12Index) || 
             (comp[i].index2 == t5Index && comp[i].index1 == t12Index)) {
                std::cout << "The angle is " << comp[i].angle << " between t5 and t12" << std::endl;
                std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
                if (comp[i].angle < 10) {
                    sagittalModifier = " - (hypokhyphosis)";
                }
                else if (comp[i].angle >= 10 && comp[i].angle <= 40) {
                    sagittalModifier = " N (normal)";
                }
                else {
                    sagittalModifier = " + (hyperkhyphosis)";
                }   
            }
          }

          std::string classification = curveType + lumbarModifier + sagittalModifier;
          std::cout << "The final lenke classification is " << classification << endl;
      }
    }
}

void ScCalc::doLenkeClassification(std::vector<SpineComponent> comp, std::vector<SpineComponent> results) {
    std::string primaryCurve = "";
    std::string lumbarModifier = "";
    std::string sagittalModifier = "";

    // Find some important indicies:
    int t5Index = 0;
    int t12Index = 0;
    int s1ClosestIndex = 0;
    bool s1Assigned = false;
    for(int i = 0; i < vertebraePoints.size(); i++) {
        if (vertebraePoints[i].compare("T5_center") == 0)
            t5Index = i;
        else if (vertebraePoints[i].compare("T12_center") == 0) 
            t12Index = i;
        else if (vertebraePoints[i].compare("S1_center") == 0) {
            s1ClosestIndex = i;
            s1Assigned = true;
        }
    }
    int offset = 1;
    if (!s1Assigned) {
        if (t5Index < t12Index) {
            s1ClosestIndex = vertebraePoints.size()-1;
            offset =-1;
        }
        else
            s1ClosestIndex = 0;
    }

    // Step 1: Find the apex (CORONAL)
    // Find the average X of all vertibrae, then find one with the greatest deviation from the average.
    double sum = 0;
    double avg = 0;
    for (int i = 0; i <  vertebraePoints.size(); i++) {
        sum += fit1[i][0];
    }
    avg = sum / vertebraePoints.size();

    std::string apexVert = "";
    int apexVertIndex = 0;
    for (int i = 0; i <  vertebraePoints.size(); i++) {
        if (vertebraePoints[i] == "T3_center") {
            apexVertIndex = i;
            break;
        }
    }

    for (int i = apexVertIndex; i <  vertebraePoints.size(); i++) {
        std::cout << "Distance is " << abs(avg - fit1[i][0]) << endl;
        if (abs(avg - fit1[i][0]) > abs(avg - fit1[apexVertIndex][0])) {
            apexVertIndex = i;
        }
    }
    apexVert = vertebraePoints[apexVertIndex];
    std::cout << "The apex vert is " << apexVert << endl;
    std::string curveType = "Undefined";
    // Classify to the appropriate curve type
    if (apexVert.substr(0,2) == "T3" || apexVert.substr(0,2) == "T4" || apexVert.substr(0,2) == "T5") {
        curveType = "Proximal Thoracic (Type 1)";
    } else if (apexVert.substr(0,2) == "T6" || apexVert.substr(0,2) == "T7" || apexVert.substr(0,2) == "T8" || 
        apexVert.substr(0,2) == "T9" || apexVert.substr(0,3) == "T10" || apexVert.substr(0,3) == "T11") {

        curveType = "Main Thoracic (Type 2/3/4)";
    } else if ( apexVert.substr(0,2) == "T12" || apexVert.substr(0,1) == "L") {
        curveType = "Thoracolumbar/Lumbar (TL) (Type 5/6)";
    }
    cout << "Curve type is " << curveType <<endl;

    // Step 2: Find lumbar modifier (CORONAL)
    // Find the CSVL
    sum = 0;
    avg = 0;
    int lumbarCount = 0;
    for (int i = 0; i <  vertebraePoints.size(); i++) {
        if (vertebraePoints[i].substr(0, 1) == "L"){
            sum += fit1[i][0];
            lumbarCount += 1;
        }
    }
    avg = sum / lumbarCount;

    std::cout<< fit1[s1ClosestIndex][0] << " is the xcoordinate of s1" << std::endl;
    int maxLumbarIndex = 0;
    // find the lumbar vertebrae that is furthest from the CSVL
    for(int i = 0; i < vertebraePoints.size(); i++) {
        if (vertebraePoints[i].substr(0, 1) == "L"){
            if (abs(avg - fit1[i][0]) > abs(avg - fit1[maxLumbarIndex][0])) {
                maxLumbarIndex = i;
            }
        }
    }
    int vertSize = 0;
    // Assume the size of the vertebrae is the distance between two consecutive vertebrae
    vertSize = abs(fit1[s1ClosestIndex][0] - fit1[s1ClosestIndex + offset][0]);
    // Assign the appropriate lumbar modifier
    if (fit1[maxLumbarIndex][0] - fit1[s1ClosestIndex][0] > (vertSize/2)) {
        lumbarModifier = " Modifier C";
    }
    else if (fit1[maxLumbarIndex][0] - fit1[s1ClosestIndex][0] < (vertSize/2)) {
        lumbarModifier = " Modifier A ";
    }
    else {
        lumbarModifier = " Modifier B";
    }

    // Step 3: Find sagittal modifier (SAGITTAL)
    // Measure T5 to T12
  for(int i = 0; i < comp.size(); i++ ) {
    if ((comp[i].index1 == t5Index && comp[i].index2 == t12Index) || 
     (comp[i].index2 == t5Index && comp[i].index1 == t12Index)) {
        std::cout << "The angle is " << comp[i].angle << " between t5 and t12" << std::endl;
        std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
        if (comp[i].angle < 10) {
            sagittalModifier = " - (hypokhyphosis)";
        }
        else if (comp[i].angle >= 10 && comp[i].angle <= 40) {
            sagittalModifier = " N (normal)";
        }
        else {
            sagittalModifier = " + (hyperkhyphosis)";
        }   
    }
  }
  std::string classification = curveType + lumbarModifier + sagittalModifier;
  std::cout << "The final lenke classification is " << classification << endl;
}


void ScCalc::doSimpleClassification(std::vector<SpineComponent> comp, std::vector<SpineComponent> results) {
    // Classify with Thoracic scoliosis or lumbar scoliosis

    // Measure TK and LL
    int t4Index = 0;
    int t12Index = 0;
    int s1ClosestIndex = 0;
    int l1Index = 0;
    bool s1Assigned = false;
    for(int i = 0; i < vertebraePoints.size(); i++) {
        if (vertebraePoints[i].compare("T4_center") == 0)
            t4Index = i;
        else if (vertebraePoints[i].compare("T12_center") == 0) 
            t12Index = i;
        else if (vertebraePoints[i].compare("S1_center") == 0) {
            s1ClosestIndex = i;
            s1Assigned = true;
        }
        else if (vertebraePoints[i].compare("L1_center") == 0) {
            l1Index = i;
        }
    }

   for(int i = 0; i < comp.size(); i++ ) {
    if ((comp[i].index1 == t4Index && comp[i].index2 == t12Index) || 
     (comp[i].index2 == t4Index && comp[i].index1 == t12Index)) {
        std::cout << "TK: " << comp[i].angle << " between t4 and t12" << std::endl;
        std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
    }
    if ((comp[i].index1 == s1ClosestIndex && comp[i].index2 == l1Index) || 
     (comp[i].index2 == l1Index && comp[i].index1 == s1ClosestIndex)) {
        std::cout << "LL: " << comp[i].angle << " between L1 and S1" << std::endl;
        std::cout << i + 1 << "." << " Angle: " << comp[i].angle << " between: " << comp[i].index1 << " and " << comp[i].index2 << endl;
    }
  }   

    // Measure LL


}

void ScCalc::maXanglesX(double points[][3], double angles[], int npoints){
    for (int j = 0; j < npoints; j++)
    {
        double previousAngle = 360;

        for(int i = j; i < npoints; i++)
        {
            int last = (j - 1 + npoints) % npoints;
            int next = (i + 1) % npoints;
            double x1 = points[i][0] - points[last][0];
            double z1 = points[i][2] - points[last][2];
            double x2 = points[next][0] - points[i][0];
            double z2 = points[next][2] - points[i][2];
            double theta1 = atan2(z1, x1)*180/3.1415926358979323;
            double theta2 = atan2(z2, x2)*180/3.1415926358979323;
            double thisAngle = (180 + theta1 - theta2 + 360);
            while(thisAngle>360)thisAngle-=360;

            if (previousAngle < thisAngle)
                break;
            else
                previousAngle = thisAngle;
        } 
        angles[j] = previousAngle;
    }
}

void ScCalc::anglesX(double points[][3], double angles[], int npoints){
for(int i = 0; i < npoints; i++){
    int last = (i - 1 + npoints) % npoints;
    int next = (i + 1) % npoints;
    double x1 = points[i][0] - points[last][0];
    double z1 = points[i][2] - points[last][2];
    double x2 = points[next][0] - points[i][0];
    double z2 = points[next][2] - points[i][2];
    double theta1 = atan2(z1, x1)*180/3.1415926358979323;
    double theta2 = atan2(z2, x2)*180/3.1415926358979323;
    angles[i] = (180 + theta1 - theta2 + 360);
    std::cout<<"AnglesX= "<<angles[i] << " with " << i << std::endl;
    while(angles[i]>360)angles[i]-=360;
} }

void ScCalc::anglesY(double points[][3], double angles[], int npoints){
for(int i = 0; i < npoints; i++){
    int last = (i - 1 + npoints) % npoints;
    int next = (i + 1) % npoints;
    double y1 = points[i][1] - points[last][1];
    double z1 = points[i][2] - points[last][2];
    double y2 = points[next][1] - points[i][1];
    double z2 = points[next][2] - points[i][2];
    double theta1 = atan2(z1, y1)*180/3.1415926358979323;
    double theta2 = atan2(z2, y2)*180/3.1415926358979323;
    angles[i] = (180 + theta1 - theta2 + 360);
    while(angles[i]>360)angles[i]-=360;
} }

// added this:
vector<vector<int>> ScCalc::loadAnnotationData(char* fileName)
{
    FILE *fp;
    vector<vector<int> > vec;

    string sFileName = fileName;

    size_t lastdot = sFileName.find_last_of(".");

    if (lastdot != std::string::npos) {
        sFileName = sFileName.substr(0, lastdot);
    }
    sFileName.append(".lml");

    std::cout << "The .lml file name is " << sFileName << endl;
    ifstream fileStream(sFileName);
    if (!fileStream.is_open())
    {
        cout << "Exiting unable to open file " << fileName << endl;
    }
    int counter = 0;
    string line;
    int spacingIndex = 0;
    while(getline(fileStream, line, '\n')) {
        if (counter == 0) {
            counter ++;
            continue;
        }
        stringstream ss(line);
        vector<int> numbers;
        string in_line;
        int col = 0;
        while((getline (ss, in_line, '\t') || getline(ss, in_line, '\n')) && col <=4 ) 
        {
          if (col != 0 && col != 1) {
            cout << in_line << " IS ";
            int i = stoi(in_line, 0);
            i = i/spacing1[spacingIndex%3];
            numbers.push_back(i);;
            spacingIndex++;
          }
          col ++;
        }
        if (numbers.size() > 0)
          vec.push_back(numbers);
    }

    return vec;
}

void ScCalc::loadColorMap(char* fileName) {
    FILE *fp;
    vector<vector<int> > vec;

    string sFileName = fileName;

    size_t lastdot = sFileName.find_last_of(".");

    if (lastdot != std::string::npos) {
        sFileName = sFileName.substr(0, lastdot);
    }
    sFileName.append(".lml");

    std::cout << "The .lml file name is " << sFileName << endl;
    ifstream fileStream(sFileName);
    if (!fileStream.is_open())
    {
        cout << "Exiting unable to open file " << fileName << endl;
    }
    int counter = 0;
    string line;
    int spacingIndex = 0;
    while(getline(fileStream, line, '\n')) {
        if (counter == 0) {
            counter ++;
            continue;
        }
        stringstream ss(line);
        vector<int> numbers;
        string in_line;
        int col = 0;
        while((getline (ss, in_line, '\t')) )
        {
            if (col == 1) {
                vertebraePoints.push_back(in_line);
            }
          col ++;
        }
    }
}
// //____________________________________________________________________
// void ScCalc::makeData(Double_t* x, Double_t& d, Double_t& e)
// {
//   // Make data points
//   Double_t upp[5] = { 10, 10, 10, 10,  1 };
//   Double_t low[5] = {  0,  0,  0,  0, .1 };
//   for (int i = 0; i < 2; i++)
//     x[i] = (upp[i] - low[i]) * gRandom->Rndm() + low[i];

//   d = x[0] * TMath::Sqrt(x[1] * x[1]);

//   e = gRandom->Gaus(upp[4],low[4]);
// }

// //____________________________________________________________________
// int ScCalc::CompareResults(TMultiDimFit *fit, bool doFit)
// {
//    //Compare results with reference run


//    // the right coefficients (before fit)
//   double GoodCoeffsNoFit[] = {
//   -4.37056,
//   43.1468,
//   13.432,
//   13.4632,
//   13.3964,
//   13.328,
//   13.3016,
//   13.3519,
//   4.49724,
//   4.63876,
//   4.89036,
//   -3.69982,
//   -3.98618,
//   -3.86195,
//   4.36054,
//   -4.02597,
//   4.57037,
//   4.69845,
//   2.83819,
//   -3.48855,
//   -3.97612
// };

//    // the right coefficients (after fit)
//   double GoodCoeffs[] = {
//      -4.399,
//      43.15,
//      13.41,
//      13.49,
//      13.4,
//      13.23,
//      13.34,
//      13.29,
//      4.523,
//      4.659,
//      4.948,
//      -4.026,
//      -4.045,
//      -3.939,
//      4.421,
//      -4.006,
//      4.626,
//      4.378,
//      3.516,
//      -4.111,
//      -3.823,
// };

// // Good Powers
//   int GoodPower[] = {
//   1,  1, 
//   2,  1, 
//   1,  1,  
//   1,  1, 
//   1,  2,  
//   2,  2,  
//   2,  1,  
//   2,  1, 
//   1,  1,  
//   1,  3,  
//   1,  1,  
//   1,  1,  
//   1,  2,  
//   1,  2,  
//   2,  1,  
//   2,  2,  
//   2,  1,  
//   2,  3,  
//   1,  2,  
//   2,  1,  
//   2,  2
// };

//   Int_t nc = fit->GetNCoefficients();
//   Int_t nv = fit->GetNVariables();
//   const Int_t *powers = fit->GetPowers();
//   const Int_t *pindex = fit->GetPowerIndex();
//   if (nc != 21) return 1;
//   const TVectorD *coeffs = fit->GetCoefficients();
//   int k = 0;
//   for (Int_t i=0;i<nc;i++) {
//      if (doFit) {
//         if (!TMath::AreEqualRel((*coeffs)[i],GoodCoeffs[i],1e-3)) return 2;
//      }
//      else {
//         if (TMath::Abs((*coeffs)[i] - GoodCoeffsNoFit[i]) > 5e-5) return 2;
//      }
//      for (Int_t j=0;j<nv;j++) {
//         if (powers[pindex[i]*nv+j] != GoodPower[k]) return 3;
//         k++;
//      }
//   }

//   // now test the result of the generated function
//   gROOT->ProcessLine(".L MDF.C");

//   Double_t refMDF = (doFit) ? 43.95 : 43.98;
//   // this does not work in CLing since the function is not defined
//   //Double_t x[]    = {5,5,5,5};
//   //Double_t rMDF   = MDF(x);
//   //LM:  need to return the address of the result since it is casted to a long (this should not be in a tutorial !)
//   Long_t iret = gROOT->ProcessLine(" Double_t x[] = {5,5,5,5}; double result=MDF(x); &result;");
//   Double_t rMDF = * ( (Double_t*)iret);
//   //printf("%f\n",rMDF);
//   if (TMath::Abs(rMDF -refMDF) > 1e-2) return 4;
//   return 0;
// }

// //____________________________________________________________________
// Int_t ScCalc::multidimfit(bool doFit)
// {

//   cout << "*************************************************" << endl;
//   cout << "*             Multidimensional Fit              *" << endl;
//   cout << "*                                               *" << endl;
//   cout << "* By Christian Holm <cholm@nbi.dk> 14/10/00     *" << endl;
//   cout << "*************************************************" << endl;
//   cout << endl;

//   // Initialize global TRannom object.
//   gRandom = new TRandom();

//   // Open output file
//   TFile* output = new TFile("mdf.root", "RECREATE");

//   // Global data parameters
//   Int_t nVars       = 2;
//   Int_t nData       = 500;
//   Double_t x[2];

//   // make fit object and set parameters on it.
//   TMultiDimFit* fit = new TMultiDimFit(nVars, TMultiDimFit::kMonomials,"v");

//   Int_t mPowers[]   = { 6 , 6};
//   fit->SetMaxPowers(mPowers);
//   fit->SetMaxFunctions(1000);
//   fit->SetMaxStudy(1000);
//   fit->SetMaxTerms(30);
//   fit->SetPowerLimit(1);
//   fit->SetMinAngle(10);
//   fit->SetMaxAngle(10);
//   fit->SetMinRelativeError(.01);

//   // variables to hold the temporary input data
//   Double_t d;
//   Double_t e;

//   // Print out the start parameters
//   fit->Print("p");

//   printf("======================================\n");

//   // Create training sample
//   Int_t i;
//   for (i = 0; i < nData ; i++) {

//     // Make some data
//     makeData(x,d,e);

//     // Add the row to the fit object
//     fit->AddRow(x,d,e);
//   }

//   // Print out the statistics
//   fit->Print("s");

//   // Book histograms
//   fit->MakeHistograms();

//   // Find the parameterization
//   fit->FindParameterization();

//   // Print coefficents
//   fit->Print("rc");

//   // Get the min and max of variables from the training sample, used
//   // for cuts in test sample.
//   Double_t *xMax = new Double_t[nVars];
//   Double_t *xMin = new Double_t[nVars];
//   for (i = 0; i < nVars; i++) {
//     xMax[i] = (*fit->GetMaxVariables())(i);
//     xMin[i] = (*fit->GetMinVariables())(i);
//   }

//   nData = fit->GetNCoefficients() * 100;
//   Int_t j;

//   // Create test sample
//   for (i = 0; i < nData ; i++) {
//     // Make some data
//     makeData(x,d,e);

//     for (j = 0; j < nVars; j++)
//       if (x[j] < xMin[j] || x[j] > xMax[j])
//     break;

//     // If we get through the loop above, all variables are in range
//     if (j == nVars)
//       // Add the row to the fit object
//       fit->AddTestRow(x,d,e);
//     else
//       i--;
//   }
//   //delete gRandom;

//   // Test the parameterizatio and coefficents using the test sample.
//   if (doFit)
//      fit->Fit("M");

//   // Print result
//   fit->Print("fc v");

//   // Write code to file
//   fit->MakeCode();

//   // Write histograms to disk, and close file
//   output->Write();
//   output->Close();
//   delete output;

//   // Compare results with reference run
//   Int_t compare = CompareResults(fit, doFit);
//   if (!compare) {
//      printf("\nmultidimfit ..............................................  OK\n");
//   } else {
//      printf("\nmultidimfit ..............................................  fails case %d\n",compare);
//   }

//   //What we will need:
//   //fit->Eval(x);

//   // We're done
//   delete fit;
//   return compare;
// }
