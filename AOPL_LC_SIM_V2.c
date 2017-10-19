#include <windows.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <random>

#pragma warning ( disable : 4996 ) // functions like strcpy are now deprecated for security reasons

extern "C" {
	int __declspec(dllexport) APIENTRY UserBulkDefinition(double *data);
	int __declspec(dllexport) APIENTRY UserParamNames(char *data);
	int __declspec(dllexport) APIENTRY UserModelInformation(int *settings);
}
void CrossProduct(double x1, double y1, double z1, double x2, double y2, double z2, double *x3, double *y3, double *z3);
void Normalize(double *x, double *y, double *z);

#define PI 3.14159265358979323846

BOOL WINAPI DllMain(HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved) {
	switch (ul_reason_for_call) {
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}

/* the data are stored as follows:

data[0]  = the total number of values in the passed data array
data[1]  = x position
data[2]  = y position
data[3]  = z position
data[4]  = x cosine
data[5]  = y cosine
data[6]  = z cosine
data[7] = wavelength in µm
data[8] = index of refraction

data[9] = Unscattered path length
if data[9] is longer than the randomly generated scatter position, then the DLL should scatter the ray.
if data[9] is shorter than the randomly generated scatter position, then no scattering should occur!

data[10] = 0 initially
if the DLL scatters the ray return 1 in data[10].
if the DLL returns full polarization data return 2 in data[10].

data[11] = millimeters per unit length (1.0 for mm, 25.4 for inches, 10.0 for cm and 1000.0 for meters)

data[12] = relative energy (to be computed by the dll and returned)
data[13] = output phase added to ray
data[14] = mean path value from dialog box
data[15] = angle value from dialog box
data[16] = a random value to use as a seed
data[17] = the glue distance... do not scatter a distance shorter than this

data[20] = incident Ex real
data[21] = incident Ex imaginary
data[22] = incident Ey real
data[23] = incident Ey imaginary
data[24] = incident Ez real
data[25] = incident Ez imaginary

data 35-37 need to be computed if the DLL sets data[10] = 1 or 2
data[35] = output x cosine
data[36] = output y cosine
data[37] = output z cosine

data 40-45 need to be computed if the DLL sets data[10] = 2
data[40] = output Ex real
data[41] = output Ex imaginary
data[42] = output Ey real
data[43] = output Ey imaginary
data[44] = output Ez real
data[45] = output Ez imaginary

data[51] = input parameter 1 from user
data[52] = input parameter 2 from user
etc... up to data[maxdata] where maxdata = int(data[0])

Return 0 if it works; else return -1.

*/

int __declspec(dllexport) APIENTRY UserBulkDefinition(double *data)
{
	double wave, mean, path, t, theta, phi;
	double sx, sy, sz, nx, ny, nz, tx, ty, tz, rx, ry, rz;
	double cosphi, sinphi, costheta, sintheta;
	double absCenter, abs1e, qEff;
	double alpha, alpha_exp;
	double random_number;
	double vabsCenter;
	
	/* We need some way of randomizing the random number generator */
	std::random_device rd;
	std::mt19937 gen((unsigned int)data[16]);
	std::uniform_real_distribution<> dist(0.0, 1.0); //distribution unifomly generates [0,1)

	//data[12] = 1.0; //Leave ray intensity unchanged.

	wave = data[7]; //the wavelength in microns of the input ray
	path = data[9]; //the unscattered path length
	mean = data[14]; //the mean free path at the absorption max. wavelength?

	if (path <= 0.0) return 0;			// we will not scatter
	if (mean <= 0.0) return 0;			// we will not scatter

	absCenter = data[51]; //absorption center in um
	abs1e = data[52]; //absorption 1/e width in um
	qEff = data[53]; //quantum efficiency
	vabsCenter = data[54];

	//Does Scattering occur in this segment?
get_new_random_number:;
	random_number = dist(gen);
	if (random_number == 1.0) goto get_new_random_number;
	t = -mean*log(1.0 - random_number);
	if (t <= data[17]) goto get_new_random_number; // the scatter path must be longer than the glue distance

	if (t > path) return 0; //scattering does not occur

	//Segment will 'scatter,' but will absorption occur?
	random_number = dist(gen);
	alpha_exp = ((wave - absCenter)*(wave - absCenter) - (absCenter - vabsCenter)*(absCenter - vabsCenter)) / (abs1e*abs1e);
	alpha = exp(-alpha_exp);

	if (alpha < random_number) 
	{
		data[12] = 1.0;
		return 0; //absorption does not occur
	}
	
	data[9] = t; //ray will now travel distance t, then scatter.
	data[10] = 1; //tell Zemax ray will bulk scatter.

	random_number = dist(gen);
	if (random_number < (1 - qEff)) //emission does not occur
	{
		data[12] = 0.0; //terminate this ray
		return 0;
	}

	//pick new theta
get_new_random_number2:;
	random_number = dist(gen);
	if (random_number == 0.0) goto get_new_random_number2;
	phi = 2.0*PI*random_number;

	//pick new phi
get_new_random_number3:;
	random_number = dist(gen);
	if (random_number == 0.0) goto get_new_random_number3;
	theta = acos(2.0*random_number - 1.0);

	//construct a ray that makes an angle theta with respect to the current ray
	rx = data[4]; //direction cosines (a.k.a. specular cosines)
	ry = data[5];
	rz = data[6];

	//determine vector, n, which is normal to r
	if (fabs(rz) > 0.9)
	{
		sx = 1.0;
		sy = 0.0;
		sz = 0.0;
	}
	else
	{
		sx = 0.0;
		sy = 0.0;
		sz = 1.0;
	}

	CrossProduct(sx, sy, sz, rx, ry, rz, &nx, &ny, &nz);
	Normalize(&nx, &ny, &nz);

	CrossProduct(nx, ny, nz, rx, ry, rz, &tx, &ty, &tz);
	Normalize(&tx, &ty, &tz);

	//find a random vector s that lies in the plane of n and t
	cosphi = cos(phi);
	sinphi = sin(phi);
	sx = cosphi*nx + sinphi*tx;
	sy = cosphi*ny + sinphi*ty;
	sz = cosphi*nz + sinphi*tz;

	//the scattered vector cos(theta)*r+sin(theta)*s
	costheta = cos(theta);
	sintheta = sin(theta);
	sx = costheta*rx + sintheta*sx;
	sy = costheta*ry + sintheta*sy;
	sz = costheta*rz + sintheta*sz;

	//return scattered ray cosines to ZEMAX
	data[35] = sx;
	data[36] = sy;
	data[37] = sz;

	data[12] = 1.0;
	return 0;
}

int __declspec(dllexport) APIENTRY UserParamNames(char *data)
{
	/* this function returns the name of the parameter requested */
	int i;
	i = (int)data[0];
	strcpy(data, "");
	if (i == 1) strcpy(data, "Abs. Max");
	if (i == 2) strcpy(data, "Abs. 1/e Width");
	if (i == 3) strcpy(data, "PLQY");
	if (i == 4) strcpy(data, "Virtual Abs. Max");
	return 0;
}

int __declspec(dllexport) APIENTRY UserModelInformation(int *settings)
{
	/* this function returns whether standard model parameters are required
	settings[0] = size of settings array (should be 3)
	settings[1] = mean path required? 1 = yes, 0 = no
	settings[2] = angle required? 1 = yes, 0 = no
	*/

	if (settings[0] == 3)
	{
		settings[1] = 1;
		settings[2] = 0;
	}

	return 0;
}

void CrossProduct(double x1, double y1, double z1, double x2, double y2, double z2, double *x3, double *y3, double *z3)
{
	*x3 = y1*z2 - z1*y2;
	*y3 = z1*x2 - x1*z2;
	*z3 = x1*y2 - y1*x2;
}

void Normalize(double *x, double *y, double *z)
{
	double temp;
	temp = (*x)*(*x) + (*y)*(*y) + (*z)*(*z);
	temp = sqrt(temp);
	if (temp == 0) return;
	temp = 1.0 / temp;
	*x *= temp;
	*y *= temp;
	*z *= temp;
}



