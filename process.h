#ifndef PROCESS_H
#define PROCESS_H

#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolume.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include <vtkDICOMImageReader.h>
#include "vtkSmartPointer.h" 
#include "vtkImageCast.h"
#include "vtkImageData.h"
#include "vtkMath.h"
#include "vtkImageViewer2.h"
#include "vtkImageReader.h"
#include "vtkJPEGReader.h"
#include "vtkJPEGWriter.h"
#include "vtkBMPWriter.h"
#include "vtkBMPReader.h"
#include "vtkPNGReader.h"
#include "vtkPNGWriter.h"

#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <unistd.h>

using namespace std;

class Process {
public:
    Process(int);
    Process(string);
    ~Process();
    void render();
    void readImage();
    void addGaussianNoise(double segma);
    void addImpluseNoise(double) ;
    void addRicianNoise(double level, double &level_p);
    void caluSNR_PSNR(double, string);
    void setImageData(vtkImageData *);
    vtkImageData *getImageData();
    void writePic(vtkImageData *, string);
    
private:
    //static const int upBoard = 1 << 11;
    static const int upBoard = 256; 
    string path;
    vtkImageData *initImageData;
    string getTime();
    int dim[3];
    vtkImageData *imageData;

    void debug(vtkImageData *inputData);
    double calcu_ssim(string);
};

#endif //PROCESS_H
