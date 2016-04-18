#include "process.h"

#define SMALL_DATA 0

Process::Process(int pic) {
    if(pic == 1)
        path = "../data/lena.jpg";
    else if(pic == 2)
        path = "../data/baboo.bmp";
    else if(pic == 3)
        path = "../data/peppers512.bmp";
    else if(pic == 4)
        path = "../data/bridge512.jpg";
    else if(pic == 5)
        path = "../data/boats512.bmp";
    else if(pic == 6)
        path = "../data/cell.bmp";
    else if(pic == 7)
        path = "../data/couple.bmp";
    else if(pic == 8)
        path = "../data/elaine.bmp";
    else if(pic == 9)
        path = "../data/jokul.bmp";
    else if(pic == 10)
        path = "../data/man.bmp";
    else if(pic == 11)
        path = "../data/photography.bmp";
    else if(pic == 12)
        path = "../data/satellite.bmp";
    else if(pic == 13)
        path = "../data/satellite2.bmp";
    else if(pic == 14)
        path = "../data/saturn.bmp";
    else if(pic == 15)
        path = "../data/terrace.bmp";
    else if(pic == 16)
        path = "../data/testpat.bmp";
    else if(pic == 17)
        path = "../data/barbara.png";
    else if(pic == 18)
        path = "../data/xb50.jpg";
    else if(pic == 19)
        path = "../data/IMG.jpg";
    else if(pic == 20)
        path = "../data/dwi.bmp";
}

Process::Process(string pic) {
    path = pic;
}

Process::~Process() {
    /*
    if(initImageData != NULL)
        initImageData->Delete();
    if(imageData != NULL)
        imageData->Delete();
    */
}

void Process::debug(vtkImageData *inputData) {/*{{{*/
    string writePath = "./";
    writePath += getTime();
    ofstream of(writePath.c_str(), ios::out);
    of << dim[0] << " " << dim[1] << " " << dim[2] << "!!!" << endl;
    int mx = 0;
    for(int i = 0; i < dim[0]; i ++) {
        for(int j = 0; j < dim[1]; j ++) {
            unsigned char *tmp = static_cast<unsigned char*>(inputData->GetScalarPointer(i, j, 0));
            //unsigned short *tmp = static_cast<unsigned short*>(inputData->GetScalarPointer(i, j, 0));
            of << (int)*tmp << " ";
            mx = max(mx, (int)*tmp);
        }
        of << endl;
    }
    of << "max: " << mx << endl;
    of.close();
    sleep(2);
}/*}}}*/

string Process::getTime() {/*{{{*/
    char tmp[100];
    time_t t = time(0);
    strftime( tmp, sizeof(tmp), "log%Y-%m-%d_%H:%M:%S", localtime(&t)); 
    return string(tmp);
}/*}}}*/

void Process::render() {/*{{{*/
    //undo
    vtkSmartPointer<vtkImageViewer2> imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
    //imageViewer->SetInputConnection( imageData->GetOutputPort() );
#if VTK_MAJOR_VERSION <= 5
    imageViewer->SetInput( imageData );
#else
    imageViewer->SetInputData( imageData );
#endif
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
    imageViewer->SetupInteractor(renderWindowInteractor);
    imageViewer->Render();
    imageViewer->GetRenderer()->ResetCamera();
    imageViewer->Render();
    
    renderWindowInteractor->Start();
}/*}}}*/

void Process::readImage() {/*{{{*/
    vtkImageReader2 *reader;
    if(path[path.size() - 2] == 'p')
        reader = vtkJPEGReader::New();
    else if(path[path.size() - 2] == 'n')
        reader = vtkPNGReader::New();
    else if(path[path.size() - 1] == 'p')
        reader = vtkBMPReader::New();
    else if(path[path.size() - 1] == 'm')
        reader = vtkDICOMImageReader::New();
    reader->SetFileName(path.c_str());
    reader->Update();

    initImageData = vtkImageData::New();
    initImageData = reader->GetOutput();
    initImageData->GetDimensions(dim);

    imageData = vtkImageData::New();
    if(SMALL_DATA) {
        dim[0] = dim[1] = 10, dim[2] = 1;
        imageData->SetDimensions(dim);
        imageData->SetExtent(0, 9, 0, 9, 0, 0);
    }
    else {
        imageData->SetDimensions(dim);
        imageData->SetExtent(initImageData->GetExtent());
    }
    imageData->SetSpacing(initImageData->GetSpacing());
#if VTK_MAJOR_VERSION <= 5
    imageData->SetScalarType(VTK_INT);
    imageData->SetNumberOfScalarComponents(1);
#else
    imageData->AllocateScalars(VTK_INT, 1);
#endif
    
    cout << dim[0] << " " << dim[1] << " " << dim[2] << endl;
    //debug(initImageData);
    if(SMALL_DATA) {
        for(int i = 0; i < dim[0]; i ++) {
            for(int j = 0; j < dim[1]; j ++) {
                int *tmp = static_cast<int*>(imageData->GetScalarPointer(i, j, 0));
                unsigned char *tmp2 = static_cast<unsigned char*>(initImageData->GetScalarPointer(i + 100, j + 100, 0));
                *tmp = (int)*tmp2;
            }
        }
        debug(imageData);
    }
    else {
        for(int i = 0; i < dim[0]; i ++) {
            for(int j = 0; j < dim[1]; j ++) {
                int *tmp = static_cast<int*>(imageData->GetScalarPointer(i, j, 0));
                //unsigned short *tmp2 = static_cast<unsigned short*>(initImageData->GetScalarPointer(i, j, 0));
                unsigned char *tmp2 = static_cast<unsigned char*>(initImageData->GetScalarPointer(i, j, 0));
                *tmp = (int)*tmp2;
            }
        }
    }
    //cout << writePic_ssim() << endl;

}/*}}}*/

void Process::addGaussianNoise(double segema) {/*{{{*/
    //srand(unsigned(time(0)));
    for(int i = 0; i < dim[0]; i ++) {
        for(int j = 0; j < dim[1]; j ++) {
            int *tmp =  static_cast<int*>(imageData->GetScalarPointer(i, j, 0));
            double tt = vtkMath::Gaussian(*tmp, segema);
            //tt = (tt >= upBoard ? upBoard : tt);
            //tt = (tt >= 0 ? tt :  0);
            *tmp = int(tt);
        }
    }
}/*}}}*/

void Process::addImpluseNoise(double p) {/*{{{*/
    //srand(unsigned(time(0)));
    for(int i = 0; i < dim[0]; i ++) {
        for(int j = 0; j < dim[1]; j ++) {
            int *tmp =  static_cast<int*>(imageData->GetScalarPointer(i, j, 0));
            int x = rand() % 100;
            if(x >= 0 && x < 100 * p)
                *tmp = rand() % upBoard;
        }
    }

    if(SMALL_DATA)
        debug(imageData);
}/*}}}*/

void Process::addRicianNoise(double level, double &level_p) {
    /*level = 1,2,..15*/
    int mx = 0;
    for(int i = 0; i < dim[0]; i ++)
        for(int j = 0; j < dim[1]; j ++) {
            int *tmp =  static_cast<int*>(imageData->GetScalarPointer(i, j, 0));
            mx = max(mx, *tmp);
        }
    level_p = level * mx / 100;
    
    for(int i = 0; i < dim[0]; i ++)
        for(int j = 0; j < dim[1]; j ++) {
            int *tmp =  static_cast<int*>(imageData->GetScalarPointer(i, j, 0));
            double tmp1 = (*tmp / sqrt(2.0) + level_p * vtkMath::Gaussian(0, 1));
            double tmp2 = (*tmp / sqrt(2.0) + level_p * vtkMath::Gaussian(0, 1));
            double tt = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
            *tmp = int(tt);
        }
}



void Process::caluSNR_PSNR(double m, string filename) {
    if(SMALL_DATA)
        debug(imageData);
    double ISum = 0, IMinusKSum = 0;
    for(int i = 0; i < dim[0]; i ++) {
        for(int j = 0; j < dim[1]; j ++) {
            unsigned char *node0 = NULL;
            if(SMALL_DATA)
                node0 =  static_cast<unsigned char*>(initImageData->GetScalarPointer(i + 100, j + 100, 0));
            else
                node0 =  static_cast<unsigned char*>(initImageData->GetScalarPointer(i, j, 0));
            int *node =  static_cast<int *>(imageData->GetScalarPointer(i, j, 0));
            
            int node0_int = *node0;
            ISum += (node0_int) * (node0_int);
            IMinusKSum += (node0_int - *node) * (node0_int - *node);
        }
    }
    double SNR = 10 * log(ISum / IMinusKSum) / log(10.0);
    double MSE = IMinusKSum / (dim[0] * dim[1]);
    double PSNR = 10 * log(upBoard * upBoard * 1.0 / MSE) / log(10.0);

    //string path2 = path;
    //path2.insert(path2.find_last_of("."), "_new");
    string path2 = "new.png";
    writePic(imageData, path2);
    double ssim = calcu_ssim(path2);

    ofstream of(filename.c_str(), ios::out | ios::app);
    of << "SNR:" << SNR << " PSNR:" << PSNR << " SSIM:" << ssim << endl;
    of.close();
}

void Process::setImageData(vtkImageData *inputData) {
    //pay attention, not delete old imageData
    if(imageData != NULL)
        imageData->Delete();
    //imageData = vtkImageData::New();
    imageData = inputData;
}

vtkImageData *Process::getImageData() {
    return imageData;
}

void Process::writePic(vtkImageData *imageData, string path2) {
    for(int i = 0; i < dim[0]; i ++) {
        for(int j = 0; j < dim[1]; j ++) {
            int *tmp =  static_cast<int*>(imageData->GetScalarPointer(i, j, 0));
            double tt = *tmp;
            tt = (tt >= upBoard ? upBoard - 1 : tt);
            tt = (tt >= 0 ? tt :  0);
            /*
            if(*tmp < 0)
                cout << *tmp << "->" << tt << endl;
            if(*tmp >= 256)
                cout << *tmp << "->" << tt << endl;
                */
            *tmp = int(tt);
        }
    }

    /*cast and write*/
    vtkSmartPointer<vtkImageCast> castFilter = vtkSmartPointer<vtkImageCast>::New();
    castFilter->SetOutputScalarTypeToUnsignedChar();
#if VTK_MAJOR_VERSION <= 5
    castFilter->SetInput(imageData);
#else
    castFilter->SetInputData(imageData);
#endif
    castFilter->Update();

    /*there is a bug in vtkJPEGWriter, so i use the vtkBMPWriter*/
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(path2.c_str());
    writer->SetInputConnection(castFilter->GetOutputPort());
    writer->Write();
}

double Process::calcu_ssim(string path2) {
    /*compute the ssim*/
    string cmd = "";
    cmd = cmd + "pyssim " + path + " " + path2;
    FILE * fp;
    char buffer[80];
    fp = popen(cmd.c_str(), "r");
    fgets(buffer, sizeof(buffer), fp);
    pclose(fp);

    return atof(buffer); 
}
