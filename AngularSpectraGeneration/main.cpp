// c++ program to compute the C_ell's for the ISW power spectrum
// Elizabeth Brown
// Compile with: g++ -std=c++11 -lgsl -o PowerSpec -pthread -O3 main.cpp
// Run with: ./PowerSpec 40 2000 4 0.15 4 /cosma6/data/dp004/dc-brow6/InputFiles/Forecasting/parameterFile.csv None None
// or: ./PowerSpec 40 2000 4 0.15 4 /Users/cb607/Documents/Power_Spectrum_Calc/InputFiles/Forecasting/parameterFile.csv None None
#include "stdafx.h"
#include "ISWIntegrator.hpp"


///////////////
// FUNCTIONS //
///////////////

// Get length of array in a file
int getDataLength(std::string fileName) {
    std::ifstream file(fileName);
    if(!file.good()) {
        std::cerr << "File is bad" << std::endl;
        exit(EXIT_FAILURE);
    }
    int N = 0;
    while(!file.eof()) {
        std::string line;
        file >> line;
        N++;
    }
    file.close();
    N--;
    return N;
}

// read out data for 1d and 2d arrays
double * read1dDatFile(std::string fileName) {
    std::ifstream file(fileName);
    if(!file.good()) {
        std::cerr << "File is bad" << std::endl;
        exit(EXIT_FAILURE);
    }
    int N = 0;

    while(!file.eof()) {
        std::string line;
        file >> line;
        N++;
    }
    file.close();
    double * Data = new double[N-1];
    std::ifstream File(fileName);
    int i = 0;
    while(!File.eof()) {
        std::string line;
        File >> line;
        if(!line.empty()) {
            double dat = std::stod(line);
            Data[i] = dat;
        }
        i++;
    }
    File.close();
    return Data;
}

double** read2dDatFile(std::string fileName, int a_length, int k_length) {
    std::ifstream file;
    file.open(fileName);
    if(!file.good()) {
        std::cerr << "File is bad." << std::endl;
        exit(EXIT_FAILURE);
    }
    // This part needs changing if file is not 354x217
    double ** Data = new double* [a_length];
    for (int i = 0; i<a_length; i++) {
        Data[i] = new double[k_length];
    }
    double * data = new double [a_length*k_length];
    std::cout<<a_length*k_length<<std::endl;
    int i = 0;
    while(!file.eof()) {
        std::string line;
        file >> line;
        if(!line.empty()) {
            double dat = std::stod(line);
            data[i] = dat;
            i++;
        }
    }
    file.close();
    for(int i = 0; i<a_length; i++) {
        for(int j = 0; j<k_length; j++) {
            Data[i][j] = data[j+(i*k_length)];
        }
    }
    delete[] data;
    return Data;
}

void adjustArray(double* x, double* y, int length, double limit, int index) {
    std::pair<gsl_spline *, gsl_interp_accel *> spline;
    spline = spline1dSteffenData(x, y, length);
    y[index] = gsl_spline_eval(spline.first, limit, spline.second);
    gsl_spline_free (spline.first);
    gsl_interp_accel_free (spline.second);
    return;
}


parameterFile readParamFile(std::string fileName) {
    std::ifstream file;
    file.open(fileName);
    if(!file.good()) {
        std::cerr << "File is bad." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string lineInput, lineParamNum, lineh, lineParams, lineGalNum, lineGalBins, lineLensNum, lineLensBins; 
    std::string lineCrossCov, lineKernInst, lineNumObs, lineSpeedUp, lineObs, lineOutput, nullLine;
    std::string delimiter = ",";
    // Retrieve input file directory
    file >> lineInput;
    size_t pos = 0;
    int i = 0;
    std::string inputFileDirectory;
    while((pos=lineInput.find(delimiter)) != std::string::npos) {
        if(i==1) {
            inputFileDirectory = lineInput.substr(0, pos);
        }
        lineInput.erase(0, pos+delimiter.length());
        i++;
    }
    // Retrieve number of parameters
    file >> lineParamNum;
    int numParams = 0;
    pos = 0;
    i = 0;
    while((pos=lineParamNum.find(delimiter)) != std::string::npos) {
        if(i==1) {
            numParams = std::stoi(lineParamNum.substr(0, pos));
        }
        lineParamNum.erase(0, pos+delimiter.length());
        i++;
    }
    // Clear header line
    file >> nullLine;
    // Add parameter file name endings to array
    std::string * paramNames = new std::string[numParams];
    // Add T/F whether parameter derivative is calculated to array
    bool * calcParam = new bool[numParams];
    for(int j = 0; j<numParams; j++) {
        pos = 0;
        i = 0;
        file >> lineParams;
        while((pos=lineParams.find(delimiter)) != std::string::npos) {
            if(i==0) {
                paramNames[j] = lineParams.substr(0, pos);
            }
            else if(i==1) {
                calcParam[j] = std::stoi(lineParams.substr(0, pos));
            }
            lineParams.erase(0, pos+delimiter.length());
            i++;
        }
    }
    // Retrieve hubble parameter info
    double h; double hDev;
    pos = 0;
    i = 0;
    file>>lineh;
    while((pos=lineh.find(delimiter)) != std::string::npos) {
        if(i==1) {
            h = std::stod(lineh.substr(0, pos));
            lineh.erase(0, pos+delimiter.length());
            pos=lineh.find(delimiter);
            hDev = std::stod(lineh.substr(0, pos));
        }
        lineh.erase(0, pos+delimiter.length());
        i++;
    }
    // Retrieve number of galaxy bins
    int numGal = 0;
    pos = 0;
    i = 0;
    file >> lineGalNum;
    while((pos=lineGalNum.find(delimiter)) != std::string::npos) {
        if(i==1) {
            numGal = std::stoi(lineGalNum.substr(0, pos));
        }
        lineGalNum.erase(0, pos+delimiter.length());
        i++;
    }
    double* galFSky = new double[numGal];
    double* galNumDen = new double[numGal];
    double* galNumDenDev = new double[numGal];
    double* galZLow = new double[numGal];
    double* galZHigh = new double[numGal];
    double* galZMed = new double[numGal];
    std::string* galSuffix = new std::string[numGal];
    // Clear header line
    file >> nullLine;
    for(int k = 0; k<numGal; k++) {
        pos = 0;
        i = 0;
        file >> lineGalBins;
        while((pos=lineGalBins.find(delimiter)) != std::string::npos) {
            if(i==1) {
                galFSky[k] = std::stod(lineGalBins.substr(0, pos));
            }
            else if(i==2) {
                galNumDen[k] = std::stod(lineGalBins.substr(0, pos));
            }
            else if(i==3) {
                galNumDenDev[k] = std::stod(lineGalBins.substr(0, pos));
            }
            else if(i==4) {
                galZLow[k] = std::stod(lineGalBins.substr(0, pos));
            }
            else if(i==5) {
                galZHigh[k] = std::stod(lineGalBins.substr(0, pos));
            }
            else if(i==6) {
                galZMed[k] = std::stod(lineGalBins.substr(0, pos));
                lineGalBins.erase(0, pos+delimiter.length());
                pos=lineGalBins.find(delimiter);
                galSuffix[k] = lineGalBins.substr(0, pos);
            }
            lineGalBins.erase(0, pos+delimiter.length());
            i++;
        }
    }
    // Retrieve number of lensing bins
    int numLens = 0;
    pos = 0;
    i = 0;
    file >> lineLensNum;
    while((pos=lineLensNum.find(delimiter)) != std::string::npos) {
        if(i==1) {
            numLens = std::stoi(lineLensNum.substr(0, pos));
        }
        lineLensNum.erase(0, pos+delimiter.length());
        i++;
    }
    double* lensFSky = new double[numLens];
    double* lensNumDen = new double[numLens];
    double* lensNumDenDev = new double[numLens];
    double* lensZLow = new double[numLens];
    double* lensZHigh = new double[numLens];
    double* lensZMed = new double[numLens];
    std::string* lensSuffix = new std::string[numLens];
    // Clear header line
    file >> nullLine;
    for(int k = 0; k<numLens; k++) {
        pos = 0;
        i = 0;
        file >> lineLensBins;
        while((pos=lineLensBins.find(delimiter)) != std::string::npos) {
            if(i==1) {
                lensFSky[k] = std::stod(lineLensBins.substr(0, pos));
            }
            else if(i==2) {
                lensNumDen[k] = std::stod(lineLensBins.substr(0, pos));
            }
            else if(i==3) {
                lensNumDenDev[k] = std::stod(lineLensBins.substr(0, pos));
            }
            else if(i==4) {
                lensZLow[k] = std::stod(lineLensBins.substr(0, pos));
            }
            else if(i==5) {
                lensZHigh[k] = std::stod(lineLensBins.substr(0, pos));
            }
            else if(i==6) {
                lensZMed[k] = std::stod(lineLensBins.substr(0, pos));
                lineLensBins.erase(0, pos+delimiter.length());
                pos=lineLensBins.find(delimiter);
                lensSuffix[k] = lineLensBins.substr(0, pos);
            }
            lineLensBins.erase(0, pos+delimiter.length());
            i++;
        }
    }
    // T/F to calculate galaxy bin cross-covariance
    bool calcCrossCov = 0;
    pos = 0;
    i = 0;
    file >> lineCrossCov;
    while((pos=lineCrossCov.find(delimiter)) != std::string::npos) {
        if(i==1) {
            calcCrossCov = std::stoi(lineCrossCov.substr(0, pos));
        }
        lineCrossCov.erase(0, pos+delimiter.length());
        i++;
    }
    // Clear header line
    file >> nullLine;
    // Assign kernels to instruments
    int* kernInst = new int[6];
    pos = 0;
    i = 0;
    file >> lineKernInst;
    while((pos=lineKernInst.find(delimiter)) != std::string::npos) {
        if(i<6) {
            kernInst[i] = std::stoi(lineKernInst.substr(0, pos));
        }
        lineKernInst.erase(0, pos+delimiter.length());
        i++;
    }
    // Retrieve number of observables
    int numObs = 0;
    pos = 0;
    i = 0;
    file >> lineNumObs;
    while((pos=lineNumObs.find(delimiter)) != std::string::npos) {
        if(i==1) {
            numObs = std::stoi(lineNumObs.substr(0, pos));
        }
        lineNumObs.erase(0, pos+delimiter.length());
        i++;
    }
    // T/F to speed up calculation
    bool calcSpeedUp = 0;
    pos = 0;
    i = 0;
    file >> lineSpeedUp;
    while((pos=lineSpeedUp.find(delimiter)) != std::string::npos) {
        if(i==1) {
            calcSpeedUp = std::stoi(lineSpeedUp.substr(0, pos));
        }
        lineSpeedUp.erase(0, pos+delimiter.length());
        i++;
    }
    // Clear header line
    file >> nullLine;
    // Add observable file suffix to array
    std::string * observableNames = new std::string[numObs];
    // Add T/F to calculate observable to array
    bool * calcObservable = new bool[numObs];
    // Array for number of kernels to be calculated
    int ** calcKernel = new int*[numObs];
    for(int j = 0; j<numObs; j++) {
        calcKernel[j] = new int[6];
        pos = 0;
        i = 0;
        file >> lineObs;
        while((pos=lineObs.find(delimiter)) != std::string::npos) {
            if(i==0) {
                observableNames[j] = lineObs.substr(0, pos);
            }
            else if(i==1) {
                calcObservable[j] = std::stoi(lineObs.substr(0, pos));
            }
            else if(i==2) {
                calcKernel[j][0] = std::stoi(lineObs.substr(0, pos));
            }
            else if(i==3) {
                calcKernel[j][1] = std::stoi(lineObs.substr(0, pos));
            }
            else if(i==4) {
                calcKernel[j][2] = std::stoi(lineObs.substr(0, pos));
            }
            else if(i==5) {
                calcKernel[j][3] = std::stoi(lineObs.substr(0, pos));
            }
            else if(i==6) {
                calcKernel[j][4] = std::stoi(lineObs.substr(0, pos));
                lineObs.erase(0, pos+delimiter.length());
                pos=lineObs.find(delimiter);
                calcKernel[j][5] = std::stoi(lineObs.substr(0, pos));
            }
            lineObs.erase(0, pos+delimiter.length());
            i++;
        }
    }
    // Retrieve output file directory
    std::string outputFileDirectory;
    pos = 0;
    i = 0;
    file >> lineOutput;
    while((pos=lineOutput.find(delimiter)) != std::string::npos) {
        if(i==1) {
            outputFileDirectory = lineOutput.substr(0, pos);
        }
        lineOutput.erase(0, pos+delimiter.length());
        i++;
    }
    file.close();
    // Add above to parameter file structure
    struct parameterFile paramFile = {inputFileDirectory, numParams, paramNames, calcParam, h, hDev, numGal,
                                        galFSky, galNumDen, galNumDenDev, galZLow, galZHigh, galZMed, galSuffix, 
                                        numLens, lensFSky, lensNumDen, lensNumDenDev, lensZLow,
                                        lensZHigh, lensZMed, lensSuffix, calcCrossCov, kernInst, numObs, calcSpeedUp,
                                        observableNames, calcObservable, calcKernel, outputFileDirectory};
    return paramFile;
}

// Function outputting data
void outputIntData(std::string fileName, int * x, double * y, int length) {
    std::ofstream file(fileName);
    for (int i = 0; i < length-1; i++) {
        file << x[i] << "," << y[i] << std::endl;
    }
    file << x[length-1] << "," << y[length-1];
    file.close();
}

void outputDoubleData(std::string fileName, double * x, double * y, int length) {
    std::ofstream file(fileName);
    for (int i = 0; i < length-1; i++) {
        file << std::setprecision(6) << x[i] << "," << y[i] << std::endl;
    }
    file << std::setprecision(6) << x[length-1] << "," << y[length-1];
    file.close();
}

void splineDataOutput(std::pair<gsl_spline *, gsl_interp_accel *> spline, int size, std::string fileName, 
                        double xLow, double xHigh) {
    double* splineData = new double[size+1];
    double* xArray = new double[size+1];
    double stepSize = (xHigh - xLow) / size;
    for(int i = 0; i < size+1; i++) {
        splineData[i] = gsl_spline_eval(spline.first, (i*stepSize) + xLow, spline.second);
        xArray[i] = (i*stepSize) + xLow;
    }
    outputDoubleData(fileName, xArray, splineData, size+1);
    delete[] splineData;
    delete[] xArray;
}


int* kernelsToCalc(parameterFile paramFile) {
    int* kernels = new int[6]();
    kernels[0]=0; kernels[1]=0; kernels[2]=0; kernels[3]=0; kernels[4] = 0; kernels[5] = 0;
    for (int g = 0; g<paramFile.numObs; g++) {
        if(paramFile.calcObs[g]) {
            for(int i = 0; i<6; i++) {
                kernels[i] = kernels[i] + paramFile.calcKernel[g][i];
            }
        }
    }
    return kernels;
}

// Makes an array for each observable, showing which instruments must be calculated
int* kernInstToCalc(parameterFile paramFile) {
    int* kernInstArray = new int[paramFile.numObs];
    for(int o = 0; o<paramFile.numObs; o++) {
        kernInstArray[o] = 0;
        // Iterate over observables
        for(int i = 0; i<6; i++) {
            // if kernel is being calculated for this observable
            if(paramFile.calcKernel[o][i]>0) {
                // Unique combinations of numbers indicate instrument calculation type
                if(paramFile.calcKernel[o][i]==1) {
                    kernInstArray[o] = kernInstArray[o] + paramFile.kernInst[i];
                }
                else {
                    kernInstArray[o] = 2*paramFile.kernInst[i];
                }
            }
        }
    }
    return kernInstArray;
}

int** observableIndices(parameterFile paramFile) {
    int** indices = new int*[paramFile.numObs];
    for(int o = 0; o<paramFile.numObs; o++) {
        indices[o] = new int[2];
        int index = 0;
        int j = 0;
        for(int i = 0; i<6; i++) {
            // if kernel is calculated for this observable
            if(paramFile.calcKernel[o][i]>0) {
                // Cross-correlation
                if(paramFile.calcKernel[o][i]==1) {
                    indices[o][j]=index;
                    j++;
                }
                // Auto-correlation
                else {
                    indices[o][0]=index;
                    indices[o][1]=index;
                }
            }
            // increase index depending on which kernel we are currently on
            if(paramFile.kernInst[i]==0) {
                index++;
            }
            else if(paramFile.kernInst[i]==1) {
                index = index+paramFile.numGal;
            }
            else {
                index = index+paramFile.numLens;
            }
        }
    }
    return indices;
}

int observablesToCalc(parameterFile paramFile) {
    int galgalObs = 0;
    for(int g = 0; g<paramFile.numGal; g++) {
        for(int h = g; h<paramFile.numGal; h++) {
            galgalObs++;
        }
    }
    int numObs = 0;
    for(int i = 0; i<paramFile.numObs; i++) {
        if(paramFile.calcObs[i]) {
            if(paramFile.calcKernel[i][1]==2) {
                numObs = numObs + galgalObs;
            }
            else if(paramFile.calcKernel[i][1]==1) {
                numObs = numObs + paramFile.numGal;
            }
            else {
                numObs++;
            }
        }
    }
    return numObs;
}


double** ISWKernel(double** densityISW, double* fOfK, double* HubbleRate, int k_length, int a_length) {
    double** kernel = new double*[k_length];
    for(int j = 0; j<k_length; j++) {
        kernel[j] = new double[a_length];
        double* densityISWSlice = invertArray(densityISW[j], a_length);
        double* invertedHubble = invertArray(HubbleRate, a_length);
        for(int i = 0; i<a_length; i++) {
            kernel[j][i]=8.*PI*fOfK[j]*invertedHubble[i]*densityISWSlice[i]/pow(SPEED_OF_LIGHT, 3.);
        }
        delete[] densityISWSlice;
        delete[] invertedHubble;
    }
    return kernel;
}

double** GalKernel(double** densityGal, double* bias, double* selectionFunc, int k_length, int a_length) {
    double** kernel = new double*[k_length];
    for(int j = 0; j<k_length; j++) {
        kernel[j] = new double[a_length];
        double* densityGalSlice = invertArray(densityGal[j], a_length);
        double* invertedBias = invertArray(bias, a_length);
        double* invertedSelectionFunc = invertArray(selectionFunc, a_length);
        for(int i = 0; i<a_length; i++) {
            kernel[j][i] = 4.*PI*invertedBias[i]*invertedSelectionFunc[i]*densityGalSlice[i];
        }
        delete[] densityGalSlice;
        delete[] invertedBias;
        delete[] invertedSelectionFunc;
    }
    return kernel;
}

double** CMBLenseKernel(double** densityLense, double* CMBLenseKernel, int k_length, int a_length) {
    double** kernel = new double*[k_length];
    for(int j = 0; j<k_length; j++) {
        kernel[j] = new double[a_length];
        double* densityLenseSlice = invertArray(densityLense[j], a_length);
        for(int i = 0; i<a_length; i++) {
            kernel[j][i] = 8.*PI*CMBLenseKernel[i]*densityLenseSlice[i];
        }
        delete[] densityLenseSlice;
    }
    return kernel;
}

double** LenseKernel(double** density, double* fOfK, double* k, double* lenseKernel, int k_length, int a_length) {
    double** kernel = new double*[k_length];
    for(int j = 0; j<k_length; j++) {
        kernel[j] = new double[a_length];
        double* densitySlice = invertArray(density[j], a_length);
        for(int i = 0; i<a_length; i++) {
            kernel[j][i] = 4.*PI*fOfK[j]*k[j]*k[j]*lenseKernel[i]*densitySlice[i];
        }
        delete[] densitySlice;
    }
    return kernel;
}

double** GalLenseKernel(double* galNumDen, double** density, double* bias, int k_length, int a_length) {
    double** kernel = new double*[k_length];
    for(int j = 0; j<k_length; j++) {
        kernel[j] = new double[a_length];
        double* densitySlice = invertArray(density[j], a_length);
        double* invertedBias = invertArray(bias, a_length);
        for(int i = 0; i<a_length; i++) {
            kernel[j][i] = 4.*PI*galNumDen[i]*invertedBias[i]*densitySlice[i];
        }
        delete[] densitySlice;
        delete[] invertedBias;
    }
    return kernel;
}

double** IntAlignKernel(double** density, double* IAKernel, double* numDen, int k_length, int a_length) {
    double** kernel = new double*[k_length];
    for(int j = 0; j<k_length; j++) {
        kernel[j] = new double[a_length];
        double* densitySlice = invertArray(density[j], a_length);
        for(int i = 0; i<a_length; i++) {
            kernel[j][i] = 4.*PI*IAKernel[i]*numDen[i]*densitySlice[i];
        }
        delete[] densitySlice;
    }
    return kernel;
}



void getSpectraODE(parameterFile paramFile, int* &ellLims, int l_max1, int l_max2, int numThreads, double aLow, double aHigh, 
                    int BesselArray_size, double* kChi, std::string inputSuffix="", std::string suffix="", std::string outFileSuffix="") {
    // Output .dat files to 1D (2D) arrays
    std::string inputDirectory = paramFile.inputDirectory;
    std::string outputDirectory = paramFile.outputDirectory;
    const int k_size = getDataLength(inputDirectory + "wavenum"+inputSuffix+suffix+".dat");
    const int a_size = getDataLength(inputDirectory + "scale"+inputSuffix+suffix+".dat");
    double * scale;
    double ** galSelec;
    double ** galSelectionFunc = new double*[paramFile.numGal];
    double ** galNumDen = new double*[paramFile.numGal];
    for(int g = 0; g<paramFile.numGal; g++) {
        galNumDen[g] = new double[a_size];
    }
    double ** lensNumDen = new double*[paramFile.numLens];
    for(int g = 0; g<paramFile.numLens; g++) {
        lensNumDen[g] = new double[a_size];
    }
    double * loga = new double [a_size];
    // a2H is actually a2H/(H_0).  a2Hinv inverts a2H and multiplies by H_0.
    double * HubbleRate = new double [a_size];
    double * a2Hinv = new double [a_size];
    double * a2H;
    // k /Mpc.
    double * k;
    // A_s(k/k*)^(n-1)
    double * P_primordial;
    // density = SUM(fT), logGrowth = 1-SUM(dlogT/dloga)
    double ** densityGal;
    double ** densityISW;
    double ** densityCMBLense;
    double ** IAKernel;
    double ** nonLinPkCorrection;
    double ** nonLinGal;
    double ** nonLinISW;
    double ** nonLinCMBLense;
    // density, logGrowth transpose
    double ** densityGalT = new double* [k_size];
    double ** densityISWT = new double* [k_size];
    double ** densityCMBLenseT = new double* [k_size];
    double ** NLGalT = new double* [k_size];
    double ** NLISWT = new double* [k_size];
    double ** NLCMBLenseT = new double* [k_size];
    double ** NLCorrectionT = new double* [k_size];
    for (int i = 0; i<k_size; i++) {
        densityGalT[i] = new double[a_size];
        densityISWT[i] = new double[a_size];
        densityCMBLenseT[i] = new double[a_size];
        NLGalT[i] = new double[a_size];
        NLISWT[i] = new double[a_size];
        NLCMBLenseT[i] = new double[a_size];
        NLCorrectionT[i] = new double[a_size];
    }
    double ** bias;
    double * F_of_k;
    // Kernel calculation T/F
    int* kernels = kernelsToCalc(paramFile);
    int* kernelInstruments = kernInstToCalc(paramFile);
    int** obsIndices = observableIndices(paramFile);
    int numObs = observablesToCalc(paramFile);
    // Read out data from files.  Change path if needed.
    scale = read1dDatFile(inputDirectory + "scale"+inputSuffix+suffix+".dat");
    k = read1dDatFile(inputDirectory + "wavenum"+inputSuffix+suffix+".dat");
    a2H = read1dDatFile(inputDirectory + "a2Hubble"+inputSuffix+suffix+".dat");
    P_primordial = read1dDatFile(inputDirectory + "P_primordial"+inputSuffix+suffix+".dat");
    bias = read2dDatFile(inputDirectory + "bias"+inputSuffix+suffix+".dat", paramFile.numGal, a_size);
    galSelec = read2dDatFile(inputDirectory + "galSelecFunc"+inputSuffix+suffix+".dat", paramFile.numGal, a_size);
    F_of_k = read1dDatFile(inputDirectory + "FofK"+inputSuffix+suffix+".dat");
    IAKernel = read2dDatFile(inputDirectory + "IAKernel"+inputSuffix+suffix+".dat", paramFile.numLens, a_size);
    densityGal = read2dDatFile(inputDirectory + "densityGal"+inputSuffix+suffix+".dat", a_size, k_size);
    densityISW = read2dDatFile(inputDirectory + "densityISW"+inputSuffix+suffix+".dat", a_size, k_size);
    densityCMBLense = read2dDatFile(inputDirectory + "densityLense"+inputSuffix+suffix+".dat", a_size, k_size);
    nonLinGal = read2dDatFile(inputDirectory + "nonLinDensityGal"+inputSuffix+suffix+".dat", a_size, k_size);
    nonLinISW = read2dDatFile(inputDirectory + "nonLinDensityISW"+inputSuffix+suffix+".dat", a_size, k_size);
    nonLinCMBLense = read2dDatFile(inputDirectory + "nonLinDensityLense"+inputSuffix+suffix+".dat", a_size, k_size);
    nonLinPkCorrection = read2dDatFile(inputDirectory + "nonLinPkCorrection"+inputSuffix+suffix+".dat", a_size, k_size);


    // Get transpose of density transfers
    for (int i = 0; i<a_size; i++) {
        for (int j = 0; j<k_size; j++) {
            densityGalT[j][i] = densityGal[i][j];
            densityISWT[j][i] = densityISW[i][j];
            densityCMBLenseT[j][i] = densityCMBLense[i][j]/(k[j]*k[j]);
            NLGalT[j][i] = nonLinGal[i][j];
            NLISWT[j][i] = nonLinISW[i][j];
            NLCMBLenseT[j][i] = nonLinCMBLense[i][j]/(k[j]*k[j]);
            NLCorrectionT[j][i] = nonLinPkCorrection[i][j];
        }
    }
    for(int i = 0; i<a_size; i++) {
        delete[] densityISW[i];
        delete[] densityGal[i];
        delete[] densityCMBLense[i];
        delete[] nonLinISW[i];
        delete[] nonLinGal[i];
        delete[] nonLinCMBLense[i];
        delete[] nonLinPkCorrection[i];
    }
    delete[] densityISW;
    delete[] densityGal;
    delete[] densityCMBLense;
    delete[] nonLinISW;
    delete[] nonLinGal;
    delete[] nonLinCMBLense;
    delete[] nonLinPkCorrection;

    // Sort suffix dependent parameters
    double hubbleParam = paramFile.h;
    if(suffix=="h+") {
        hubbleParam = hubbleParam+paramFile.hDev;
    }
    else if(suffix=="h-") {
        hubbleParam = hubbleParam-paramFile.hDev;
    }

    // calculate a spline of chi(a) using a2Hinv, and adjust scale datasets
    for (int i = 0; i < a_size; i++) {
        a2H[i] = (H_0*a2H[i]);
        loga[i] = log(scale[i]);
    }
    for (int i = 0; i < a_size; i++) {
        a2Hinv[i] = SPEED_OF_LIGHT / (a2H[i]);
        HubbleRate[i] = a2H[i]/(scale[i]*scale[i]);
    }
    for(int i = 0; i < k_size; i++) {
        F_of_k[i]=F_of_k[i]*H_0*H_0;
    }
    std::pair<gsl_spline *, gsl_interp_accel *> chi;
    chi = chiSpline(scale, loga, a2Hinv, a_size);
    double* invertedChiData = retrieveSplineArray(chi, loga, a_size);

    // Adjust arrays to avoid singular z=0
    double zLow = (1/exp(aHigh))-1;
    double zHigh = (1/exp(aLow))-1;
    double zMin = (1/scale[a_size-1])-1;
    if(zMin<Z_LOW) {
        std::cout<<"Adjusting scale arrays. z_min="<<Z_LOW<<std::endl;
        if(zLow<Z_LOW) {
            std::cout<<"Adjusting lower limit. z_low="<<10*Z_LOW<<std::endl;
            aHigh = log(1/((10*Z_LOW)+1));
        }
        adjustArray(scale, a2H, a_size, 1/(Z_LOW+1), a_size-1);
        adjustArray(scale, a2Hinv, a_size, 1/(Z_LOW+1), a_size-1);
        adjustArray(scale, HubbleRate, a_size, 1/(Z_LOW+1), a_size-1);
        for(int g = 0; g<paramFile.numGal; g++) {
            adjustArray(scale, bias[g], a_size, 1/(Z_LOW+1), a_size-1);
        }
        adjustArray(loga, invertedChiData, a_size, log(1/(Z_LOW+1)), a_size-1);
        for(int i = 0; i<k_size; i++) {
            adjustArray(scale, densityGalT[i], a_size, 1/(Z_LOW+1), a_size-1);
            adjustArray(scale, densityISWT[i], a_size, 1/(Z_LOW+1), a_size-1);
            adjustArray(scale, densityCMBLenseT[i], a_size, 1/(Z_LOW+1), a_size-1);
        }
        scale[a_size-1]=1/(Z_LOW+1);
        loga[a_size-1]=log(1/(Z_LOW+1));
    }
    gsl_spline_free (chi.first);
    gsl_interp_accel_free (chi.second);
    chi = spline1dSteffenData(loga, invertedChiData, a_size);

    double* chiData = invertArray(invertedChiData, a_size);
    std::cout<<chiData[0]<<", "<<chiData[a_size-1]<<std::endl;

    // Sort galaxy number density selection functions and number density arrays
    // Make selection function arrays

    double* normalisation = new double[paramFile.numGal];
    for(int i = 0; i<paramFile.numGal; i++) {
        normalisation[i] = renormalise(paramFile.galZMed[i]/1.2, 1./(paramFile.galZHigh[i]+1), 1./(paramFile.galZLow[i]+1)).first;
    }

    for(int g = 0; g < paramFile.numGal; g++) {
        double z0=paramFile.galZMed[g]/1.2;
        double* redshift = new double[a_size];
        for(int i = 0; i < a_size; i++) {  // Produce Z array for splining selection function
            redshift[i] = (1./scale[a_size-i-1])-1.;
        }
        std::pair<gsl_spline *, gsl_interp_accel *> SelectionFunction;
        SelectionFunction = spline1dSteffenData(redshift, galSelec[g], a_size);  // spline Selection Function
        double norm = gsl_spline_eval_integ(SelectionFunction.first, (1./(exp(aHigh)))-1, (1./(exp(aLow)))-1, SelectionFunction.second);  // Integrate over z-range
        galSelectionFunc[g] = invertArray(galSelec[g], a_size);  // invert for function to follow scale array
        for(int i = 0; i < a_size; i++) {  // normalise selection function and set as 1/N dN/dchi following scale array
            galSelectionFunc[g][i] = HubbleRate[i]*galSelectionFunc[g][i]/(SPEED_OF_LIGHT*norm);
        }
        delete[] redshift;  // Clean-up
        delete[] galSelec[g];
        gsl_spline_free (SelectionFunction.first);
        gsl_interp_accel_free (SelectionFunction.second);
        galNumDen[g] = invertArray(galSelectionFunc[g], a_size);
    }
    delete[] galSelec;

    // Sort lensing number density arrays
    // Make number density arrays
    double* redshift = new double[a_size];
    for(int i = 0; i < a_size; i++) {
        redshift[i] = (1/scale[a_size-i-1])-1;
    }
    double** numDenIntegral1 = new double*[paramFile.numLens];
    for(int g = 0; g < paramFile.numLens; g++) {
        numDenIntegral1[g] = new double[a_size];
        double z0=paramFile.lensZMed[g]/1.2;
        for(int i = 0; i < a_size; i++) {
            // Uses Euclid's parameters
            numDenIntegral1[g][i] = numDenIntegrationOriginal(redshift[i], z0, 1.0, 0.0, 0.05, 1.0, 0.1, 0.05, 0.1, 
                                                    paramFile.lensZLow[g], paramFile.lensZHigh[g], redshift, a_size);
        }
    }

    double ** lensNumDen1 = new double*[paramFile.numLens];
    for(int g = 0; g<paramFile.numLens; g++) {
        lensNumDen1[g] = new double[a_size];
    }

    double** numDenIntegral;
    numDenIntegral = read2dDatFile(inputDirectory + "sourceSelecFunc"+inputSuffix+suffix+".dat", paramFile.numLens, a_size);
    for(int g = 0; g < paramFile.numLens; g++) {
        std::pair<gsl_spline *, gsl_interp_accel *> numDenNormInt = spline1dSteffenData(redshift, numDenIntegral[g], a_size);
        double numDenNorm = gsl_spline_eval_integ(numDenNormInt.first, (1./(exp(aHigh)))-1, (1./(exp(aLow)))-1, numDenNormInt.second);
        gsl_spline_free (numDenNormInt.first);
        gsl_interp_accel_free (numDenNormInt.second);
        std::pair<gsl_spline *, gsl_interp_accel *> numDenNormInt1 = spline1dSteffenData(redshift, numDenIntegral1[g], a_size);
        double numDenNorm1 = gsl_spline_eval_integ(numDenNormInt1.first, (1./(exp(aHigh)))-1, (1./(exp(aLow)))-1, numDenNormInt1.second);
        gsl_spline_free (numDenNormInt1.first);
        gsl_interp_accel_free (numDenNormInt1.second);
        for(int i = 0; i < a_size; i++) {
            lensNumDen[g][i] = HubbleRate[a_size-i-1]*(numDenIntegral[g][i])/(numDenNorm*SPEED_OF_LIGHT);
            lensNumDen1[g][i] = HubbleRate[a_size-i-1]*(numDenIntegral1[g][i])/(numDenNorm1*SPEED_OF_LIGHT);
        }
        delete[] numDenIntegral[g];
        delete[] numDenIntegral1[g];
    }
    delete[] numDenIntegral;
    delete[] numDenIntegral1;
    delete[] redshift;
    std::cout<<"lensNumDen: "<<lensNumDen[0][10]<<"; lensNumDen1: "<<lensNumDen1[0][10]<<std::endl;
    for(int g = 0; g<paramFile.numLens; g++) {
        delete[] lensNumDen1[g];
    }
    delete[] lensNumDen1;
    
    double chiLow = gsl_spline_eval(chi.first, aHigh, chi.second);
    std::cout << "Maximum linear l = " << chiLow*K_NON_LIN - 0.5 << std::endl;
    double chiHigh = gsl_spline_eval(chi.first, aLow, chi.second);
    double chiCMB = gsl_spline_eval(chi.first, log(1./(1101.)), chi.second);
    double* CMBLensingKernel = new double[a_size];
    for(int i = 0; i < a_size; i++) {
        CMBLensingKernel[i] = (chiCMB-chiData[i])/(chiCMB*chiData[i]);
    }
    // Produce gravitational lensing source upper limits
    double* lensingSourceChi = new double[paramFile.numLens];
    double** galaxyBinChi = new double*[paramFile.numGal];
    double* invertedScale = invertArray(scale, a_size);
    for(int g = 0; g<paramFile.numLens; g++) {
        double chiSource = gsl_spline_eval(chi.first, log(1./(paramFile.lensZHigh[g]+1.)), chi.second);
        lensingSourceChi[g] = std::min(chiSource, chiHigh);
    }
    for(int g = 0; g<paramFile.numGal; g++) {
        galaxyBinChi[g] = new double[2];
        double chiSourceLow = gsl_spline_eval(chi.first, log(1./(paramFile.galZLow[g]+1.)), chi.second);
        double chiSourceHigh = gsl_spline_eval(chi.first, log(1./(paramFile.galZHigh[g]+1.)), chi.second);
        galaxyBinChi[g][0] = std::max(chiSourceLow, chiLow);
        std::cout<<"Galaxy "<<g<<", l_max = "<<galaxyBinChi[g][0]*K_LIMIT-0.5<<std::endl;
        galaxyBinChi[g][1] = std::min(chiSourceHigh, chiHigh);
    }
    

    // Produce gravitational lensing kernels
    double** lensingKernel = new double*[paramFile.numLens];
    for(int g = 0; g<paramFile.numLens; g++) {
        lensingKernel[g] = new double[a_size];
        for(int i = 0; i<a_size; i++) {
            double* lenseIntegrand = new double[a_size];
            bool firstTime = true;  // D. Bacon function
            double chiWeight = 0;  // D. Bacon function
            for(int j = 0; j<a_size; j++) {
                lenseIntegrand[j] = lensNumDen[g][j]*(chiData[j]-chiData[i])/(chiData[j]);
                if(firstTime && 1/(scale[a_size-j-1])-1>2.0) {  // D. Bacon function. z=2
                    chiWeight = (chiData[j]-chiData[i])/(chiData[j]);
                    firstTime = false;
                }
            }
            std::pair<gsl_spline *, gsl_interp_accel *> lenseInt = spline1dSteffenData(chiData, lenseIntegrand, a_size);
            delete[] lenseIntegrand;
            lensingKernel[g][i] = chiData[i]*gsl_spline_eval_integ(lenseInt.first, chiData[i], chiData[a_size-1], lenseInt.second)/invertedScale[i];
            //lensingKernel[g][i] = chiData[i]*chiWeight/invertedScale[i];  // D. Bacon function
            gsl_spline_free (lenseInt.first);
            gsl_interp_accel_free (lenseInt.second);
        }
    }

    delete[] invertedScale;
    

    // Initialise Limber power spectra arrays
    double** matterSpectrum = new double*[k_size];
    for(int i = 0; i<k_size; i++) {
        matterSpectrum[i] = new double[a_size];
    }
    double*** corrFuncs = new double**[2+(2*paramFile.numGal)+(2*paramFile.numLens)];
    for (int i = 0; i<2+(2*paramFile.numGal)+(2*paramFile.numLens); i++) {
        corrFuncs[i] = new double*[k_size];
        for(int j = 0; j<k_size; j++) {
            corrFuncs[i][j] = new double[a_size];
        }
    }

    // Calculate total number of spectra to be calculated
    int totalObs = 0;
    for(int obsType=0; obsType<paramFile.numObs; obsType++) {
        if(paramFile.calcObs[obsType]) {
            // Galaxy instrument kernel auto-spectra
            if(kernelInstruments[obsType]==2) {
                for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numGal; g2++) {  // column galaxies
                        if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                            if(g1==g2) {
                                totalObs++;
                            }
                        }
                        else {  // if calculate full covar
                            totalObs++;
                        }
                    }
                }
            }
            // Lensing instrument kernel auto-spectra
            else if(kernelInstruments[obsType]==6) {
                for(int g1=0; g1<paramFile.numLens; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numLens; g2++) {  // column galaxies
                        if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                            if(g1==g2) {
                                totalObs++;
                            }
                        }
                        else {  // if calculate full covar
                            totalObs++;
                        }
                    }
                }
            }
            // Lensing-Galaxy instrument kernels cross-spectra
            else if(kernelInstruments[obsType]==4) {
                for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numLens; g2++) {  // column lens
                        totalObs++;
                    }
                }
            }
            // Galaxy instrument cross-spectra (non-Lensing)
            else if(kernelInstruments[obsType]==1) {
                for(int g=0; g<paramFile.numGal; g++) {
                    totalObs++;
                }
            }
            // Lensing instrument cross-spectra (non-Galaxy)
            else if(kernelInstruments[obsType]==3) {
                for(int g=0; g<paramFile.numLens; g++) {
                    totalObs++;
                }
            }
            // Non-galaxy, non-lensing spectra
            else {
                totalObs++;
            }
        }
    }

    // Limber correlation functions
    matterPowerSpectrum(matterSpectrum, P_primordial, k, a_size, k_size);
    ISWPowerSpectrumCorr(corrFuncs[0], F_of_k, HubbleRate, densityISWT, a_size, k_size);
    CMBLensePowerSpectrumCorr(corrFuncs[1], CMBLensingKernel, densityCMBLenseT, a_size, k_size);
    for(int g = 0; g<paramFile.numGal; g++) {
        galPowerSpectrumCorr(corrFuncs[2+g], bias[g], galSelectionFunc[g], densityGalT, a_size, k_size);
        GalLensePowerSpectrumCorr(corrFuncs[2+paramFile.numGal+g], galNumDen[g], bias[g], densityGalT, a_size, k_size);
    }
    for(int g = 0; g<paramFile.numLens; g++) {
        LensePowerSpectrumCorr(corrFuncs[2+(2*paramFile.numGal)+g], lensingKernel[g], densityGalT, F_of_k, k, a_size, k_size);
        IntAlignPowerSpectrumCorr(corrFuncs[2+(2*paramFile.numGal)+paramFile.numLens+g], IAKernel[g], lensNumDen[g], 
                                    densityGalT, a_size, k_size);
    }

    // Limber splines, linear spectra
    std::pair<gsl_spline *, gsl_interp_accel *>** LimberSpectra = new std::pair<gsl_spline *, gsl_interp_accel *>*[totalObs];
    double** integrationLimitsLimber = new double*[totalObs];
    int* l_lims = new int[totalObs];
    for(int i = 0; i < totalObs; i++) {
        integrationLimitsLimber[i] = new double[2];
    }
    totalObs = 0;
    for(int obsType=0; obsType<paramFile.numObs; obsType++) {  // Generate Limber power spectra
        if(paramFile.calcObs[obsType]) {
            // Galaxy auto-spectra
            if(kernelInstruments[obsType]==2) {
                for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numGal; g2++) {  // column galaxies
                        if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                            if(g1==g2) {
                                LimberSpectra[totalObs] = powerSpectrum(matterSpectrum, corrFuncs[obsIndices[obsType][0]+g1], corrFuncs[obsIndices[obsType][1]+g1], k, a_size, k_size);
                                integrationLimitsLimber[totalObs][0] = std::max(log(1./(paramFile.galZHigh[g1]+1)), aLow);
                                integrationLimitsLimber[totalObs][1] = std::min(log(1./(paramFile.galZLow[g1]+1)), aHigh);
                                l_lims[totalObs] = int(round((K_LIMIT*gsl_spline_eval(chi.first, integrationLimitsLimber[totalObs][1], chi.second))-0.5));
                                totalObs++;
                            }
                        }
                        else {  // if calculate full covar
                            LimberSpectra[totalObs] = powerSpectrum(matterSpectrum, corrFuncs[obsIndices[obsType][0]+g1], corrFuncs[obsIndices[obsType][1]+g2], k, a_size, k_size);
                            integrationLimitsLimber[totalObs][0] = std::max(log(1./(std::min(paramFile.galZHigh[g1], paramFile.galZHigh[g2])+1)), aLow);
                            integrationLimitsLimber[totalObs][1] = std::min(log(1./(std::max(paramFile.galZLow[g1], paramFile.galZLow[g2])+1)), aHigh);
                            if(integrationLimitsLimber[totalObs][0] > integrationLimitsLimber[totalObs][1]) { // non-overlapping redshift bins
                                integrationLimitsLimber[totalObs][0] = integrationLimitsLimber[totalObs][1];
                            }
                            // Set ell-limit by theoretical integral lower bound
                            double lowestZ = std::min(log(1./(std::min(paramFile.galZLow[g1], paramFile.galZLow[g2])+1)), aHigh);
                            l_lims[totalObs] = int(round((K_LIMIT*gsl_spline_eval(chi.first, lowestZ, chi.second))-0.5));
                            totalObs++;
                        }
                    }
                }
            }
            // Lensing auto-spectra
            else if(kernelInstruments[obsType]==6) {
                for(int g1=0; g1<paramFile.numLens; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numLens; g2++) {  // column galaxies
                        if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                            if(g1==g2) {
                                LimberSpectra[totalObs] = powerSpectrum(matterSpectrum, corrFuncs[obsIndices[obsType][0]+g1], corrFuncs[obsIndices[obsType][1]+g1], k, a_size, k_size);
                                integrationLimitsLimber[totalObs][0] = std::max(log(1./(paramFile.lensZHigh[g1]+1)), aLow);
                                integrationLimitsLimber[totalObs][1] = aHigh;
                                l_lims[totalObs] = int(round((K_LIMIT_SHEAR*gsl_spline_eval(chi.first, integrationLimitsLimber[totalObs][1], chi.second))-0.5));
                                totalObs++;
                            }
                        }
                        else {  // if calculate full covar
                            LimberSpectra[totalObs] = powerSpectrum(matterSpectrum, corrFuncs[obsIndices[obsType][0]+g1], corrFuncs[obsIndices[obsType][1]+g2], k, a_size, k_size);
                            integrationLimitsLimber[totalObs][0] = std::max(log(1./(std::min(paramFile.lensZHigh[g1], paramFile.lensZHigh[g2])+1)), aLow);
                            integrationLimitsLimber[totalObs][1] = aHigh;
                            l_lims[totalObs] = int(round((K_LIMIT_SHEAR*gsl_spline_eval(chi.first, integrationLimitsLimber[totalObs][1], chi.second))-0.5));
                            totalObs++;
                        }
                    }
                }
            }
            // Lensing-Galaxy cross-spectra
            else if(kernelInstruments[obsType]==4) {
                for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numLens; g2++) {  // column lensing
                        LimberSpectra[totalObs] = powerSpectrum(matterSpectrum, corrFuncs[obsIndices[obsType][0]+g1], corrFuncs[obsIndices[obsType][1]+g2], k, a_size, k_size);
                        integrationLimitsLimber[totalObs][0] = std::max(log(1./(std::min(paramFile.galZHigh[g1], paramFile.lensZHigh[g2])+1)), aLow);
                        integrationLimitsLimber[totalObs][1] = std::min(log(1./(paramFile.galZLow[g1]+1)), aHigh);
                        if(integrationLimitsLimber[totalObs][0] > integrationLimitsLimber[totalObs][1]) { // Source galaxies are in front of lensing galaxies
                            integrationLimitsLimber[totalObs][0] = integrationLimitsLimber[totalObs][1];
                        }
                        if(paramFile.galZHigh[g1] > paramFile.lensZHigh[g2]) { // Lens galaxy bin extends beyond source galaxy bin
                            integrationLimitsLimber[totalObs][0] = integrationLimitsLimber[totalObs][1];
                        }
                        l_lims[totalObs] = int(round((K_LIMIT*gsl_spline_eval(chi.first, integrationLimitsLimber[totalObs][1], chi.second))-0.5));
                        totalObs++;
                    }
                }
            }
            // Galaxy cross-spectra (no lensing)
            else if(kernelInstruments[obsType]==1) {
                for(int g=0; g<paramFile.numGal; g++) {
                    LimberSpectra[totalObs] = powerSpectrum(matterSpectrum, corrFuncs[obsIndices[obsType][0]], corrFuncs[obsIndices[obsType][1]+g], k, a_size, k_size);
                    integrationLimitsLimber[totalObs][0] = std::max(log(1./(paramFile.galZHigh[g]+1)), aLow);
                    integrationLimitsLimber[totalObs][1] = std::min(log(1./(paramFile.galZLow[g]+1)), aHigh);
                    l_lims[totalObs] = int(round((K_LIMIT*gsl_spline_eval(chi.first, integrationLimitsLimber[totalObs][1], chi.second))-0.5));
                    totalObs++;
                }
            }
            // Lensing cross-spectra (no galaxy)
            else if(kernelInstruments[obsType]==3) {
                for(int g=0; g<paramFile.numLens; g++) {
                    LimberSpectra[totalObs] = powerSpectrum(matterSpectrum, corrFuncs[obsIndices[obsType][0]], corrFuncs[obsIndices[obsType][1]+g], k, a_size, k_size);
                    integrationLimitsLimber[totalObs][0] = std::max(log(1./(paramFile.lensZHigh[g]+1)), aLow);
                    integrationLimitsLimber[totalObs][1] = aHigh;
                    l_lims[totalObs] = int(round((K_LIMIT_SHEAR*gsl_spline_eval(chi.first, integrationLimitsLimber[totalObs][1], chi.second))-0.5));
                    totalObs++;
                }
            }
            // Non-galaxy, non-lensing spectra
            else {
                LimberSpectra[totalObs] = powerSpectrum(matterSpectrum, corrFuncs[obsIndices[obsType][0]], corrFuncs[obsIndices[obsType][1]], k, a_size, k_size);
                integrationLimitsLimber[totalObs][0] = aLow;
                integrationLimitsLimber[totalObs][1] = aHigh;
                l_lims[totalObs] = int(round((K_LIMIT_SHEAR*gsl_spline_eval(chi.first, integrationLimitsLimber[totalObs][1], chi.second))-0.5));
                totalObs++;
            }
        }
    }
    for(int i = 0; i<2+(2*paramFile.numGal)+(2*paramFile.numLens); i++) {  // Memory clear-up of correlation functions
        for(int j = 0; j<k_size; j++) {
            delete[] corrFuncs[i][j];
        }
        delete[] corrFuncs[i];
    }
    delete[] corrFuncs;

    // Setting dataset lengths by mean result
    if(suffix=="Mean") {
        // Set l_lim in struct
        int* ellLimits = new int[totalObs];
        for(int i=0; i<totalObs; i++) {
            l_lims[i] = std::min(l_lims[i], l_max2+1);
            ellLimits[i] = l_lims[i];
        }
        ellLims = ellLimits;
    }
    else {
        // Retrieve l_lim from struct
        for(int i=0; i<totalObs; i++) {
            l_lims[i] = ellLims[i];
        }
    }
    
    // Non-linear correlation functions
    double*** corrNLFuncs = new double**[2+(2*paramFile.numGal)+(2*paramFile.numLens)];
    for (int i = 0; i<2+(2*paramFile.numGal)+(2*paramFile.numLens); i++) {
        corrNLFuncs[i] = new double*[k_size];
        for(int j = 0; j<k_size; j++) {
            corrNLFuncs[i][j] = new double[a_size];
        }
    }
    ISWPowerSpectrumCorr(corrNLFuncs[0], F_of_k, HubbleRate, NLISWT, a_size, k_size);
    CMBLensePowerSpectrumCorr(corrNLFuncs[1], CMBLensingKernel, NLCMBLenseT, a_size, k_size);
    for(int g = 0; g<paramFile.numGal; g++) {
        galPowerSpectrumCorr(corrNLFuncs[2+g], bias[g], galSelectionFunc[g], NLGalT, a_size, k_size);
        GalLensePowerSpectrumCorr(corrNLFuncs[2+paramFile.numGal+g], galNumDen[g], bias[g], NLGalT, a_size, k_size);
    }
    for(int g = 0; g<paramFile.numLens; g++) {
        LensePowerSpectrumCorr(corrNLFuncs[2+(2*paramFile.numGal)+g], lensingKernel[g], NLGalT, F_of_k, k, a_size, k_size);
        IntAlignPowerSpectrumCorr(corrNLFuncs[2+(2*paramFile.numGal)+paramFile.numLens+g], IAKernel[g], lensNumDen[g], 
                                    NLGalT, a_size, k_size);
    }

    // Limber splines, non-linear spectra
    std::pair<gsl_spline *, gsl_interp_accel *>** NLLimberSpectra = new std::pair<gsl_spline *, gsl_interp_accel *>*[totalObs];
    totalObs = 0;
    for(int obsType=0; obsType<paramFile.numObs; obsType++) {  // Generate Limber power spectra
        if(paramFile.calcObs[obsType]) {
            // Galaxy auto-spectra
            if(kernelInstruments[obsType]==2) {
                for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numGal; g2++) {  // column galaxies
                        if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                            if(g1==g2) {
                                NLLimberSpectra[totalObs] = powerSpectrum(NLCorrectionT, corrNLFuncs[obsIndices[obsType][0]+g1], corrNLFuncs[obsIndices[obsType][1]+g1], k, a_size, k_size);
                                totalObs++;
                            }
                        }
                        else {  // if calculate full covar
                            NLLimberSpectra[totalObs] = powerSpectrum(NLCorrectionT, corrNLFuncs[obsIndices[obsType][0]+g1], corrNLFuncs[obsIndices[obsType][1]+g2], k, a_size, k_size);
                            totalObs++;
                        }
                    }
                }
            }
            // Lensing auto-spectra
            else if(kernelInstruments[obsType]==6) {
                for(int g1=0; g1<paramFile.numLens; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numLens; g2++) {  // column galaxies
                        if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                            if(g1==g2) {
                                NLLimberSpectra[totalObs] = powerSpectrum(NLCorrectionT, corrNLFuncs[obsIndices[obsType][0]+g1], corrNLFuncs[obsIndices[obsType][1]+g1], k, a_size, k_size);
                                totalObs++;
                            }
                        }
                        else {  // if calculate full covar
                            NLLimberSpectra[totalObs] = powerSpectrum(NLCorrectionT, corrNLFuncs[obsIndices[obsType][0]+g1], corrNLFuncs[obsIndices[obsType][1]+g2], k, a_size, k_size);
                            totalObs++;
                        }
                    }
                }
            }
            // Lensing-Galaxy cross-spectra
            else if(kernelInstruments[obsType]==4) {
                for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numLens; g2++) {  // column lensing
                        NLLimberSpectra[totalObs] = powerSpectrum(NLCorrectionT, corrNLFuncs[obsIndices[obsType][0]+g1], corrNLFuncs[obsIndices[obsType][1]+g2], k, a_size, k_size);
                        totalObs++;
                    }
                }
            }
            // Galaxy cross-spectra (no lensing)
            else if(kernelInstruments[obsType]==1) {
                for(int g=0; g<paramFile.numGal; g++) {
                    NLLimberSpectra[totalObs] = powerSpectrum(NLCorrectionT, corrNLFuncs[obsIndices[obsType][0]], corrNLFuncs[obsIndices[obsType][1]+g], k, a_size, k_size);
                    totalObs++;
                }
            }
            // Lensing cross-spectra (no galaxy)
            else if(kernelInstruments[obsType]==3) {
                for(int g=0; g<paramFile.numLens; g++) {
                    NLLimberSpectra[totalObs] = powerSpectrum(NLCorrectionT, corrNLFuncs[obsIndices[obsType][0]], corrNLFuncs[obsIndices[obsType][1]+g], k, a_size, k_size);
                    totalObs++;
                }
            }
            // Non-galaxy, non-lensing spectra
            else {
                NLLimberSpectra[totalObs] = powerSpectrum(NLCorrectionT, corrNLFuncs[obsIndices[obsType][0]], corrNLFuncs[obsIndices[obsType][1]], k, a_size, k_size);
                totalObs++;
            }
        }
    }
    for(int i = 0; i<2+(2*paramFile.numGal)+(2*paramFile.numLens); i++) {  // Memory clear-up of correlation functions
        for(int j = 0; j<k_size; j++) {
            delete[] corrNLFuncs[i][j];
        }
        delete[] corrNLFuncs[i];
    }
    delete[] corrNLFuncs;




    if(chiLow<chiData[1]) {
        chiLow=chiData[1];
    }
    double* shotConstant = new double[paramFile.numGal];
    for (int g1 = 0; g1<paramFile.numGal; g1++) {
        double numDen = 0;
        if(suffix=="numDen"+std::to_string(g1)+"+") {
            numDen = paramFile.galNumDen[g1]+paramFile.galNumDenDev[g1];
        }
        else if(suffix=="numDen"+std::to_string(g1)+"-") {
            numDen = paramFile.galNumDen[g1]-paramFile.galNumDenDev[g1];
        }
        else {
            numDen = paramFile.galNumDen[g1];
        }
        double meanNumDen = paramFile.galNumDen[g1];
        double meanAngNumDen = meanNumDen*(180*180)/(PI*PI);
        double angNumDen = numDen*(180*180)/(PI*PI);
        shotConstant[g1] = (1./angNumDen)-(1./meanAngNumDen);
        std::cout<<"Galaxy bins "<<g1+1<<" shot noise: "<<shotConstant[g1]<<std::endl;
    }
    delete[] invertedChiData;

    // give chi limits, and k limits for inspection
    std::cout << gsl_spline_eval (chi.first, loga[0], chi.second) << ", " << gsl_spline_eval (chi.first, loga[a_size - 2], chi.second) << std::endl;
    std::cout << k[0] << ", " << k[k_size - 1] << std::endl;

    // Full sky integration kernels with limits
    double*** densityKernels = new double**[2+(2*paramFile.numGal)+(2*paramFile.numLens)];
    densityKernels[0] = ISWKernel(densityISWT, F_of_k, HubbleRate, k_size, a_size);
    densityKernels[1] = CMBLenseKernel(densityCMBLenseT, CMBLensingKernel, k_size, a_size);
    for(int g = 0; g<paramFile.numGal; g++) {
        densityKernels[2+g] = GalKernel(densityGalT, bias[g], galSelectionFunc[g], k_size, a_size);
        densityKernels[2+paramFile.numGal+g] = GalLenseKernel(galNumDen[g], densityGalT, bias[g], k_size, a_size);
    }
    for(int g = 0; g<paramFile.numLens; g++) {
        densityKernels[2+(2*paramFile.numGal)+g] = LenseKernel(densityGalT, F_of_k, k, lensingKernel[g], k_size, a_size);
        densityKernels[2+(2*paramFile.numGal)+paramFile.numLens+g] = IntAlignKernel(densityGalT, IAKernel[g], lensNumDen[g], k_size, a_size);
    }


    double** AngularPowerSpectra = new double*[totalObs];
    for(int i = 0; i<totalObs; i++) {
        AngularPowerSpectra[i] = new double[l_max2+1];
    }
    double** NLAngularPowerSpectra = new double*[totalObs];
    for(int i = 0; i<totalObs; i++) {
        NLAngularPowerSpectra[i] = new double[l_max2+1];
    }

    // Calculate integrals for all l values
    std::cout << "Beginning C_ell computation. l_max = " << l_max2 << "." << std::endl;
    std::thread* Threads = new std::thread[numThreads];
    for(int i = 0; i<numThreads; i++) {
        Threads[i] = std::thread(&angularPowerSpectrumCurvedSky, AngularPowerSpectra, densityKernels,
                                    paramFile, P_primordial, kernels, kernelInstruments, obsIndices, k, kChi, chiData,
                                    chiLow, chiHigh, lensingSourceChi, galaxyBinChi, k_size, a_size, l_max1, numThreads, i,
                                    BesselArray_size);
    }
    for(int i = 0; i<numThreads; i++) {
        Threads[i].join();
    }
    for(int j = 0; j<totalObs; j++) {
        for(int i = 0; i<numThreads; i++) {
            Threads[i] = std::thread(&threadFunctionLimber, AngularPowerSpectra[j], LimberSpectra[j], scale, loga,
                                HubbleRate, a_size, chi, l_lims[j], l_max1+1, l_max2, numThreads, i, k_size,
                                integrationLimitsLimber[j][0], integrationLimitsLimber[j][1], k[k_size-1]);
        }
        for(int i = 0; i<numThreads; i++) {
            Threads[i].join();
        }
    }
    for(int j = 0; j<totalObs; j++) {
        for(int i = 0; i<numThreads; i++) {
            Threads[i] = std::thread(&threadFunctionLimber, NLAngularPowerSpectra[j], NLLimberSpectra[j], scale, loga,
                                HubbleRate, a_size, chi, l_lims[j], 0, l_max2, numThreads, i, k_size,
                                integrationLimitsLimber[j][0], integrationLimitsLimber[j][1], k[k_size-1]);
        }
        for(int i = 0; i<numThreads; i++) {
            Threads[i].join();
        }
    }
    delete[] Threads;
    for(int j = 0; j<totalObs; j++) {
        for(int i = 0; i < l_max2 + 1; i++) {
            AngularPowerSpectra[j][i] = AngularPowerSpectra[j][i] + NLAngularPowerSpectra[j][i];
        }
    }
    std::cout<<"Calculation complete."<<std::endl;

    // Generate ell array
    int * ell = new int[l_max2 + 1];
    for (int i = 0; i < l_max2 + 1; i++) {
        ell[i] = i;
    }
    // Output spectra
    totalObs = 0;
    for(int obsType=0; obsType<paramFile.numObs; obsType++) {
        if(paramFile.calcObs[obsType]) {
            // Galaxy auto-spectra
            if(kernelInstruments[obsType]==2) {
                for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numGal; g2++) {  // column galaxies
                        if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                            if(g1==g2) {
                                for(int i = 0; i < l_lims[totalObs]; i++) {
                                    AngularPowerSpectra[totalObs][i] = AngularPowerSpectra[totalObs][i] + shotConstant[g1];
                                    NLAngularPowerSpectra[totalObs][i] = NLAngularPowerSpectra[totalObs][i];
                                }
                                outputIntData(outputDirectory + "C_ell_"+paramFile.obsName[obsType]+paramFile.galSuffix[g1]+paramFile.galSuffix[g2]+suffix+outFileSuffix+".dat", ell, AngularPowerSpectra[totalObs], l_max2 + 1);
                                outputIntData(outputDirectory + "C_ell_NL"+paramFile.obsName[obsType]+paramFile.galSuffix[g1]+paramFile.galSuffix[g2]+suffix+outFileSuffix+".dat", ell, NLAngularPowerSpectra[totalObs], l_max2 + 1);
                                totalObs++;
                            }
                        }
                        else {  // if calculate full covar
                            if(g1==g2) {
                                for(int i = 0; i < l_lims[totalObs]; i++) {
                                    AngularPowerSpectra[totalObs][i] = AngularPowerSpectra[totalObs][i] + shotConstant[g1];
                                    NLAngularPowerSpectra[totalObs][i] = NLAngularPowerSpectra[totalObs][i];
                                }
                            }
                            if(integrationLimitsLimber[totalObs][0] == integrationLimitsLimber[totalObs][1]) { // non-overlapping redshift bins
                                for(int i = 0; i < l_max2 + 1; i++) {
                                    AngularPowerSpectra[totalObs][i] = 0;
                                }
                            }
                            outputIntData(outputDirectory + "C_ell_"+paramFile.obsName[obsType]+paramFile.galSuffix[g1]+paramFile.galSuffix[g2]+suffix+outFileSuffix+".dat", ell, AngularPowerSpectra[totalObs], l_max2 + 1);
                            outputIntData(outputDirectory + "C_ell_NL"+paramFile.obsName[obsType]+paramFile.galSuffix[g1]+paramFile.galSuffix[g2]+suffix+outFileSuffix+".dat", ell, NLAngularPowerSpectra[totalObs], l_max2 + 1);
                            totalObs++;
                        }
                    }
                }
            }
            // Lensing auto-spectra
            else if(kernelInstruments[obsType]==6) {
                for(int g1=0; g1<paramFile.numLens; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numLens; g2++) {  // column galaxies
                        if(paramFile.calcCrossCov==false) {  // if just calculate diagonal elements of covar
                            if(g1==g2) {
                                outputIntData(outputDirectory + "C_ell_"+paramFile.obsName[obsType]+paramFile.lensSuffix[g1]+paramFile.lensSuffix[g2]+suffix+outFileSuffix+".dat", ell, AngularPowerSpectra[totalObs], l_max2 + 1);
                                outputIntData(outputDirectory + "C_ell_NL"+paramFile.obsName[obsType]+paramFile.lensSuffix[g1]+paramFile.lensSuffix[g2]+suffix+outFileSuffix+".dat", ell, NLAngularPowerSpectra[totalObs], l_max2 + 1);
                                totalObs++;
                            }
                        }
                        else {  // if calculate full covar
                            outputIntData(outputDirectory + "C_ell_"+paramFile.obsName[obsType]+paramFile.lensSuffix[g1]+paramFile.lensSuffix[g2]+suffix+outFileSuffix+".dat", ell, AngularPowerSpectra[totalObs], l_max2 + 1);
                            outputIntData(outputDirectory + "C_ell_NL"+paramFile.obsName[obsType]+paramFile.lensSuffix[g1]+paramFile.lensSuffix[g2]+suffix+outFileSuffix+".dat", ell, NLAngularPowerSpectra[totalObs], l_max2 + 1);
                            totalObs++;
                        }
                    }
                }
            }
            // Lensing-Galaxy cross-spectra
            else if(kernelInstruments[obsType]==4) {
                for(int g1=0; g1<paramFile.numGal; g1++) {  // row galaxies
                    for(int g2=0; g2<paramFile.numLens; g2++) {  // column galaxies
                        if(integrationLimitsLimber[totalObs][0] == integrationLimitsLimber[totalObs][1]) { // non-overlapping redshift bins
                            for(int i = 0; i < l_max2 + 1; i++) {
                                AngularPowerSpectra[totalObs][i] = 0;
                            }
                        }
                        outputIntData(outputDirectory + "C_ell_"+paramFile.obsName[obsType]+paramFile.galSuffix[g1]+paramFile.lensSuffix[g2]+suffix+outFileSuffix+".dat", ell, AngularPowerSpectra[totalObs], l_max2 + 1);
                        outputIntData(outputDirectory + "C_ell_NL"+paramFile.obsName[obsType]+paramFile.galSuffix[g1]+paramFile.lensSuffix[g2]+suffix+outFileSuffix+".dat", ell, NLAngularPowerSpectra[totalObs], l_max2 + 1);
                        totalObs++;
                    }
                }
            }
            // Galaxy cross-spectra
            else if(kernelInstruments[obsType]==1) {
                for(int g=0; g<paramFile.numGal; g++) {
                    outputIntData(outputDirectory + "C_ell_"+paramFile.obsName[obsType]+paramFile.galSuffix[g]+suffix+outFileSuffix+".dat", ell, AngularPowerSpectra[totalObs], l_max2 + 1);
                    outputIntData(outputDirectory + "C_ell_NL"+paramFile.obsName[obsType]+paramFile.galSuffix[g]+suffix+outFileSuffix+".dat", ell, NLAngularPowerSpectra[totalObs], l_max2 + 1);
                    totalObs++;
                }
            }
            // Lensing cross-spectra
            else if(kernelInstruments[obsType]==3) {
                for(int g=0; g<paramFile.numLens; g++) {
                    outputIntData(outputDirectory + "C_ell_"+paramFile.obsName[obsType]+paramFile.lensSuffix[g]+suffix+outFileSuffix+".dat", ell, AngularPowerSpectra[totalObs], l_max2 + 1);
                    outputIntData(outputDirectory + "C_ell_NL"+paramFile.obsName[obsType]+paramFile.lensSuffix[g]+suffix+outFileSuffix+".dat", ell, NLAngularPowerSpectra[totalObs], l_max2 + 1);
                    totalObs++;
                }
            }
            // Non-galaxy spectra
            else {
                outputIntData(outputDirectory + "C_ell_"+paramFile.obsName[obsType]+suffix+outFileSuffix+".dat", ell, AngularPowerSpectra[totalObs], l_max2 + 1);
                outputIntData(outputDirectory + "C_ell_NL"+paramFile.obsName[obsType]+suffix+outFileSuffix+".dat", ell, NLAngularPowerSpectra[totalObs], l_max2 + 1);
                totalObs++;
            }
        }
    }
    std::cout<<"Files generated."<<std::endl;

    // Memory clean-up
    gsl_spline_free (chi.first);
    gsl_interp_accel_free (chi.second);
    delete[] scale;
    delete[] chiData;
    delete[] loga;
    delete[] a2Hinv;
    delete[] HubbleRate;
    delete[] a2H;
    delete[] k;
    delete[] P_primordial;
    delete[] F_of_k;
    delete[] kernels;
    delete[] kernelInstruments;
    delete[] CMBLensingKernel;
    delete[] lensingSourceChi;
    delete[] ell;
    for (int i = 0; i<k_size; i++) {
        delete[] densityISWT[i];
        delete[] densityGalT[i];
        delete[] densityCMBLenseT[i];
        delete[] matterSpectrum[i];
        delete[] NLISWT[i];
        delete[] NLGalT[i];
        delete[] NLCMBLenseT[i];
        delete[] NLCorrectionT[i];
    }
    delete[] densityISWT;
    delete[] densityGalT;
    delete[] densityCMBLenseT;
    delete[] matterSpectrum;
    delete[] NLISWT;
    delete[] NLGalT;
    delete[] NLCMBLenseT;
    delete[] NLCorrectionT;
    for(int g = 0; g<paramFile.numGal; g++) {
        delete[] galSelectionFunc[g];
        delete[] galNumDen[g];
        delete[] galaxyBinChi[g];
        delete[] bias[g];
    }
    for(int g = 0; g<paramFile.numLens; g++) {
        delete[] lensNumDen[g];
        delete[] lensingKernel[g];
        delete[] IAKernel[g];
    }
    delete[] galSelectionFunc;
    delete[] galNumDen;
    delete[] normalisation;
    delete[] lensNumDen;
    delete[] galaxyBinChi;
    delete[] l_lims;
    delete[] bias;
    delete[] IAKernel;
    delete[] lensingKernel;
    delete[] shotConstant;
    for(int i = 0; i<2+(2*paramFile.numGal)+(2*paramFile.numLens); i++) {
        for(int j = 0; j<k_size; j++) {
            delete[] densityKernels[i][j];
        }
        delete[] densityKernels[i];
    }
    delete[] densityKernels;
    for(int i = 0; i<totalObs; i++) {
        delete[] AngularPowerSpectra[i];
        delete[] NLAngularPowerSpectra[i];
        for(int j = 0; j<a_size; j++) {
            gsl_spline_free (LimberSpectra[i][j].first);
            gsl_interp_accel_free (LimberSpectra[i][j].second);
            gsl_spline_free (NLLimberSpectra[i][j].first);
            gsl_interp_accel_free (NLLimberSpectra[i][j].second);
        }
        delete[] LimberSpectra[i];
        delete[] NLLimberSpectra[i];
        delete[] integrationLimitsLimber[i];
    }
    delete[] AngularPowerSpectra;
    delete[] LimberSpectra;
    delete[] NLAngularPowerSpectra;
    delete[] NLLimberSpectra;
    delete[] integrationLimitsLimber;
    for(int i = 0; i<paramFile.numObs; i++) {
        delete[] obsIndices[i];
    }
    delete[] obsIndices;

    return;
}

void forecastSpectraODE(parameterFile paramFile, int l_max1, int l_max2, int numThreads, double aLow, double aHigh, 
                        std::string inputSuffix, std::string suffix) {
    
    int numParams = paramFile.numParams;
    int* l_lims = new int[1]();

    int BesselArray_size = kChiLength();
    double* kChi = kChiArray(BesselArray_size);
    std::cout<<kChi[BesselArray_size-1]<<std::endl;

    // Always compute mean spectra in order to set l_lims
    std::cout<<"Getting mean spectra."<<std::endl;
    getSpectraODE(paramFile, l_lims, l_max1, l_max2, numThreads, aLow, aHigh, BesselArray_size, kChi, inputSuffix, "Mean", suffix);

    for(int i = 0; i<numParams; i++) {
        if(paramFile.calcParam[i]) {
            std::cout<<"Getting "<<paramFile.paramNames[i]<<" spectra."<<std::endl;
            getSpectraODE(paramFile, l_lims, l_max1, l_max2, numThreads, aLow, aHigh, BesselArray_size, kChi, inputSuffix, paramFile.paramNames[i]+"+", suffix);
            getSpectraODE(paramFile, l_lims, l_max1, l_max2, numThreads, aLow, aHigh, BesselArray_size, kChi, inputSuffix, paramFile.paramNames[i]+"-", suffix);
        }
    }
    // Memory clean-up
    delete[] paramFile.paramNames;
    delete[] paramFile.calcParam;
    delete[] paramFile.galFSky;
    delete[] paramFile.galNumDen;
    delete[] paramFile.galNumDenDev;
    delete[] paramFile.galZLow;
    delete[] paramFile.galZHigh;
    delete[] paramFile.galZMed;
    delete[] paramFile.galSuffix;
    delete[] paramFile.lensFSky;
    delete[] paramFile.lensNumden;
    delete[] paramFile.lensNumDenDev;
    delete[] paramFile.lensZLow;
    delete[] paramFile.lensZHigh;
    delete[] paramFile.lensZMed;
    delete[] paramFile.lensSuffix;
    delete[] paramFile.kernInst;
    delete[] paramFile.obsName;
    delete[] paramFile.calcObs;
    for(int i = 0; i<paramFile.numObs; i++) {
        delete[] paramFile.calcKernel[i];
    }
    delete[] paramFile.calcKernel;
    delete[] l_lims;
    delete[] kChi;
    return;
}


//////////
// MAIN //
//////////

int main(int argc, char **argv) {
    // To compile, use
    // g++ -std=c++11 -lgsl -o PowerSpec -pthread -O3 main.cpp stdafx.cpp
    // To execute for exact sltn from ell=0 to ell=30 and limber from ell=31 to ell=100, from z=0 to z=21,
    // with no suffix
    // ./PowerSpec ell1 ell2 #threads zlow zhigh /Path/To/parameterFile.csv inputSuffix, suffix

    int l_max1, l_max2;
    int numThreads;
    double aLow, aHigh;
    std::string paramFile;
    std::string inputSuffix;
    std::string suffix;

    if (argc!=9) {
        std::cout<<"Arguments are missing."<<std::endl;
        return 0;
    }
    else {
        l_max1 = std::atoi(argv[1]);
        l_max2 = std::atoi(argv[2]);
        numThreads = std::atoi(argv[3]);
        aHigh = log(1./(1+std::stod(argv[4])));
        aLow = log(1./(1+std::stod(argv[5])));
        paramFile = argv[6];
        inputSuffix = argv[7];
        if(inputSuffix=="None" || inputSuffix=="NONE") {
            inputSuffix = "";
        }
        suffix = argv[8];
        if(suffix=="None" || suffix=="NONE") {
            suffix = "";
        }
    }
    forecastSpectraODE(readParamFile(paramFile), l_max1, l_max2, numThreads, aLow, aHigh, inputSuffix, suffix);
    
    return 0;
}
