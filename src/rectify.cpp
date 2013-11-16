/*
 * Rectify: Converge towards a target image using a Monte Carlo algorithm.
 *          Run "rectify -h" for help.
 *
 * Copyright (C) 2013  Lester Hedges <lester.hedges+rectify@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include "lodepng.h"
#include "MersenneTwister.h"

using namespace std;

// FUNCTION PROTOTYPES
void printHelpMessage();
void parseCommandLineArguments(int, char**, ostringstream&, ostringstream&,
    ostringstream&, long&, unsigned int&, double&, double&, bool&, bool&, bool&);
void generateSamplePoints(vector <long>&, long&, unsigned int&, bool&);
void decodeImage(ostringstream&, vector <unsigned char>&, unsigned int&, unsigned int&);
void encodeImage(unsigned int&, ostringstream&, vector <unsigned char>&, unsigned int&, unsigned int&);
void applyTrialBrushStroke(vector <unsigned char>&, unsigned int&,
    unsigned int&, MTRand&, unsigned int&, int&, int&, int&, int&);
void blendBlock(vector <unsigned char>&, unsigned int&, unsigned int&, double, int&, int&, int&, int&);
void copyBrushStroke(vector <unsigned char>&, vector <unsigned char>&, unsigned int&, int&, int&, int&, int&);
double computeBlockValue(vector <unsigned char>&, vector <unsigned char>&, unsigned int&, int&, int&, int&, int&);
void initializeImage(vector <unsigned char>&, unsigned int&, unsigned int&);
void convertToGrayscale(vector <unsigned char>&, unsigned int&, unsigned int&);
void resetTallyCounter(vector <vector <long> >&, unsigned int&);

// BEGIN MAIN FUNCTION
int main (int argc, char** argv)
{
    // counters
    long i;                             // trial iteration counter
    long nAccepted = 0;                 // number of accepted trials
    unsigned int samples = 0;           // number of samples
    unsigned int j,k;                   // extra variables for inner loops

    // file name buffers
    ostringstream targetFileName;       // file name of target image
    ostringstream startingFileName;     // file name of starting image

    // defaults
    long iterations = 1000000;          // total number of iterations
    unsigned int frames = 1000;         // number of images to sample
    double temperature = 0.1;           // temperature of thermal bath
    double maxStrokeFraction = 0.1;     // maximum extent of stroke (in fraction of target image width)
    bool isLogarithmic = false;         // whether sampling is performed at logarithmic intervals
    bool isCanvas = false;              // whether a starting canvas is specified
    bool isMonitor = false;             // whether to monitor acceptance statistics

    // absoulte pixel value of image block before and after trial
    double currentBlockValue, trialBlockValue, blockValueChange;

    // bounday values of stroke block
    int x1, x2, y1, y2;

    // width and height of image in pixels
    unsigned int width, height;

    // distance from target image
    double error;

    // random number generator
    MTRand rng;

    // output directory name (none selected if empty)
    ostringstream directoryName;
    directoryName.str("");

    // parse command-line arguments
    parseCommandLineArguments(argc, argv, targetFileName, startingFileName, directoryName,
        iterations, frames, maxStrokeFraction, temperature, isLogarithmic, isCanvas, isMonitor);

    // file name variables for stat logging
    ostringstream fileName, fileName2, fileName3;

    // file pointer
    FILE *pFile;

    if (isMonitor)
    {
        fileName.str("");
        fileName2.str("");
        fileName3.str("");

        if (strcmp(directoryName.str().c_str(), "") == 0)
        {
            fileName << "log.txt";
            fileName2 << "hist.txt";
            fileName3 << "hist_xy.txt";
        }
        else
        {
            fileName << directoryName.str() << "/log.txt";
            fileName2 << directoryName.str() << "/hist.txt";
            fileName3 << directoryName.str() << "/hist_xy.txt";
        }

        // wipe existing log and histogram files
        pFile = fopen(fileName.str().c_str(), "w");
        fclose(pFile);
        pFile = fopen(fileName2.str().c_str(), "w");
        fclose(pFile);
        pFile = fopen(fileName3.str().c_str(), "w");
        fclose(pFile);
    }

    // array of sampling check points
    vector <long> samplePoints;

    // generate sampling points
    generateSamplePoints(samplePoints, iterations, frames, isLogarithmic);

    // set inverse temperature
    double beta = 1.0/temperature;

    // initialize and read target image from file
    vector <unsigned char> targetImage, trialImage, currentImage;

    // decode target image from file
    decodeImage(targetFileName, targetImage, width, height);

    // set maximum extent of stroke (in pixels)
    unsigned int maxStrokeWidth = maxStrokeFraction*width;

    // set total number of pixels in canvas
    unsigned int pixels = width*height;

    // tally arrays for accepted stroke widths and heights
    vector <vector <long> > strokeTally(maxStrokeWidth, vector <long> (maxStrokeWidth));

    // clear tally counter
    resetTallyCounter(strokeTally, maxStrokeWidth);

    // acceptance flag
    bool isAccepted;

    // canvas image is specified
    if (isCanvas)
    {
        unsigned int w, h;

        // decode canvas image from file
        decodeImage(startingFileName, currentImage, w, h);

        bool wrongSize = false;
        if (w != width) wrongSize = true;
        else if (h != height) wrongSize = true;
        if (wrongSize)
        {
            cout << "Warning: starting image does not have same dimensions as target. Using blank canvas!" << endl;
            currentImage.resize(targetImage.size());
            initializeImage(currentImage, width, height);
        }

        if (currentImage[0] != currentImage[1])
        {
            if (currentImage[0] != currentImage[2])
            {
                cout << "Starting image is not grayscale... converting." << endl;
                convertToGrayscale(currentImage, width, height);
            }
        }
    }
    else
    {
        currentImage.resize(targetImage.size());
        initializeImage(currentImage, width, height);
    }

    trialImage = currentImage;

    // set canvas boundaries
    x1 = 0; x2 = width-1; y1 = 0; y2 = width-1;

    // work out current distance from target
    error = computeBlockValue(currentImage, targetImage, width, x1, y1, x2, y2);

    cout << "Starting image generation..." << endl;

    // MAIN LOOP
    for (i=1;i<=iterations;i++)
    {
        // reset acceptance flag
        isAccepted = false;

        // make trial brush stroke
        applyTrialBrushStroke(trialImage, width, height, rng, maxStrokeWidth, x1, y1, x2, y2);

        // compute difference in block values between target and evolved image
        currentBlockValue = computeBlockValue(currentImage, targetImage, width, x1, y1, x2, y2);

        // compute difference in block values between target and trial image
        trialBlockValue = computeBlockValue(trialImage, targetImage, width, x1, y1, x2, y2);

        // evaluate change in block value
        blockValueChange = trialBlockValue - currentBlockValue;

        // check whether move is accepted
        if (blockValueChange < 0) isAccepted = true;
        else
        {
            if (rng() < exp(-beta*blockValueChange)) isAccepted = true;
        }

        // copy block from trial to current image if accepted, do the opposite if rejected
        if (isAccepted)
        {
            nAccepted++;
            copyBrushStroke(trialImage, currentImage, width, x1, y1, x2, y2);
            error += blockValueChange;

            // update tally counter (max sure block extrema are different)
            if (isMonitor)
            {
                if (x2 != x1)
                {
                    if (y2 != y1) strokeTally[x2-x1-1][y2-y1-1]++;
                }
            }
        }
        else copyBrushStroke(currentImage, trialImage, width, x1, y1, x2, y2);

        // encode current image and write stats to stdout and file
        if (i >= samplePoints[samples])
        {
            samples++;

            // write to stdout
            printf("Frame: %4d, error: %5.4f, acceptance: %5.4f\n", samples, error/pixels, ((double) nAccepted/i));

            // write acceptance statistics to file
            if (isMonitor)
            {
                // write to log file
                pFile = fopen(fileName.str().c_str(), "a");
                fprintf(pFile, "%d %5.4f %5.4f\n", samples, error/pixels, ((double) nAccepted/i));
                fclose(pFile);

                // write to histogram data file
                pFile = fopen(fileName2.str().c_str(), "w");
                for (j=0;j<maxStrokeWidth;j++)
                {
                    for (k=0;k<maxStrokeWidth;k++)
                    {
                        fprintf(pFile, "%d %d %5.4f\n", j+1, k+1, ((double) strokeTally[j][k]/nAccepted));
                    }
                }
                fclose(pFile);

                // write to histogram matrix file
                pFile = fopen(fileName3.str().c_str(), "w");

                for (j=0;j<maxStrokeWidth;j++)
                {
                    fprintf(pFile, "%5.4f", ((double) strokeTally[j][0]/nAccepted));
                    for (k=1;k<maxStrokeWidth;k++)
                    {
                        fprintf(pFile, " %5.4f", ((double) strokeTally[j][k]/nAccepted));
                    }
                    fprintf(pFile, "\n");
                }
                fclose(pFile);
            }

            // write current evolution
            encodeImage(samples, directoryName, currentImage, width, height);
        }
    }
    // END MAIN LOOP

    cout << "Complete!" << endl;

    return 0;
}
// END MAIN FUNCTION

// FUNCTION DECLARATIONS

// Print help message to stdout
void printHelpMessage()
{
    puts("Rectify version 0.1.0, by Lester O. Hedges\n\n"
         "Syntax: rectify image [options]\n\n"

         "Available options:\n"
         " -h/--help                 : Print this help information\n"
         " -i/--iterations <int>     : Number of MC iterations\n"
         " -f/--frames <int>         : Number of frames to encode\n"
         " -s/--stroke <double>      : Maximum stroke size (fraction of image width)\n"
         " -c/--canvas <string>      : File name of starting canvas\n"
         " -t/--temperature <double> : Temperature of thermal bath\n"
         " -l/--log                  : Sample at logarithmic intervals\n"
         " -d/--directory            : Name of output directory (will be created)\n"
         " -m/--monitor              : Monitor acceptance statistics\n");
}

// Parse arguments from command-line
void parseCommandLineArguments(int argc, char **argv, ostringstream &targetFileName,
    ostringstream &startingFileName, ostringstream &directoryName, long &iterations, unsigned int &frames,
        double &maxStrokeFraction, double &temperature, bool &isLogarithmic, bool &isCanvas, bool &isMonitor)
{
    int i = 1;
    string s;

    // check whether only arg is requesting help
    if (argc == 2)
    {
        if(strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0)
        {
            printHelpMessage();
            exit(0);
        }
    }

    // must have at least 3 arguments
    if (argc < 2)
    {
        printHelpMessage();
        exit(1);
    }
    else
    {
        s = argv[i];
        targetFileName.str(s);
        i++;

        while (i < argc)
        {
            if(strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--iterations") == 0)
            {
                i++;
                iterations = atoi(argv[i]);
            }

            else if(strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--frames") == 0)
            {
                i++;
                frames = atoi(argv[i]);
            }

            else if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--stroke") == 0)
            {
                i++;
                maxStrokeFraction = atof(argv[i]);
            }

            else if(strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"--canvas") == 0)
            {
                i++;
                s = argv[i];
                startingFileName.str(s);
                isCanvas = true;
            }

            else if(strcmp(argv[i],"-d") == 0 || strcmp(argv[i],"--directory") == 0)
            {
                i++;
                s = argv[i];
                directoryName.str(s);

                struct stat st;
                string systemCommand;

                // check to see if folder already exists
                if (stat(directoryName.str().c_str(), &st) != 0)
                {
                    // create top level directory for output
                    systemCommand = "mkdir " + directoryName.str();
                    system(systemCommand.c_str());

                    cout << "Created output directory: " << directoryName.str() << "/\n" << endl;
                }
            }

            else if(strcmp(argv[i],"-t") == 0 || strcmp(argv[i],"--temperature") == 0)
            {
                i++;
                temperature = atof(argv[i]);
            }

            else if(strcmp(argv[i],"-l") == 0 || strcmp(argv[i],"--logarithmic") == 0)
            {
                isLogarithmic = true;
            }

            else if(strcmp(argv[i],"-m") == 0 || strcmp(argv[i],"--monitor") == 0)
            {
                isMonitor = true;
            }

            else
            {
                cerr << "Error: unknown command-line parameter " << argv[i] << endl;
                cerr << "For help run \"rectify -h\"" << endl;
                exit(1);
            }

            i++;
        }
    }
}

// Generate array of linear or logarithmic sampling points
void generateSamplePoints(vector <long> &samplePoints, long &iterations, unsigned int &frames, bool &isLogarithmic)
{
    samplePoints.resize(frames);

    unsigned int i;
    if (isLogarithmic)
    {
        for (i=0;i<frames;i++)
        {
            samplePoints[i] = pow(iterations, ((double) i / (frames-1)));
        }
    }
    else
    {
        unsigned int sampleInterval = iterations/frames;

        for (i=0;i<frames;i++) samplePoints[i] = (i+1)*sampleInterval;
    }
}

// Decode png image from disk
void decodeImage(ostringstream &fileName, vector <unsigned char> &image, unsigned int &width, unsigned int &height)
{
    // decode
    unsigned error = lodepng::decode(image, width, height, fileName.str().c_str());

    //if there's an error, display it
    if(error) cout << "Decoder error " << error << ": " << lodepng_error_text(error) << endl;

    // N.B. lodepng reads in 4 bytes per pixel, i.e. data is stored as RGBARGBA...
}

// Encode png image to disk
void encodeImage(unsigned int &frame, ostringstream &directoryName, vector <unsigned char> &image, unsigned int &width, unsigned int &height)
{
    ostringstream fileName,num;
    num.str("");
    num.width(5);
    num.fill('0');
    num << right << frame;

    fileName.str("");

    if (strcmp(directoryName.str().c_str(), "") == 0)
    {
        fileName << "Frame_" << num.str() << ".png";
    }
    else fileName << directoryName.str() << "/Frame_" << num.str() << ".png";

    // encode the image
    unsigned error = lodepng::encode(fileName.str().c_str(), image, width, height);

    // if there's an error, display it
    if(error) cout << "Encoder error " << error << ": "<< lodepng_error_text(error) << endl;
}

// Apply a random brush stroke to trial image
void applyTrialBrushStroke(vector <unsigned char> &trialImage, unsigned int &width,
    unsigned int &height, MTRand &rng, unsigned int &maxStrokeWidth, int &x1, int &y1, int &x2, int &y2)
{
    // choose random stroke length and width
    int strokeLength = 1 + rng.randInt(maxStrokeWidth-1);
    int strokeWidth = 1 + rng.randInt(maxStrokeWidth-1);

    unsigned int color;

    // choose black or white stroke
    if (rng() < 0.5) color = 255;
    else color = 0;

    // choose random starting point for stroke
    x1 = rng.randInt(width);
    y1 = rng.randInt(height);

    // shift to avoid bias to right of image
    x1 -= 0.5*strokeLength;
    y1 -= 0.5*strokeWidth;

    x2 = x1 + strokeLength;
    y2 = y1 + strokeWidth;

    // enforce boundary restrictions
    if (x1 < 0) x1 = 0;
    else if (x1 > int(width)-1) x1 = width-1;

    if (y1 < 0) y1 = 0;
    else if (y1 > int(height)-1) y1 = height-1;

    if (x2 < 0) x2 = 0;
    else if (x2 > int(width)-1) x2 = width-1;

    if (y2 < 0) y2 = 0;
    else if (y2 > int(height)-1) y2 = height-1;

    // blend block
    blendBlock(trialImage, width, color, rng(), x1, y1, x2, y2);
}

// Blend stroke with existing image block
void blendBlock(vector <unsigned char> &image, unsigned int &width,
    unsigned int &color, double opacity, int &x1, int &y1, int &x2, int &y2)
{
    int x,y;
    unsigned char blendedColour;

    for (y=y1;y<=y2;y++)
    {
        for (x=x1;x<=x2;x++)
        {
            // work out new blended color
            blendedColour = char(opacity*color + (1.0-opacity)*image[4 * width * y + 4 * x + 0]);

            image[4 * width * y + 4 * x + 0] = blendedColour;
            image[4 * width * y + 4 * x + 1] = blendedColour;
            image[4 * width * y + 4 * x + 2] = blendedColour;
            image[4 * width * y + 4 * x + 3] = 255;
        }
    }
}

// Copy stroke region between images
void copyBrushStroke(vector <unsigned char> &from,
    vector <unsigned char> &to, unsigned int &width, int &x1, int &y1, int &x2, int &y2)
{
    int x,y;

    // copy pixel block
    for (y=y1;y<=y2;y++)
    {
        for (x=x1;x<=x2;x++)
        {
            to[4 * width * y + 4 * x + 0] = from[4 * width * y + 4 * x + 0];
            to[4 * width * y + 4 * x + 1] = from[4 * width * y + 4 * x + 1];
            to[4 * width * y + 4 * x + 2] = from[4 * width * y + 4 * x + 2];
            to[4 * width * y + 4 * x + 3] = from[4 * width * y + 4 * x + 3];
        }
    }
}

// Compute absolute difference between two pixel blocks
double computeBlockValue(vector <unsigned char> &image,
    vector <unsigned char> &target, unsigned int &width, int &x1, int &y1, int &x2, int &y2)
{
    int x,y;
    double blockValue = 0;

    // sum pixel block
    for (y=y1;y<=y2;y++)
    {
        for (x=x1;x<=x2;x++)
        {
            // RGB values are the same, only need to consider one vector component
            blockValue += fabs(int(image[4 * width * y + 4 * x + 0]) - int(target[4* width * y + 4 * x + 0]));
        }
    }

    return blockValue/255;
}

// Initialise blank image
void initializeImage(vector <unsigned char> &image, unsigned int &width, unsigned int &height)
{
    unsigned int x,y;

    // create solid black canvas
    for (y=0;y<height;y++)
    {
        for (x=0;x<width;x++)
        {
            image[4 * width * y + 4 * x + 0] = 0;
            image[4 * width * y + 4 * x + 1] = 0;
            image[4 * width * y + 4 * x + 2] = 0;
            image[4 * width * y + 4 * x + 3] = 255;
        }
    }
}

// Convert color image to grayscale
void convertToGrayscale(vector <unsigned char> &image, unsigned int &width, unsigned int &height)
{
    unsigned int x,y;
    unsigned int color;

    for (y=0;y<height;y++)
    {
        for (x=0;x<width;x++)
        {
            color = int((double)image[4 * width * y + 4 * x + 0]*0.2989
                    + (double)image[4 * width * y + 4 * x + 1]*0.5870
                    + (double)image[4 * width * y + 4 * x + 2]*0.1140);

            image[4 * width * y + 4 * x + 0] = color;
            image[4 * width * y + 4 * x + 1] = color;
            image[4 * width * y + 4 * x + 2] = color;
            image[4 * width * y + 4 * x + 3] = 255;
        }
    }
}

// Clear tally counter entries
void resetTallyCounter(vector <vector <long> > &tallyCounter, unsigned int &size)
{
    unsigned int i,j;

    for (i=0;i<size;i++)
    {
        for (j=0;j<size;j++)
        {
            tallyCounter[i][j] = 0;
        }
    }
}
