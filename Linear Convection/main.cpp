//Coding for solving the convection linear equation from Navier-Stokes
#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
#include <numeric>

using namespace std;

int main()
{
    std::string line;               // A line of key/values from text
    std::string key;                // Temporary for our key
    std::string value;              // Temporary for our value
    std::ifstream stream(Parameters);     // Load the file stream
    std::stringstream splitter;     // Prepare a stringstream as a splitter (splits on spaces) for reading key/values from a line

    // Make sure we can read the stream
    if (stream) {
        // As long as there are lines of data, we read the file
        while (std::getline(stream, line)) {
            splitter << line;                                   // Load line into splitter
            splitter >> key;                                    // Read the key back into temporary
            splitter >> value;                                  // Read the value back into temporary
            splitter.clear();                                   // Clear for next line
            variables[key] = value;                             // Store the key/value pair in our variable map.
        }
    }
    else {
        // The file was not found or locked, etc...
        std::cout << "Unable to open file: " << path << std::endl;
    }
    // Defining some initial parameters
    float endtime=0.3; float dx=0.025; float xmax=2; float c=1; float C=0.3;
    float dt=C*dx/c; int nx=xmax/dx + 1; int nt=endtime/dt;
    float u[nx] = { }; //defining an initial velocity array.
    float un[nx] = { };

    for(int i = 0; i <= 9; i++)
    {
        u[i] = 2;
    }
    for(int n=1; n<=nt; n++)
    {
        std::copy(u, u + nx, un);
        for(int i=1; i<=nx; i++)
        {
            u[i] = un[i] - c*(dt/dx)*(un[i]-un[i-1]);

        }
    }
}
