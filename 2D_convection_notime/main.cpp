//Code for solving the 2D diffusion Equation with no time term

#include <iostream>
#include <cmath>

using namespace std;

float nu = 1; //viscosity
const int nx = 31; //number of points in x
const int ny = 31; //number of points in y
float xmin = 0.0; //domain size in x
float ymin = 0.0; //domain size in y 
float xmax = 1.0; //domain size in x
float ymax = 1.0; //domain size in y 
float dx = (xmax-xmin)/(nx - 1); //length of the element in x
float dy = (ymax-ymin)/(ny - 1); //length of the element in y

float l2_norm(float matA[nx][ny], float matB[nx][ny]);
float L_inf(float matA[nx][ny], float matB[nx][ny]);
float mat_display(float matA[nx][ny]);
float Ue(float xc, float yc);
float sterm(float xc, float yc);
float neumann(float xc, float yc);

int main() 
{
    //WITH NO SOURCE TERM

    //Initializing domain
    float x[nx], y[ny];
    float u[nx][ny] = {}; //current velocity
    float un[nx][ny];  //laststep velocity
    float u_diff[nx][ny]; //difference between current velocity and laststep velocity
    
    int cc[4];
    cout << "Type the contour conditions for each wall (1 - Dirichlet, 2 - Neumann)" << endl;
    cout << "Left Wall: ";
    cin >> cc[0];
    cout << "Right Wall: ";
    cin >> cc[1];
    cout << "Bottom Wall: ";
    cin >> cc[2];
    cout << "Top Wall: ";
    cin >> cc[3];

    for (int i = 0; i<nx; i++) 
    {
        x[i] = xmin + i*dx;
        for (int j=0; j<ny; j++)
        {
            y[j] = ymin + j*dy;
            if (i==0) //left wall 
            {
                if (cc[0]==1)
                {
                    u[i][j] = Ue(xmin,y[j]);
                }
                else if (cc[0]==2)
                {
                    u[i][j] = neumann(xmin,y[j]);
                }
            }
            else if (i==nx-1) //right wall 
            {
                if (cc[1]==1)
                {
                    u[i][j] = Ue(xmax,y[j]);
                }
                else if (cc[1]==2)
                {
                    u[i][j] = neumann(xmax,y[j]);
                }
            }             
            else if (j==0) //bottom wall
            {
                if (cc[2]==1)
                {
                    u[i][j] = Ue(x[i],ymin);
                }
                else if (cc[2]==2)
                {
                    u[i][j] = neumann(x[i],ymin);
                }                

            }
            else if (j==ny-1) //top wall
            {
                if (cc[3]==1)
                {
                    u[i][j] = Ue(x[i],ymax);
                }
                else if (cc[3]==2)
                {
                    u[i][j] = neumann(x[i],ymax);
                }
            }
            else 
            {
                u[i][j] = 1; //inside domain
            }
        }
    }
    
    // cout << "Initial conditions and contour conditions: " << endl;
    // mat_display(u);
 
    cout <<endl;

    float e = 1000; //initializing error

    while (e>=0.00000000001)
    {
        std::copy(&u[0][0], &u[0][0] + nx*ny, &un[0][0]); //copying the current velocity to the laststep velocity matrix
 
        //Calculating the current velocity
        for (int j = 1; j<ny-1; j++)
        {
            for (int i = 1; i<nx-1; i++)
            { 
                u[i][j] = (sterm(x[i],y[j]) + un[i+1][j]/pow(dx,2) + un[i-1][j]/pow(dx,2) + un[i][j+1]/pow(dy,2) + un[i][j-1]/pow(dy,2))/((2/pow(dx,2)) + (2/pow(dy,2)));
            }
        }
         
        e = l2_norm(u, un);

    }

    //WITH SOURCE TERM
    float analytical[nx][ny] = {};

    //Calculating the analytical velocity
    for (int i = 0; i<nx; i++) 
    {
        x[i] = xmin + i*dx;
        for (int j=0; j<ny; j++)
        {
            y[j] = ymin + j*dy;
            if (i==0) //left wall 
            {
                if (cc[0]==1)
                {
                    analytical[i][j] = Ue(xmin,y[j]);
                }
                else if (cc[0]==2)
                {
                    analytical[i][j] = neumann(xmin,y[j]);
                }
            }
            else if (i==nx-1) //right wall 
            {
                if (cc[1]==1)
                {
                    analytical[i][j] = Ue(xmax,y[j]);
                }
                else if (cc[1]==2)
                {
                    analytical[i][j] = neumann(xmax,y[j]);
                }
            }             
            else if (j==0) //bottom wall
            {
                if (cc[2]==1)
                {
                    analytical[i][j] = Ue(x[i],ymin);
                }
                else if (cc[2]==2)
                {
                    analytical[i][j] = neumann(x[i],ymin);
                }                

            }
            else if (j==ny-1) //top wall
            {
                if (cc[3]==1)
                {
                    analytical[i][j] = Ue(x[i],ymax);
                }
                else if (cc[3]==2)
                {
                    analytical[i][j] = neumann(x[i],ymax);
                }
            }
            else 
            {
                analytical[i][j] = Ue(x[i],y[j]); //inside domain
            }
        }
    } 
    
    cout << "Mean error by L-2 norm: " << l2_norm(u, analytical) << endl;
    cout << "Error by L_inf: " << L_inf(u, analytical) << endl;

    //mat_display(u);
    //mat_display(analytical);

}  

//----------Function for computing L2-Norm-----------//
float l2_norm(float matA[nx][ny], float matB[nx][ny])
{
    float diff_mat[nx][ny];
    for (int j = 0; j<ny; j++)
    {
        for (int i = 0; i<nx; i++)
        { 
             diff_mat[i][j]=pow(matA[i][j] - matB[i][j],2);
        }
    }

    float sum_diff = 0; //defining a variable for the sum of all elements of the u_diff matrix
    for(int j=0; j<ny; j++)
    {
     	for(int i=0; i<nx; i++)
        {
     		sum_diff=(diff_mat[i][j] + sum_diff);
 		}
 	}

    float mean_diff = sum_diff/(nx*ny); //taking the mean 
        
    float error = pow(mean_diff, 0.5); //aplying the square root of the L2-norm

    return error; 
    
}

//-----------------------------------------------//


//----------Function for computing Linf error----------//
float L_inf(float matA[nx][ny], float matB[nx][ny])
{

    float diff_mat[nx][ny];
    for (int j = 0; j<ny; j++)
    {
        for (int i = 0; i<ny; i++)
        { 
             diff_mat[i][j]=(matA[i][j] - matB[i][j]);
        }
    }

    float maxElement = 0; 
    for (int i = 0; i<nx; i++) { 
        for (int j = 0; j<ny; j++) { 
            if (diff_mat[i][j] > maxElement) { 
                maxElement = diff_mat[i][j]; 
            } 
        } 
    }

    return maxElement;

}
//-----------------------------------------------------//


//------Function for displaying matrices-------//
float mat_display(float matA[nx][ny]){

    for(int x=0;x<nx;x++)  // loop 3 times for three lines
    {
        for(int y=0;y<ny;y++)  // loop for the three elements on the line
        {
            cout<<matA[x][y];  // display the current element out of the array
        }
    cout<<endl;  // when the inner loop is done, go to a new line
    }

}
//--------------------------------------------//

//----------Function for computing analytical function Ue ----------//
float Ue(float xc, float yc)
{
    float fxy;
    fxy = sin(xc)*cos(yc);

    return fxy;

}
//-----------------------------------------------------------------//

//----------Function for computing the source term ----------------//
float sterm(float xc, float yc)
{
    float fxy;
    fxy = 2.0*sin(xc)*cos(yc);

    return fxy;

}
//-----------------------------------------------------------------//

//-------Function for computing CC with Neumann Condition----------//
float neumann(float xc, float yc)
{
    float neu_u;
    if(xc==xmin || xc==xmax)
    {
        neu_u = cos(xc)*cos(yc);    
    }
    else if(yc==ymin || yc==ymax)
    {
        neu_u = -sin(xc)*sin(yc);
    }

    return neu_u;
}
//-----------------------------------------------------------------//