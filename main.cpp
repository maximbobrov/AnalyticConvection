
#include <stdio.h>
#include <stdlib.h>

#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
//#include <GL/glut.h>
#include  <math.h>
#include <time.h>

#define NR 100
#define NZ 100

double Om[NR][NZ];
double Uz[NR][NZ];
double Ur[NR][NZ];
double B[NR][NZ];

double rhsOm[NR][NZ];
double rhsUz[NR][NZ];
double rhsB[NR][NZ];
double bcOm[NR];
double bcB[NR];

double dr = 0.1;
double dz = 0.1;

bool showOm = false;
bool showUz = false;
bool showB = true;
bool go = false;

double sc=0.5;
double t_x=15.0;
double t_y=15.0;

double nu,k,N;

void display(void);

void init();


void sweep(int itn, double alpha);

double solve_incompr(int itn);


double ck=1.0;//0.250;
double nu0=0.001;
int j_curr=0;

//bcType ;//0 fixed value, 1 fixed dif (a1-a0)=value;
double sweep_imp_gs(int i,double arr[NR][NZ],int itn,double rhs_[NR][NZ],double alpha, int bcType,double bcVal) //universal implicit sweep(gauss seidel)
{
    int j,n;
    double a,b;

    a=(2.0);
    b=1.0;

    //a*x[i]=b*(x[i-1]+x[i+1])+rhs
    //  a*x[i]-b*(x[i-1])=b*(x[i+1])+rhs
    //  a*x[i]-b*(-value+(x[i]-x[i-1])+x[i-1])=b*(x[i+1])+rhs
    //    a*x[i]-b*x[i]=b*(x[i+1]-value)+rhs

    for (n=0;n<itn;n++)
    {

      /*  j=1;
        if (bcType==1)
         arr[i][j]=arr[i][j]*(1.0-alpha)+ alpha*(b*(arr[i][j+1]-bcVal) +  rhs_[i][j])/(a-b);
        else
         arr[i][j]=arr[i][j]*(1.0-alpha)+ alpha*(b*(arr[i][j+1]+bcVal) +  rhs_[i][j])/a;*/

        for (j=1;j<NZ-1;j++)
        {
            arr[i][j]=arr[i][j]*(1.0-alpha)+ alpha*(b*(arr[i][j+1]+arr[i][j-1]) +  rhs_[i][j])/a;
        }
        arr[i][NZ-1]=arr[i][NZ-2];
    }

    double res=0.0;
    for (j=1;j<NZ-1;j++)
    {
        res+=alpha*(b*(arr[i][j+1]+arr[i][j-1]) +  rhs_[i][j])-a*arr[i][j];
    }

    return res; //residual
}


void get_ur(int i)
{
    Ur[i][0]=0.0;
    for (int j=1;j<NZ-1;j++)
    {
        Ur[i][j]=Ur[i][j-1]+0.5*dz*(Om[i][j-1]+Om[i][j]);
    }

}

void get_ur_2D()
{
    //Ur[i][0]=0.0;
    for (int j=1;j<NZ-1;j++)
    {
        Ur[0][j]=0.0;
        for(int i=1; i<NR-1; i++)
        {
        Ur[i][j]=Ur[i-1][j]+0.5*dr*((Uz[i][j+1]-Uz[i][j-1])/(2.0*dz)+(Uz[i-1][j+1]-Uz[i-1][j-1])/(2.0*dz));
        }
    }

}


void display(void)
{
    int i,j;
    /* clear window */

    glClear(GL_COLOR_BUFFER_BIT);

    glLoadIdentity();

    glScalef(sc,sc,sc);
    glTranslatef(t_x,t_y,0);
    glColor3f(1,1,1);

    glLineWidth(1.0);
    glColor3f(1,0,0);
    glColor3f(0,1,0);
    glPointSize(2.0);


    for(i=0; i<NR-1; i++)
    {
            glBegin(GL_TRIANGLE_STRIP);
        for (j=0; j<NZ; j++)
        {
            if(showB)
                glColor3f(ck * B[i][j],-ck * B[i][j],-ck * B[i][j]);
            if(showOm)
                glColor3f(-ck * Om[i][j],ck * Om[i][j],-ck * Om[i][j]);
            if(showUz)
                glColor3f(ck * Ur[i][j],-ck * Ur[i][j],-ck * Ur[i][j]);

            glVertex3f( i * dr, j * dz, 0.0);

            if(showB)
                glColor3f(ck * B[i+1][j],-ck * B[i+1][j],-ck * B[i+1][j]);
            if(showOm)
                glColor3f(-ck * Om[i+1][j],ck * Om[i+1][j],-ck * Om[i+1][j]);
            if(showUz)
                glColor3f(ck * Ur[i+1][j],-ck * Ur[i+1][j],-ck * Ur[i+1][j]);

            glVertex3f( (i+1) * dr, j * dz, 0.0);
        }
        glEnd();
    }


    glBegin(GL_LINES);
    glLineWidth(2);
    glColor3f(1.0,1.0,1.0);
    glVertex3f( -0.001, 0.0, 0.0);
    glVertex3f( -0.001, NZ * dz, 0.0);
    glEnd();

    if(go)
    {
        for (int i = 0; i< 10 ; i++)
            sweep(5,0.5);
        glutPostRedisplay();
    }

    glutSwapBuffers();
}

double amp=1.0;
void kb(unsigned char key, int x, int y)
{

    if (key=='1')
    {
        showB = true;
        showOm = false;
        showUz = false;
    }
    if (key=='2')
    {
        showB = false;
        showOm = true;
        showUz = false;
    }
    if (key=='3')
    {
        showB = false;
        showOm = false;
        showUz = true;
    }
    if (key=='q')
    {
        sc*=1.01;
    }
    if (key=='e')
    {
        sc/=1.01;
    }
    if (key=='w')
    {
        t_y-=5.0;
    }
    if (key=='s')
    {
        t_y+=5.0;
    }
    if (key=='a')
    {
        t_x+=5.0;
    }
    if (key=='d')
    {
        t_x-=5.0;
    }
    if (key==' ')
    {
        printf(" sweep \n");
        for (int i=0;i<50;i++)
        sweep(1,0.001);

       // go =! go;
    }
    if (key=='.')
    {
        ck*=1.3;
        printf("ck=%e \n",ck);
    }
    if (key==',')
    {
        ck/=1.3;
        printf("ck=%e \n",ck);
    }
    glutPostRedisplay();
}

/*double timeStep()
{
    for (int i=0; i<NR; i++)
    {
        for (int j=0; j<NZ; j++)
        {
            UrPrev[i][j]=Ur[i][j];
            UzPrev[i][j]=Uz[i][j];
            BPrev[i][j]=B[i][j];
        }
    }
}

double solveUr(int itn)
{
    int i,j,n;
    double res=0;
    static double rhs[NR+1][NZ+1];
    double nu = 1e-2;

    for (i=1; i<NR-1; i++)
    {
        Ur[i][0] = 0.0;
        Ur[i][1] = 0.0;
        Ur[i][NZ - 2] = 0.0;
        Ur[i][NZ - 1] = 0.0;
    }

    for (j=2; j<NZ-2; j++)
    {
        Ur[0][j] = Ur[1][j];
        Ur[NR-1][j] = Ur[NR-2][j];
    }

    for (i=1; i<NR-1; i++)
    {
        for (j=2; j<NZ-2; j++)
        {
            rhs[i][j] = UrPrev[i][j] + dt * ( BPrev[i + 1][j] - BPrev[i - 1][j]) / (2 * dr);
        }
    }

    for(n=0;n<itn;n++)
    {
        for (i=1; i<NR-1; i++)
        {
            for (j=2; j<NZ-2; j++)
            {
                Ur[i][j]=(rhs[i][j] + dt * nu * (-0.5 * Ur[i][j - 2] + Ur[i][j - 1] - Ur[i][j + 1] + 0.5 * Ur[i][j + 2]) / (dz * dz * dz));
            }
        }
    }

    for (j=2; j<NZ-2; j++)
    {
        Ur[0][j] = Ur[1][j];
        Ur[NR-1][j] = Ur[NR-2][j];
    }
    return res;
}

double solveUz(int itn)
{
    int i,j,n;
    double res=0;
    static double rhs[NR+1][NZ+1];

    for (i=1; i<NR-1; i++)
    {
        for (j=1; j<NZ-1; j++)
        {
            rhs[i][j]=0.0;
        }
    }

    for(n=0;n<itn;n++)
    {
        for (i=1; i<NR-1; i++)
        {
            for (j=1; j<NZ-1; j++)
            {
                Ur[i][j]=(rhs[i][j] - 0.0);
            }
        }
    }
    return res;
}

double solveB(int itn)
{
    int i,j,n;
    double res=0;
    static double rhs[NR+1][NZ+1];

    for (i=1; i<NR-1; i++)
    {
        //B[i][0] = 0.0;
        B[i][1] = B[i][0];
        B[i][NZ - 2] = 0.0;
        B[i][NZ - 1] = 0.0;
    }

    for (j=2; j<NZ-2; j++)
    {
        B[0][j] = B[1][j];
        B[NR-1][j] = B[NR-2][j];
    }

    double N = 0.001;
    double k = 1e-2;
    for (i=1; i<NR-1; i++)
    {
        for (j=2; j<NZ-2; j++)
        {
            rhs[i][j] = BPrev[i][j] - dt * N * N / (dr * i) * ((dr * (i + 1))  * UrPrev[i + 1][j] - (dr * (i - 1)) * UrPrev[i - 1][j]) / (2 * dr);
        }
    }

    for(n=0;n<itn;n++)
    {
        for (i=1; i<NR-1; i++)
        {
            for (j=2; j<NZ-2; j++)
            {
                B[i][j] = (rhs[i][j] + dt * k * (-0.5 * B[i][j - 2] + B[i][j - 1] - B[i][j + 1] + 0.5 * B[i][j + 2]) / (dz * dz * dz));
            }
        }
    }

    for (j=2; j<NZ-2; j++)
    {
        B[0][j] = B[1][j];
        B[NR-1][j] = B[NR-2][j];
    }
    return res;
}
*/
void sweep(int itn,double alpha)
{
    int nn;



    //for (nn=0;nn<itn;nn++)
    {
        for (int i=0;i<NR;i++)
        {
            if (i==0)
            {
                for (int j=1;j<NZ;j++)
                {
                    rhsOm[i][j]=0.0;
                }

                sweep_imp_gs(i,Om,itn,rhsOm,alpha,0,bcOm[i]);

               //  get_ur(i);
                for (int j=1;j<NZ;j++)
                {
                    rhsUz[i][j]=-Om[i][j]/dr*dz*dz -20.0*dz*dz*Ur[i][j]/((i+1)*dr);

                }
                sweep_imp_gs(i,Uz,itn,rhsUz,alpha,0,0.0);

                for (int j=1;j<NZ;j++)
                {
                    rhsB[i][j]=-Uz[i][j]*dz*dz*N*N/k;// -Ur[i][j]/((i+1)*dr);
                }
                sweep_imp_gs(i,B,itn,rhsB,alpha,1,bcB[i]);
            }else
            {
                for (int j=1;j<NZ;j++)
                {
                    rhsOm[i][j]=-(B[i][j]-B[i-1][j])/dr*dz*dz/nu;
                }
                sweep_imp_gs(i,Om,itn,rhsOm,alpha,0,bcOm[i]);
               // get_ur(i);
                for (int j=1;j<NZ;j++)
                {
                    rhsUz[i][j]=-(Om[i][j]-Om[i-1][j])/dr*dz*dz - 20.0*dz*dz*Ur[i][j]/((i+1)*dr);
                }
                sweep_imp_gs(i,Uz,itn,rhsUz,alpha,0,0.0);

                for (int j=1;j<NZ;j++)
                {
                    rhsB[i][j]=-Uz[i][j]*dz*dz*N*N/k;
                }
                sweep_imp_gs(i,B,itn,rhsB,alpha,0,bcB[i]);
            }
        }
        get_ur_2D();
    }
    //timeStep();
}

void init()
{

     N = 0.01;
     k = 1e-2;
     nu = 1e-2;

    for (int i=0; i<NR; i++)
    {
        for (int j=0; j<NZ; j++)
        {
            B[i][j]=0.0;
            Om[i][j]=0.0;
            Uz[i][j]=0.0;
            if(j == 0)
            {
                B[i][j] = 100.0 * 1.0 / (fabs(i - NR / 2) + 1);
                bcB[i]=100.0 * 1.0 / (fabs(i - NR / 2) + 1);
                bcOm[i]=0.02*pow(nu,-1.0/3.0)*pow(bcB[i],2.0/3.0);
                Om[i][j] = bcOm[j];

            }

            //B[i][j] += 10.0 * rand() / RAND_MAX;
        }
    }
    /* set clear color to black */
    glClearColor (0.0, 0.1, 0.0, 0.0);

    /* set fill color to white */
    glColor3f(1.0, 1.0, 1.0);

    /* set up standard orthogonal view with clipping */
    /* box as cube of side 2 centered at origin */
    /* This is default view and these statement could be removed */
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();

    glOrtho(-1.1,  dr*NR, -1.1, dz*NZ, -10.0, 10.0);

    glMatrixMode (GL_MODELVIEW);

}
int main(int argc, char** argv)
{
    /* Initialize mode and open a window in upper left corner of screen */
    /* Window title is name of program (arg[0]) */


    srand(time(NULL));
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    //  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(800,800);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
   // glutMotionFunc(m_m);
   // glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();
}