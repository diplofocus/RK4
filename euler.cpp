#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <time.h>

using namespace std;

const int elements = 2;
const double G = 6.67E-11;
const double dt = 2 * 3600;
const double tmax = 100 * 3.15569e7;
int snapshot = 30;

 class pState
{
        public:
        double x,y;
	double x1,x2,x3,x4;
	double y1,y2,y3,y4;
        double vx,vy;
	double vx1,vx2,vx3,vx4;
	double vy1,vy2,vy3,vy4;
        double ax,ay,a;
	double ax1,ax2,ax3,ax4;
	double ay1,ay2,ay3,ay4;
	double a1,a2,a3,a4;
        double m;

	pState()
	{
        x = rand()/ ((double)RAND_MAX) * (150e3) - 150E3;
        y = rand()/ ((double)RAND_MAX) * (150e3) - 150E3;
        vx = rand()/ ((double)RAND_MAX) * (3E2) - 3E3;
        vy = rand()/ ((double)RAND_MAX) * (3E2) - 3E3;
	m = rand()/ ((double)RAND_MAX) * 7E10;
	}

};

double sqr(double a)
{
    return a * a;
}

double cube(double a)
{
	return a * a * a;
}

double dist(pState A, pState B)
{
    return sqrt(sqr(B.x - A.x) + sqr(B.y - A.y));
}

double dist(double x1, double y1, double x2, double y2)
{
	return sqrt(sqr(x2 - x1) + sqr(y2 - y1));
}

double FGrav(pState A, pState B)
{
	double F = - G * A.m * B.m / cube(dist(A, B));
	return F;
}

int main()
{
	ofstream rout;
	ofstream vout;
	ofstream aout;
	ofstream simout;
	simout.open("out.txt");
	rout.open("rout.txt");
	vout.open("vout.txt");
	aout.open("aout.txt");

    pState particle[elements];

    for(int i = 0; i < elements; i++)
    {
	particle[i] = pState();
    }

	particle[0].x = 0;
	particle[0].y = 0;
	particle[0].vx = 0;
	particle[0].vy = 0;
	particle[0].m = 1.9E30;

	particle[1].x = 150E9;
	particle[1].y = 0;
	particle[1].vx = 0;
	particle[1].vy = 30E3;
	particle[1].m = 6E24;

/*	particle[2].x = 0;
	particle[2].y = 778E9;
	particle[2].vx = 13E3;
	particle[2].vy = 0;
	particle[2].m = 1.9E27;
*/
	int frame = snapshot;
	double perc = 0;
	int currPerc = 0;
	for(double t = 0; t <= tmax; t += dt)
	{
		if(perc >= tmax / 100)
		{
			cout.flush();
			currPerc++;
			cout << currPerc << endl;
			perc=0;
		}


		if(frame % snapshot == 0)
		{
			for(int i = 0; i < elements; i++)
			simout << "#Telo#" << i << endl << particle[i].x << ", " << particle[i].y << endl;
		}

		perc += dt;
	for(int i = 0; i < elements; i++)
	{
		particle[i].ax1 = 0;
		particle[i].ay1 = 0;
		particle[i].ax2 = 0;
		particle[i].ay2 = 0;
		particle[i].ax3 = 0;
		particle[i].ay3 = 0;
		particle[i].ax4 = 0;
		particle[i].ay4 = 0;
	}

	//Determining acceleration
	    for(int i = 0; i < elements; i++)
	    {
	        for(int j = 0; j < elements; j++)
	        {
			if(j==i)
			{
			continue;
			}

			particle[i].ax1 += (-G * particle[j].m * (particle[i].x - particle[j].x))
							 / cube(dist(particle[i], particle[j]));

			particle[i].ay1 += (-G * particle[j].m * (particle[i].y - particle[j].y))
							 / cube(dist(particle[i], particle[j]));

	        }
			particle[i].vx1 = particle[i].vx;
			particle[i].vy1 = particle[i].vy;
			particle[i].x1 = particle[i].x;
			particle[i].y1 = particle[i].y;

	    }

//k2

	    for(int i = 0; i < elements; i++)
	    {
		particle[i].vx2 = particle[i].vx1 + particle[i].ax1 * dt/2.0;
		particle[i].vy2 = particle[i].vy1 + particle[i].ay1 * dt/2.0;

		particle[i].x2 = particle[i].x1 + particle[i].vx2 * dt/2.0;
		particle[i].y2 = particle[i].y1 + particle[i].vy2 * dt/2.0;
	    }

	     for(int i = 0; i < elements; i++)
	     {
	        for(int j = 0; j < elements; j++)
	        {
			if(j==i)
			{
			continue;
			}
//			cout << particle[i].x2 << ", " << particle[i].y2 << ", " << particle[j].x2 << ", " << particle[j].y2 << endl;


			particle[i].ax2 += (-G * particle[j].m * (particle[i].x2 - particle[j].x2))
						 / cube(dist(particle[i].x2, particle[i].y2, particle[j].x2, particle[j].y2));
			particle[i].ay2 += (-G * particle[j].m * (particle[i].y2 - particle[j].y2)) 
						 / cube(dist(particle[i].x2, particle[i].y2, particle[j].x2, particle[j].y2));

		}
	    }
//k3
	    for(int i = 0; i < elements; i++)
	    {
		particle[i].vx3 = particle[i].vx1 + particle[i].ax2 * dt/2.0;
		particle[i].vy3 = particle[i].vy1 + particle[i].ay2 * dt/2.0;

		particle[i].x3 = particle[i].x1 + particle[i].vx3 * dt/2.0;
		particle[i].y3 = particle[i].y1 + particle[i].vy3 * dt/2.0;
	    }
		for(int i = 0; i < elements; i++)
	    {

	        for(int j = 0; j < elements; j++)
	        {
			if(j==i)
			{
			continue;
			}

			particle[i].ax3 += (-G * particle[j].m * (particle[i].x3 - particle[j].x3)) / cube(sqrt(sqr(particle[i].x3 - particle[j].x3) + sqr(particle[i].y3 - particle[j].y3)));
			particle[i].ay3 += (-G * particle[j].m * (particle[i].y3 - particle[j].y3)) / cube(sqrt(sqr(particle[i].x3 - particle[j].x3) + sqr(particle[i].y3 - particle[j].y3)));
		}
	   }


//k4

	    for(int i = 0; i < elements; i++)
	    {
		particle[i].vx4 = particle[i].vx1 + particle[i].ax3 * dt;
		particle[i].vy4 = particle[i].vy1 + particle[i].ay3 * dt;

		particle[i].x4 = particle[i].x1 + particle[i].vx4 * dt;
		particle[i].y4 = particle[i].y1 + particle[i].vy4 * dt;
	    }
	    for(int i = 0; i < elements; i++)
	   {
	        for(int j = 0; j < elements; j++)
	        {
			if(j==i)
			{
			continue;
			}

			particle[i].ax4 += (-G * particle[j].m * (particle[i].x4 - particle[j].x4)) / cube(sqrt(sqr(particle[i].x4 - particle[j].x4) + sqr(particle[i].y4 - particle[j].y4)));
			particle[i].ay4 += (-G * particle[j].m * (particle[i].y4 - particle[j].y4)) / cube(sqrt(sqr(particle[i].x4 - particle[j].x4) + sqr(particle[i].y4 - particle[j].y4)));
		}
	   }


	for(int i = 0; i < elements; i++)
		{
		particle[i].vx += dt/6.0 * (particle[i].ax1 + 2*particle[i].ax2 + 2*particle[i].ax3 + particle[i].ax4);
		particle[i].vy += dt/6.0 * (particle[i].ay1 + 2*particle[i].ay2 + 2*particle[i].ay3 + particle[i].ay4);

		particle[i].x += dt/6.0 * (particle[i].vx1 + 2*particle[i].vx2 + 2*particle[i].vx3 + particle[i].vx4);
		particle[i].y += dt/6.0 * (particle[i].vy1 + 2*particle[i].vy2 + 2*particle[i].vy3 + particle[i].vy4);

//                particle[i].x += particle[i].vx * dt;
//                particle[i].y += particle[i].vy * dt;
//		simout << "#Telo#" << i << endl << particle[i].x << ", " << particle[i].y << endl;

//		aout << "#Telo#" << i << endl << particle[i].ax1 << ", " << particle[i].ax2 << ", " << particle[i].ax3 << ", " << particle[i].ax4 << endl << particle[i].ay1 << ", " << particle[i].ay2 << ", " << particle[i].ay3 << ", " << particle[i].ay4 << endl;
//		rout << "#Telo#" << i << endl << particle[i].x1 << ", " << particle[i].x2 << ", " << particle[i].x3 << ", " << particle[i].x4 << endl << particle[i].y1 << ", " << particle[i].y2 << ", " << particle[i].y3 << ", " << particle[i].y4 << endl;
//		vout << "#Telo#" << i << endl << particle[i].vx1 << ", " << particle[i].vx2 << ", " << particle[i].vx3 << ", " << particle[i].vx4 << endl << particle[i].vy1 << ", " << particle[i].vy2 << ", " << particle[i].vy3 << ", " << particle[i].vy4 << endl;
		}
	frame++;
  }
	rout.close();
	vout.close();
	aout.close();
	simout.close();
return 0;
}
