/*
MatzJB 21/4 2019
 2D Incompressible fluid simulation 
 
 */

final int maxDensity = 40;
final int maxVelocity = 100;
final int N = 200; //resolution of board
final int iteration = 16; //number of iterations per time step
Fluid fluid;
final int scale = 5; //size of pixels

void setup() 
{
  fluid = new Fluid(0.003, 0, 0.0);   
  frameRate(30);
  settings();
}

void settings()
{
  size(scale*N, scale*N);
}


void mouseDragged()
{
  int size =3;
  int x = mouseX/scale;
  int y = mouseY/scale;
  for(int i=x-size; i<x+size; i++)
    for(int j=y-size; j<y+size; j++)
    {
      fluid.AddDensity(x, y, maxDensity);
      float diffX = mouseX-pmouseX;
      float diffY = mouseY-pmouseY;
      fluid.AddVelocity(x, y, diffX, diffY);
    }
   
  
}


void draw() 
{
  background(0);
//  for(int i=0; i<5; i++)
  {
    fluid.AddDensity(frameCount%10, N/2, 1000);
    
   
      fluid.AddVelocity(frameCount%10, N/2, 5, 0);
    
}
  //fluid.AddVelocity(0, N/2, 500, 0);
  

  fluid.AdvanceSimulation();
  //fluid.Fade(0.99);
  //fluid.AddWallDensity();
  //fluid.AddObstruction();
  fluid.Render();
}
