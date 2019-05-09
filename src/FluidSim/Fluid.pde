

int flattenIndex(int x, int y) 
{
  x = constrain(x, 0, N-1);
  y = constrain(y, 0, N-1);
  return x + y*N;
}

class Fluid
{
    //int size;
    float dt;
    float diffusion;
    float viscocity;
    
    float[] s;
    float[] density;
    
    float[] Vx;
    float[] Vy;
    
    float[] Vx0;
    float[] Vy0;

  
  Fluid(float dt, float diffusion, float viscocity)
  {
      this.dt = dt;
      this.diffusion = diffusion;
      this.viscocity = viscocity;
      
      this.s = new float[N*N];
      this.density = new float[N*N];
      
      this.Vx = new float[N*N];
      this.Vy = new float[N*N];
           
      this.Vx0 = new float[N*N];
      this.Vy0 = new float[N*N];
  }
  
  // saturate the fluid
  void AddDensity(int x, int y, float amount)
  {
    this.density[flattenIndex(x, y)] += amount;
    this.density[flattenIndex(x, y)] = constrain(this.density[flattenIndex(x, y)], 0, maxDensity);
  }
  
  //todo: replace with vectors
  void AddVelocity(int x, int y, float amountX, float amountY)
  {
    int index = flattenIndex(x, y);
    
    this.Vx[index] += amountX;
    this.Vy[index] += amountY;
  }
  

  /* called "step" in the paper 
  Take a step in the simulation */
  void AdvanceSimulation()
  {

    float viscocity = this.viscocity;
    float diffusion = this.diffusion;
    float dt= this.dt;
    float[] Vx = this.Vx;
    float[] Vy =  this.Vy;
    
    float[] Vx0 = this.Vx0;
    float[] Vy0 = this.Vy0;
    
    float[] s = this.s;
    float[] density = this.density;
        
      diffuse(1, Vx0, Vx, viscocity, dt);
      diffuse(2, Vy0, Vy, viscocity, dt);
      
      project(Vx0, Vy0, Vx, Vy);
      
      advect(1, Vx, Vx0, Vx0, Vy0, dt);
      advect(2, Vy, Vy0, Vx0, Vy0, dt);
      
      project(Vx, Vy, Vx0, Vy0); //vy0 instead of vy!!
      
      diffuse(0, s, density, diffusion, dt);
      advect(0, density, s, Vx, Vy, dt);
  }
   

// Visualizing density matrix
void Render()
{
  
  float m = millis()*0.5;
  
  for (int i = 0; i < N-1; i++) 
  {
    for (int j = 0; j < N-1; j++) 
    {
      float x = i*scale;
      float y = j*scale;
      float density = 255*this.density[flattenIndex(i, j)]/maxDensity;
    
      float vx = 255*constrain(this.Vx[flattenIndex(i, j)], 0, 1);
      float vy = 255*constrain(this.Vy[flattenIndex(i, j)], 0, 1);
   
      fill(vx, vy, (vx*vy), (vx*vy));
      noStroke();
      square(x, y, scale);
    }
  }
}


void AddObstruction()
{
  for(int i=0; i < N-1; i++)
  {
    for(int j=0; j < N-1; j++)
    {
      this.density[flattenIndex(i, j)] = (sqrt(N/2+i) + sqrt(N/2+j))/30;
    }
  }
}


void AddWallDensity()
{
  float d = 1000;
    for (int j = 0; j < N-1; j++) 
    {
      this.Vx[flattenIndex(0, j)] = -abs(this.Vx[flattenIndex(0, j)]);
      //this.Vx[flattenIndex(N-1, j)] = -abs(this.Vx[flattenIndex(N-1, j)]);
       
 //     this.Vy[flattenIndex(j, 1)] = -abs(this.Vy[flattenIndex(j, 1)]);
  //    this.Vy[flattenIndex(j, N-1)] = -abs(this.Vy[flattenIndex(j, N-1)]);
     // this.density[flattenIndex(1, j)] = 100;
      /*
      this.density[flattenIndex(N-10, j)] = 0;
      this.density[flattenIndex(j, 0)] = 0;
      this.density[flattenIndex(j, N-10)] = 0;
      */
    }
}

void Fade(float dampening)
{
  
  for(int i=0; i<N-1; i++)
  {
    for (int j = 0; j < N-1; j++) 
    {
      this.density[flattenIndex(i, j)] *= dampening;
    }
  }
}


}//class


// following functions are generic and are outside of the class definition

// figure out velocity by stepping backwards
void advect(int b, float[] d, float[] d0,  float[] velocX, float[] velocY, float dt)
{
    float Nfloat = N;  
    float i0, i1, j0, j1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float ifloat, jfloat;
    int i, j;
    
        for(j = 1, jfloat = 1; j <= N; j++, jfloat++) 
        { 
            for(i = 1, ifloat = 1; i <= N; i++, ifloat++) 
            {
                tmp1 = dtx * velocX[flattenIndex(i, j)];
                tmp2 = dty * velocY[flattenIndex(i, j)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5f) 
                  x = 0.5f;
                  
                if(x > Nfloat + 0.5f) 
                  x = Nfloat + 0.5f;
                  
                i0 = floor(x); 
                i1 = i0 + 1.0f;
                
                if (y < 0.5f) 
                  y = 0.5f; 
                  
                if(y > Nfloat + 0.5f) 
                  y = Nfloat + 0.5f; 
                
                j0 = floor(y);
                j1 = j0 + 1.0f; 
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                
                int i0i = int(i0);
                int i1i = int(i1);
                int j0i = int(j0);
                int j1i = int(j1);

                // u1=1-u0
                d[flattenIndex(i, j)] = 
                  s0*(t0*d0[flattenIndex(i0i, j0i)] + t1*d0[flattenIndex(i0i, j1i)]) + 
                  s1*(t0*d0[flattenIndex(i1i, j0i)] + t1*d0[flattenIndex(i1i, j1i)]);
            }
        }
    
    setBounds(b, d);
}

    
void setBounds(int b, float[] x)
  {
      for(int j = 1; j < N - 1; j++) {
          for(int i = 1; i < N - 1; i++) {
              x[flattenIndex(i, j)] = b == 3 ? -x[flattenIndex(i, j)] : x[flattenIndex(i, j )];
              x[flattenIndex(i, j)] = b == 3 ? -x[flattenIndex(i, j)] : x[flattenIndex(i, j)];
          }
      }
    
      for(int k = 1; k < N - 1; k++) {
          for(int j = 1; j < N - 1; j++) {
              x[flattenIndex(0, j)] = b == 1 ? -x[flattenIndex(1, j)] : x[flattenIndex(1, j)];
              x[flattenIndex(N-1, j)] = b == 1 ? -x[flattenIndex(N-2, j)] : x[flattenIndex(N-2, j)];
          }
      }

      x[flattenIndex(0, N-1)]     = 0.33f * (x[flattenIndex(1, N-1)]
                                    + x[flattenIndex(0, N-2)]
                                    + x[flattenIndex(0, N-1)]);
      x[flattenIndex(0, 0)]     = 0.33f * (x[flattenIndex(1, 0)]
                                    + x[flattenIndex(0, 1)]
                                    + x[flattenIndex(0, 0)]);

      x[flattenIndex(N-1, 0)]     = 0.33f * (x[flattenIndex(N-2, 0)]
                                    + x[flattenIndex(N-1, 1)]
                                    + x[flattenIndex(N-1, 0)]);
      x[flattenIndex(N-1, N-1)]   = 0.33f * (x[flattenIndex(N-2, N-1)]
                                    + x[flattenIndex(N-1, N-2)]
                                    + x[flattenIndex(N-1, N-1)]);   
}



void project(float[] velocX, float[] velocY, float[] p, float[] div)
{
  for (int j = 1; j < N - 1; j++) 
  {
      for (int i = 1; i < N - 1; i++) 
      {
          div[flattenIndex(i, j)] = -0.5f*(
                   velocX[flattenIndex(i+1, j)]
                  -velocX[flattenIndex(i-1, j)]
                  +velocY[flattenIndex(i, j+1)]
                  -velocY[flattenIndex(i, j-1)]
              )/N;
          p[flattenIndex(i, j)] = 0;
      }
  }
    
    setBounds(0, div); 
    setBounds(0, p);
    
    linearSolve(0, p, div, 1, 6);
    
     for (int j = 1; j < N - 1; j++) 
     {
        for (int i = 1; i < N - 1; i++) 
        {
            velocX[flattenIndex(i, j)] -= 0.5f * (  p[flattenIndex(i+1, j)]
                                            -p[flattenIndex(i-1, j)]) * N;
            velocY[flattenIndex(i, j)] -= 0.5f * (  p[flattenIndex(i, j+1)]
                                            -p[flattenIndex(i, j-1)]) * N;
        }
      }
    
    setBounds(1, velocX);
    setBounds(2, velocY);
}


//generic function
void diffuse(int b, float[] x, float[] x0, float diff, float dt)
{
    float a = dt * diff * (N - 2) * (N - 2);
    linearSolve(b, x, x0, a, 1 + 6 * a);
}


// generic function
void linearSolve(int b, float[] x, float[] x0, float a, float c)
{
  float cRecip = 1.0 / c;
  for (int k = 0; k < iteration; k++) 
  {
    for (int j = 1; j < N - 1; j++) 
    {
        for (int i = 1; i < N - 1; i++) 
        {
          x[flattenIndex(i, j)] = (x0[flattenIndex(i, j)] 
          + a*(x[flattenIndex(i+1, j)]
          +x[flattenIndex(i-1, j)]
          +x[flattenIndex(i, j+1)]
          +x[flattenIndex(i, j-1)]
          +x[flattenIndex(i, j)]
          +x[flattenIndex(i, j)]
           )) * cRecip;
        }
    }
      setBounds(b, x);
  }
}
