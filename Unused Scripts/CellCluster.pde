/*
class CellCluster {
  SimpleCell[] cells;
  int inputs;
  //  float Wx, dWx, inputSum;
  float inputSum;
  float[] Wx, dWx;

  float[]Sh, dSh;
  float[] Wc, dWc;
  float Bc, dBc;
  float Bx, dBx;
  float Sx, Sy;
  float[] dSx;
  boolean active = true;
  float[]S0;
  float threshold = .2f;
  CellCluster(int in, int clusterCells) {
    inputs  = in;

    cells   = new SimpleCell[clusterCells];
    Sh      = new float[clusterCells];
    dSh      = new float[clusterCells];
    dSx     = new float[in];
    Wc      = new float[clusterCells];
    dWc      = new float[clusterCells];
    dWx = new float[in];
    Wx = new float[in];
    S0 = new float[in];
    for (int i=0; i<cells.length; i++) {
      cells[i] = new SimpleCell(in);
    }

    for (int i=0; i<Wx.length; i++) {
      Wx[i] = random(-1,1);
    }
  }


  //------------------------------
  float forward(float[] in) {

    arrayCopy(in, S0);
    //S0=in;
    Sx = 0;
    Sy = 0;
    active = true;
    for (int i=0; i<in.length; i++) {
      Sx += in[i]*Wx[i];
    }
    Sx += Bx;
    //Sx = softSign(Sx);
    if (abs(Sx) <= threshold) {
      active=false;
      Sy = 0;
    } else {
      Sy = 0;
      for (int i=0; i<cells.length; i++) {
        Sy += cells[i].forward(in);
      }

      //Sy+=Bc;
      //Sy = softSign(Sy);
    }
    //print((active ? "X": "O")+" ");
    return Sy;
  }
  //------------------------------
  // returns array of errors
  float[] backward(float dy) {
    dSx = new float[dSx.length]; // NOT EFFFICIENT
    float dW0;


    dW0 = (abs(Sx) <= threshold ? dy:-dy);

    for (int i=0; i<Wx.length; i++)
      dWx[i]+= (S0[i]*dW0);
    dBx += dW0;
    if (!active)
      for (int i=0; i<dSx.length; i++)
        dSx[i] =(Wx[i]*dW0);

    if (active) {
      float[] tmp = new float[dSx.length]; // EVEN LESS EFFFICIENT
      for (int i=0; i<cells.length; i++) {
        dSx = cells[i].backward(dy);
      }
    }
    return dSx;
  }
  //------------------------------

  void update(float lr) {
    float thisLR = lr*0.1f;
    if (active) {
      for (int i=0; i<cells.length; i++) {
        cells[i].update(lr);
        Wc[i] += dWc[i]*thisLR;
      }
      dWc =new float[dWc.length];
    }
    for (int i=0; i<Wx.length; i++) {
      Wx[i] += (dWx[i])*thisLR;
      //if (Wx[i]>=threshold)
      // Wx[i] -= (lr*1e-3);
      //Wx[i] *= (1f-lr*1e-4);



      // print(Wx[i]+" ");
    }
    dWx =new float[dWx.length];

    Bc +=dBc*thisLR;
    dBc = 0;

    Bx += dBx*thisLR;
    dBx = 0;
  }
  //------------------------------
  void resetStates() {
    for (int i=0; i<cells.length; i++) {
      cells[i].resetStates();
    }
  }
  //------------------------------
  float softSign(float x) {
    return 2.0 * x / (1.0 + abs(x));
  }
  //------------------------------
  float deSoftSign(float y) {
    return (2.0 / sq(1.0+abs(y)));
  }
}
*/
